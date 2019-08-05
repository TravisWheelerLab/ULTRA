//
//  FASTAReader.cpp
//  ultraP
//


#include "FASTAReader.hpp"

bool FASTAReader::ReadSpecialLine() {
    if (line[0] == '>' || line[0] == ';') {
        sequenceName =  line.substr(1, std::string::npos);
        sequenceID++;
        symbolsReadInSeq = 0;
        overlapLength = 0;
        linePlace = READ_ALL;
        return true;
    }
    
    return false;
}

bool FASTAReader::CopyOverlapBufferFromWindow(SequenceWindow* window, unsigned long long ol) {
    
    if (ol > window->length)
        ol = window->length;
    
    overlapLength = ol;
    
    symbol *ov = &window->newSeq[window->length] - ol;
    
    for (unsigned long long i = 0; i < ol; ++i) {
        // printf("%i vs %i\n", i, ol);
        overlapBuffer[i] = ov[i];
    }
    
    
    return true;
}

bool FASTAReader::FillWindows() {
    
    // printf("************\nFILLING\n************\n");
    if (doneReadingFile) {
        //   printf("*****************************\nDONE.\n*********");
        if (file.is_open())
            file.close();
        return true;
    }
    
    // If readWholeFile is set to true, then this function is only called once....
    if (readWholeFile) {
        
        SequenceWindow *win;
        
        // Read all the sequences
        while (!doneReadingFile || linePlace != READ_ALL) {
            if (waitingWindows.empty()) {
                win = new SequenceWindow(maxSeqLength, maxOverlapLength);
            }
            
            else {
                win = GetWaitingWindow();
            }
            
            if (ReadWindow(win) == false) {
                return false;
            }
            
            AddReadyWindow(win);
            
        }
        
        
    }
    
    else {
        while (!waitingWindows.empty() && (!doneReadingFile || linePlace != READ_ALL)) {
            
            
            SequenceWindow* win = GetWaitingWindow();
            if (ReadWindow(win) == false) {
                return false;
            }
            
            if (win->length > 0)
                AddReadyWindow(win);
            else
                AddWaitingWindow(win);
            
        }
    }
    
    
    return true;
}

bool FASTAReader::ReadWindow(SequenceWindow *window) {
    
    bool resetSymbolCount = false;
    
    
    
    window->PrepareWindow(sequenceName, sequenceID, symbolsReadInSeq, overlapLength);
    
    
    while (true) {
        
        if (linePlace >= 0 && linePlace < READ_ALL) {
            linePlace = window->ReadLine(line, linePlace);
            
            // Check to see if we can't read anything else into the window
            if (linePlace >= 0 && linePlace < READ_ALL) {
                break;
            }
        }
        
        if (linePlace == -1) {
            ReadSpecialLine();
            
            // We're done
            if (window->length > 0) {
                resetSymbolCount = true;
                break;
            }
            
            // We're just beginning
            else {
                symbolsReadInSeq = 0;
                window->PrepareWindow(sequenceName, sequenceID, symbolsReadInSeq, overlapLength);
            }
        }
        
        
        if (linePlace >= READ_ALL) {
            if (std::getline(file, line, '\n')) {
                linePlace = 0;
            }
            
            else {
                doneReadingFile = true;
                file.close();
                break;
            }
        }
    }
    
    symbolsReadInSeq += window->length;
    
    if (resetSymbolCount)
        symbolsReadInSeq = 0;
    
    if (symbolsReadInSeq >= maxOverlapLength)
        CopyOverlapBufferFromWindow(window, maxOverlapLength);
    
    window->readID = readID++;
    window->CopyOverlap(overlapBuffer);
    
    
    
    return true;
}



SequenceWindow *FASTAReader::GetReadyWindow()
{
    
    SequenceWindow *retval = NULL;
    if (multithread)
        pthread_mutex_lock(&readyLock);
    
    if (readyWindows.size() > 0) {
        
        retval = readyWindows.front();
        std::pop_heap(readyWindows.begin(), readyWindows.end(), CompareSequenceWindows());
        readyWindows.pop_back();
        
    }
    
    if (multithread)
        pthread_mutex_unlock(&readyLock);
    return retval;
}

SequenceWindow *FASTAReader::GetWaitingWindow()
{
    SequenceWindow *retval = NULL;
    if (multithread)
        pthread_mutex_lock(&waitingLock);
    
    if (waitingWindows.size() > 0) {
        
        retval = waitingWindows.back();
        waitingWindows.pop_back();
        
    }
    
    if (multithread)
        pthread_mutex_unlock(&waitingLock);
    return retval;
}

bool FASTAReader::AddReadyWindow(SequenceWindow *window)
{
    if (multithread)
        pthread_mutex_lock(&readyLock);
    readyWindows.push_back(window);
    std::push_heap(readyWindows.begin(), readyWindows.end(), CompareSequenceWindows());
    
    if (multithread)
        pthread_mutex_unlock(&readyLock);
    return true;
}

bool FASTAReader::AddWaitingWindow(SequenceWindow *window)
{
    if (multithread)
        pthread_mutex_lock(&waitingLock);
    
    
    waitingWindows.push_back(window);
    
    if (multithread)
        pthread_mutex_unlock(&waitingLock);
    
    return true;
}



/*
 SequenceWindow* FASTAReader::GetWindow() {
 
 if (readyWindows.size() > 0) {
 SequenceWindow* window = readyWindows.top();
 readyWindows.pop();
 
 return window;
 }
 
 return NULL;
 }
 
 void FASTAReader::ReturnWindow(SequenceWindow* window) {
 // Right now we don't do any checks, we just add it to the waitingWindows
 waitingWindows.push_back(window);
 }
 */
FASTAReader::FASTAReader(std::string         filePath,
                         unsigned long long  mxWindows,
                         unsigned long long  mxSeqLength,
                         unsigned long long  mxOverlapLength)
{
    
    overlapLength = 0;
    sequenceName = "";
    sequenceID = 0;
    readID = 0;
    symbolsReadInSeq = 0;
    
    isReading = false;
    
    readWholeFile = false;
    
    if (pthread_mutex_init(&waitingLock, NULL) != 0) {
        printf("Failed to create waiting mutex. Exiting.\n");
        exit(-1);
    }
    
    if (pthread_mutex_init(&readyLock, NULL) != 0) {
        printf("Failed to create ready mutex. Exiting.\n");
        exit(-1);
    }
    
    
    if (mxWindows == -1) {
        readWholeFile = true;
        
    }
    
    else {
        windows.reserve(mxWindows + 1);
        waitingWindows.reserve(mxWindows + 1);
        readyWindows.reserve(mxWindows + 1);
    }
    
    doneReadingFile = false;
    
    maxWindows = mxWindows;
    maxSeqLength = mxSeqLength;
    maxOverlapLength = mxOverlapLength;
    
    overlapBuffer = (symbol*)malloc(sizeof(symbol) * (maxOverlapLength + 1));
    
    
    
    //waitingWindows(mxWindows + 1);
    //readyWindows = std::priority_queue<SequenceWindow*, std::deque<SequenceWindow*>, CompareSequenceWindows>(mxWindows + 1);
    
    
    for (int i = 0; i < maxWindows; ++i) {
        SequenceWindow *newWindow = new SequenceWindow(maxSeqLength, maxOverlapLength);
        newWindow->windowID = i;
        
        AddWaitingWindow(newWindow);
        windows.push_back(newWindow);
    }
    
    file.open(filePath, std::ios::in);
    if (!file.is_open()) {
        printf("Unable to open file: %s\n", filePath.c_str());
        exit(-1);
    }
    
    line.reserve(4086);
    line = "";
    linePlace = 0;
}

FASTAReader::~FASTAReader() {
    if (file.is_open())
        file.close();
    
    free(overlapBuffer);
    for (int i = 0; i < maxWindows; ++i) {
        delete windows[i];
    }
    
}
