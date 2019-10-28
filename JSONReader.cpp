//
//  JSONReader.cpp
//  ultrax
//
//

#include "JSONReader.hpp"

void JSONReader::AddWindowsForRepeat(JSONRepeat *r) {
    
    
    SequenceWindow *window = NULL;
    
    int p = 0;
    for (unsigned long i = 0; i < r->length; ++i) {
        if (i % maxSeqLength == 0) {
            if (window != NULL) {
                window->length = p - window->overlap;
                AddReadyWindow(window);
                readyWindows.push_back(window);
                
                
            }
            
            window = new SequenceWindow(maxSeqLength, maxOverlapLength);
            
            window->start = r->start + i;
            
            window->readID = readID++;
            window->seqID = sequenceID;
            
            window->jrepeat = r;
            p=0;
            if (i > maxOverlapLength) {
                for (int x = 0; x < maxOverlapLength; ++x) {
                    window->seq[p++] = SymbolForChar(r->sequence[i - maxOverlapLength + x]);
                }
                window->overlap = p;
            }
        }
        window->seq[p++] = SymbolForChar(r->sequence[i]);
    }
    
    if (p > 0) {
        window->length = p - window->overlap;
        AddReadyWindow(window);
        windows.push_back(window);
    }
}

void JSONReader::CreateWindows() {
    
    for (int i = 0; i < repeats.size(); ++i) {
        if (i > 0) {
            if (repeats[i]->sequenceName != repeats[i - 1]->sequenceName) {
                ++sequenceID;
            }
        }
        
        AddWindowsForRepeat(repeats[i]);
        ++readID;
        
    }
}

bool JSONReader::ReadFile() {
    
    std::ifstream t(filePath);
    std::string contents((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());
    std::string er;
    
    json11::Json repeatFile = json11::Json::parse(contents, er);
    if (er.empty() == false) {
        printf("Error while reading file '%s': %s\n", filePath.c_str(), er.c_str());
        t.close();
        exit(-1);
        return false;
    }
    
    if (!repeatFile.is_object()) {
        printf("JSON File '%s' failed to be read\n", filePath.c_str());
        t.close();
        exit(-1);
        return false;
    }
    
    auto passArray = repeatFile["Passes"];
    if (passArray.type() != json11::Json::ARRAY) {
        printf("JSON File '%s' unable to find passes\n", filePath.c_str());
        t.close();
        exit(-1);
        return false;
    }
    
    for (int i = 0; i < passArray.array_items().size(); ++i) {
        auto item = passArray[i];

        JSONPass *pass = new JSONPass;
        if (!pass->InterpretPass(item)) {
            printf("JSON File '%s' unable to interpret pass %i\n", filePath.c_str(), i);
            t.close();
            exit(-1);
            return false;
        }
        
        passes.push_back(pass);
    }
    
    auto repeatArray = repeatFile["Repeats"];
    
    if (repeatArray.type() != json11::Json::ARRAY) {
        printf("JSON File '%s' unable to find repeats\n", filePath.c_str());
        t.close();
        exit(-1);
        return false;
    }
    
    for (int i = 0; i < repeatArray.array_items().size(); ++i) {
        auto item = repeatArray[i];
        JSONRepeat *repeat = new JSONRepeat;
        
        if (!repeat->InterpretRepeat(item)) {
            printf("JSON File '%s' unable to interpret repeat %i\n", filePath.c_str(), i);
            t.close();
            exit(-1);
            return false;
        }
        
        repeats.push_back(repeat);
    }
    
    // Read all repeats into windows...
    
    CreateWindows();
    
    t.close();
    return true;
}

SequenceWindow* JSONReader::GetReadyWindow() {
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

bool JSONReader::AddReadyWindow(SequenceWindow* window) {
    if (multithread)
        pthread_mutex_lock(&readyLock);
    readyWindows.push_back(window);
    std::push_heap(readyWindows.begin(), readyWindows.end(), CompareSequenceWindows());
    
    if (multithread)
        pthread_mutex_unlock(&readyLock);
    return true;
}

SequenceWindow* JSONReader::GetWaitingWindow() {
    return NULL;
}

bool JSONReader::AddWaitingWindow(SequenceWindow* window) {
  
    return true;
}

JSONReader::JSONReader(std::string fp,
                       unsigned long maxSeqLeng,
                       unsigned long maxOverlapLeng) {
    filePath = fp;
    maxSeqLength = maxSeqLeng;
    maxOverlapLength = maxOverlapLeng;
    
    if (pthread_mutex_init(&waitingLock, NULL) != 0) {
        printf("Failed to create waiting mutex. Exiting.\n");
        exit(-1);
    }
    
    if (pthread_mutex_init(&readyLock, NULL) != 0) {
        printf("Failed to create ready mutex. Exiting.\n");
        exit(-1);
    }
    
    windows.reserve(256);
    readyWindows.reserve(256);
}
