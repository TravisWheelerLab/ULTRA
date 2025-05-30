//
//  FASTAReader.cpp
//  ultraP
//

#include "FASTAReader.hpp"

bool FASTAReader::ReadSpecialLine() {
  if (line[0] == '>') {
    sequenceName = line.substr(1, std::string::npos);
    sequenceID++;
    symbolsReadInSeq = 0;
    overlapLength = 0;
    linePlace = READ_ALL;
    return true;
  }

  return false;
}

bool FASTAReader::CopyOverlapBufferFromWindow(SequenceWindow *window,
                                              unsigned long ol) {

  if (ol > window->length)
    ol = window->length;

  overlapLength = ol;

  symbol *ov = &window->newSeq[window->length] - ol;

  for (unsigned long i = 0; i < ol; ++i) {
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

  // If readWholeFile is set to true, then this function is only called
  // once....
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

      if (!ReadWindow(win)) {
        return false;
      }

      AddReadyWindow(win);
    }

  }

  else {
    while (!waitingWindows.empty() &&
           (!doneReadingFile || linePlace != READ_ALL)) {

      SequenceWindow *win = GetWaitingWindow();
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

/*

void FASTAReader::CreateRandomWindow(SequenceWindow *window) {

    double prob[4] = {A_pctg, T_pctg, C_pctg, G_pctg};
    symbol sym[4] = {N_A, N_T, N_C, N_G};
    for (int i = 0; i < window->maxLength; ++i) {


        double d = (double) (rand() % 10000000);
        d = d / 10000000.0;
        double tot = 0.0;

        symbol s = N_N;

        for (int j = 0; j < 4; ++j) {
            tot += prob[j];

            if (d < tot) {
                s = sym[j];
                break;
            }
        }

        window->seq[i] = s;
    }

    window->length = window->maxLength;

}

 */

bool FASTAReader::ReadWindow(SequenceWindow *window) {

  bool resetSymbolCount = false;

  window->PrepareWindow(sequenceName, sequenceID, symbolsReadInSeq,
                        overlapLength);

  if (overlapLength > 0)
    window->CopyOverlap(overlapBuffer);

  while (true) {
    if (linePlace >= 0 && linePlace < READ_ALL) {
      linePlace = window->ReadLine(line, linePlace, total_seq_length);

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
        window->PrepareWindow(sequenceName, sequenceID, symbolsReadInSeq,
                              overlapLength);
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

  if (shuffle) {
    ShuffleSequenceWindow(window);
  }

  if (symbolsReadInSeq >= maxOverlapLength)
    CopyOverlapBufferFromWindow(window, maxOverlapLength);

  window->readID = readID++;

  return true;
}

SequenceWindow *FASTAReader::GetReadyWindow() {

  SequenceWindow *retval = NULL;
  if (multithread)
    pthread_mutex_lock(&readyLock);

  if (readyWindows.size() > 0) {

    retval = readyWindows.front();
    std::pop_heap(readyWindows.begin(), readyWindows.end(),
                  CompareSequenceWindows());
    readyWindows.pop_back();
  }

  if (multithread)
    pthread_mutex_unlock(&readyLock);
  return retval;
}

SequenceWindow *FASTAReader::GetWaitingWindow() {
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

bool FASTAReader::AddReadyWindow(SequenceWindow *window) {
  if (multithread)
    pthread_mutex_lock(&readyLock);
  readyWindows.push_back(window);
  std::push_heap(readyWindows.begin(), readyWindows.end(),
                 CompareSequenceWindows());

  if (multithread)
    pthread_mutex_unlock(&readyLock);
  return true;
}

bool FASTAReader::AddWaitingWindow(SequenceWindow *window) {
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
FASTAReader::FASTAReader(std::string filePath, int mxWindows,
                         unsigned long mxSeqLength,
                         unsigned long mxOverlapLength) {

  total_seq_length = 0;

  overlapLength = 0;
  sequenceName = "";
  sequenceID = 0;
  readID = 0;
  symbolsReadInSeq = 0;

  isReading = false;

  readWholeFile = false;

  A_pctg = 0.3;
  T_pctg = 0.3;
  C_pctg = 0.2;
  G_pctg = 0.2;

  if (pthread_mutex_init(&waitingLock, NULL) != 0) {
    fprintf(stderr, "Failed to create waiting mutex. Exiting.\n");
    exit(-1);
  }

  if (pthread_mutex_init(&readyLock, NULL) != 0) {
    fprintf(stderr, "Failed to create ready mutex. Exiting.\n");
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

  overlapBuffer = (symbol *)malloc(sizeof(symbol) * (maxOverlapLength + 1));

  // waitingWindows(mxWindows + 1);
  // readyWindows = std::priority_queue<SequenceWindow*,
  // std::deque<SequenceWindow*>, CompareSequenceWindows>(mxWindows + 1);

  for (int i = 0; i < maxWindows; ++i) {
    SequenceWindow *newWindow =
        new SequenceWindow(maxSeqLength, maxOverlapLength);
    newWindow->windowID = i;

    AddWaitingWindow(newWindow);
    windows.push_back(newWindow);
  }

  file.open(filePath, std::ios::in);
  if (!file.is_open()) {
    fprintf(stderr, "Unable to open file: %s\n", filePath.c_str());
    exit(-1);
  }

  line.reserve(4086);
  line = "";
  linePlace = 0;
}

FASTAReader::FASTAReader(unsigned long rn, int maxWindows,
                         unsigned long maxSeqLength,
                         unsigned long maxOverlapLength) {
  doneReadingFile = false;
  symbolsReadInSeq = 0;

  total_seq_length = 0;

  readID = 0;

  A_pctg = 0.3;
  T_pctg = 0.3;
  C_pctg = 0.2;
  G_pctg = 0.2;

  unsigned int seed = (unsigned int)(time(NULL) + getpid());
  sequenceName = std::to_string(seed);
  srand(seed);

  if (pthread_mutex_init(&waitingLock, NULL) != 0) {
    fprintf(stderr, "Failed to create waiting mutex. Exiting.\n");
    exit(-1);
  }

  if (pthread_mutex_init(&readyLock, NULL) != 0) {
    fprintf(stderr, "Failed to create ready mutex. Exiting.\n");
    exit(-1);
  }

  if (maxWindows == -1) {
    readWholeFile = true;

  }

  else {
    windows.reserve(maxWindows + 1);
    waitingWindows.reserve(maxWindows + 1);
    readyWindows.reserve(maxWindows + 1);
  }

  doneReadingFile = false;

  maxWindows = maxWindows;
  maxSeqLength = maxSeqLength;
  maxOverlapLength = 0;

  // waitingWindows(mxWindows + 1);
  // readyWindows = std::priority_queue<SequenceWindow*,
  // std::deque<SequenceWindow*>, CompareSequenceWindows>(mxWindows + 1);

  for (int i = 0; i < maxWindows; ++i) {
    SequenceWindow *newWindow =
        new SequenceWindow(maxSeqLength, maxOverlapLength);
    newWindow->windowID = i;

    AddWaitingWindow(newWindow);
    windows.push_back(newWindow);
  }
}

FASTAReader::~FASTAReader() {
  if (file.is_open())
    file.close();

  free(overlapBuffer);
  for (int i = 0; i < maxWindows; ++i) {
    delete windows[i];
  }
  windows.clear();
  readyWindows.clear();
  waitingWindows.clear();
}
