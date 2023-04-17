//
//  FASTAReader.hpp
//  ultraP
//

#ifndef FASTAReader_hpp
#define FASTAReader_hpp

#include <stdio.h>

#include <deque>
#include <fstream>
#include <iostream>
#include <pthread.h>
#include <queue>
#include <sstream>
#include <stack>
#include <string>
#include <unistd.h>

#include "SequenceWindow.hpp"

class FASTAReader {
public:
  std::ifstream file;

  std::string line;
  long long linePlace;

  // change to queues
  std::vector<SequenceWindow *> windows;
  std::vector<SequenceWindow *> readyWindows;
  std::vector<SequenceWindow *> waitingWindows;

  // If readWholeFile is set to true, then maxWindows is ignored
  bool readWholeFile;
  unsigned long maxWindows;
  unsigned long maxSeqLength;
  unsigned long maxOverlapLength;

  double A_pctg;
  double T_pctg;
  double C_pctg;
  double G_pctg;

  symbol *overlapBuffer;          //    = NULL;
  unsigned long overlapLength;    //    = 0;

  std::string sequenceName;       // = "";
  unsigned long sequenceID;       // = 0;
  unsigned long readID;           // = 0; // read id's may not be contiguous
  unsigned long symbolsReadInSeq; // = 0;

  bool doneReadingFile;           //  = false;
  bool isReading;

  bool CopyOverlapBufferFromWindow(SequenceWindow *window,
                                   unsigned long overlapLength);
  bool FillWindows();
  bool ReadWindow(SequenceWindow *window);

  bool multithread = true;
  pthread_mutex_t readyLock;
  pthread_mutex_t waitingLock;

  bool ReadSpecialLine(); // Returns true if it is a new sequence

  SequenceWindow *GetReadyWindow();
  bool AddReadyWindow(SequenceWindow *window);

  SequenceWindow *GetWaitingWindow();
  bool AddWaitingWindow(SequenceWindow *window);

  FASTAReader(std::string filePath, unsigned long maxWindows = 100,
              unsigned long maxSeqLength = 10000,
              unsigned long maxOverlapLength = 100);

  FASTAReader(unsigned long randomWindows, unsigned long maxWindows = 100,
              unsigned long maxSeqLength = 10000,
              unsigned long maxOverlapLength = 100);
  ~FASTAReader();
};

#endif /* FASTAReader_hpp */
