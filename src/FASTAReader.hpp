//
//  FASTAReader.hpp
//  ultraP
//


#ifndef FASTAReader_hpp
#define FASTAReader_hpp

#include <stdio.h>

#include <queue>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>
#include <stack>
#include <pthread.h>

#include "SequenceWindow.hpp"

class FASTAReader {
public:
    
    std::ifstream       file;
    std::string         line;
    long long  linePlace;
    
    
    
    // yucky
    
    /*std::deque <SequenceWindow*> rwdq;
    std::vector <SequenceWindow*>wwvec;
    std::priority_queue<SequenceWindow*, std::deque<SequenceWindow*>, CompareSequenceWindows>    readyWindows;
    std::vector<SequenceWindow*>             waitingWindows;
    std::vector<SequenceWindow*>            windows;
    */
    
    std::vector<SequenceWindow*>windows;
    std::vector<SequenceWindow*>readyWindows;
    std::vector<SequenceWindow*>waitingWindows;
    
    // If readWholeFile is set to true, then maxWindows is ignored
    bool                readWholeFile;
    unsigned long long  maxWindows;
    unsigned long long  maxSeqLength;
    unsigned long long  maxOverlapLength;
    
    symbol*             overlapBuffer;//    = NULL;
    unsigned long long  overlapLength;//    = 0;
    
    std::string         sequenceName;       // = "";
    unsigned long long  sequenceID;         // = 0;
    unsigned long long  readID;             // = 0; // read id's may not be contiguous
    unsigned long long  symbolsReadInSeq;   // = 0;
    
    bool doneReadingFile;  //  = false;
    bool isReading;
    
    bool CopyOverlapBufferFromWindow(SequenceWindow* window, unsigned long long overlapLength);
    bool FillWindows();
    bool ReadWindow(SequenceWindow* window);
    
    bool multithread = true;
    pthread_mutex_t readyLock;
    pthread_mutex_t waitingLock;
    
    bool ReadSpecialLine(); // Returns true if it is a new sequence
    
    SequenceWindow* GetReadyWindow();
    bool            AddReadyWindow(SequenceWindow* window);
    
    SequenceWindow* GetWaitingWindow();
    bool            AddWaitingWindow(SequenceWindow* window);
    
    
    FASTAReader(std::string         filePath,
                unsigned long long  maxWindows=100,
                unsigned long long  maxSeqLength=10000,
                unsigned long long  maxOverlapLength=100);
    
    ~FASTAReader();
    
    
};

#endif /* FASTAReader_hpp */
