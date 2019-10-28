//
//  JSONReader.hpp
//  ultrax
//
//

#ifndef JSONReader_hpp
#define JSONReader_hpp

#include <stdio.h>
#include "json11.hpp"
#include "FASTAReader.hpp"
#include "JSONPass.hpp"
#include "JSONRepeat.hpp"
#include <vector>
#include "Symbol.hpp"

class JSONReader {
public:
    std::string filePath;
    unsigned long maxSeqLength = 10000;
    unsigned long maxOverlapLength = 100;
    
    int readID = 0;
    int sequenceID = 0;
    
    std::vector<JSONPass *>passes;
    std::vector<JSONRepeat *>repeats;
    
    std::vector<SequenceWindow*>windows;
    std::vector<SequenceWindow*>readyWindows;
    
    
    bool multithread = true;
    pthread_mutex_t readyLock;
    pthread_mutex_t waitingLock;
    
    
    SequenceWindow* GetReadyWindow();
    bool            AddReadyWindow(SequenceWindow* window);
    
    SequenceWindow* GetWaitingWindow();
    bool            AddWaitingWindow(SequenceWindow* window);
    
    bool ReadFile();
    bool CreateWindowsForPasses(const std::vector<int>& passes);
    
    void AddWindowsForRepeat(JSONRepeat *r);
    void CreateWindows();
    
    JSONReader(std::string filePath,
               unsigned long maxSeqLeng=10000,
               unsigned long maxOverlapLeng = 100);
    
};

#endif /* JSONReader_hpp */
