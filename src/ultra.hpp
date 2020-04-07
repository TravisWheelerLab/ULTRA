//
//  ultra.hpp
//  ultraP
//


#ifndef ultra_hpp
#define ultra_hpp

#include <stdio.h>
#include <vector>

#include "umodel.hpp"
#include "umatrix.hpp"
#include "FASTAReader.hpp"
#include "JSONReader.hpp"
#include "FileReader.hpp"

#include "repeat.hpp"
#include <algorithm>
#include <stdio.h>
//#include <omp.h>
#include <pthread.h>
#include "settings.hpp"

class Ultra;

typedef struct s_uthread {
    int id;
    pthread_t  p_thread;
    Ultra *ultra;
    UModel *model;
    std::vector<RepeatRegion *>repeats{};
} uthread;

class Ultra {
public:

    Settings*       settings = NULL;
    
    FileReader*     reader = NULL;
    //FASTAReader*    reader  = NULL;
    
    int numberOfThreads = 1;
    int minReaderSize = 100;
    int primaryThread = 0;
    
    FILE *out;
    
    double scoreThreshold = 0.0;
    
    bool outputRepeatSequence   = true;
    bool outputReadID           = false;
    bool multithreading         = false;
    bool AnalyzingJSON          = false;
    //bool outputSequenceSummary  = true;
    
    
    bool canOutput = false;
    bool firstRepeat = true;
    unsigned long repeatBuffer = 2000;
    
    bool storeTraceAndSequence = false;
    
    
    std::vector<RepeatRegion *>outRepeats{};
    std::vector<JSONRepeat *>outJRepeats{};
    
    std::vector<UModel *>models{};
    int count = 0;
    int passID = 0;
    std::vector<uthread *> threads{};
    
    void AnalyzeFile();
    void AnalyzeFileWithThread(void *tid);
    void AnalyzeSequenceWindow(SequenceWindow *sequence, uthread *thrd);
    void OutputRepeats(bool flush = false);
    void OutputRepeat(RepeatRegion *r, bool isSubRep = false);
    void OutputJSONRepeats();
    void OutputJSONKey(std::string key);
    void OutputJSONStart();
    bool FixRepeatOverlap();
    
    double Log2PvalForScore(double score, double period);
    
    SequenceWindow *GetSequenceWindow(SequenceWindow *seq);
    
    void SortRepeatRegions();
    
    pthread_mutex_t outerLock;
    pthread_mutex_t innerLock;
    pthread_mutex_t repeatLock;
    
    Ultra(Settings* settings, int numberOfThreads);
};

class CompareRepeatOrder {
public:
    bool operator() (RepeatRegion *lhs, RepeatRegion *rhs);
};



void *UltraThreadLaunch(void *dat);


#endif /* ultra_hpp */
