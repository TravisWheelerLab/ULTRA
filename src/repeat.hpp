//
//  repeat.hpp
//  ultraP
//

#ifndef repeat_hpp
#define repeat_hpp

#include <stdio.h>
#include <string>

#include "Symbol.hpp"
#include "SequenceWindow.hpp"
#include "umatrix.hpp"


class RepeatRegion {
public:
    
    std::string sequenceName;
    
    std::string sequence;
    
    unsigned long long sequenceID;
    unsigned long long readID;
    
    unsigned long long  windowStart;
    unsigned long long  sequenceStart;
    
    unsigned long long  repeatLength;
    
    int                 repeatPeriod;
    
    double              startScore;
    double              endScore;
    double              regionScore;
    
    int*                logoMemory;
    int**               logo;
    symbol*             consensus;
    
    int mismatches;
    int insertions;
    int deletions;
    
    
    std::string GetConsensus();
    void CreateLogo(SequenceWindow *window, UMatrix *matrix);
    void CreateConsensusFromLogo();
    void StoreSequence(SequenceWindow *window);
    
    RepeatRegion();
    ~RepeatRegion();
    
};

RepeatRegion *GetNextRepeat(SequenceWindow *window, UMatrix *matrix, int *i);

#endif /* repeat_hpp */
