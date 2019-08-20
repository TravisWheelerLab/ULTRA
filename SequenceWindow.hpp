//
//  SequenceWindow.hpp
//  ultraP
//


#ifndef SequenceWindow_hpp
#define SequenceWindow_hpp

#include <stdio.h>
#include <string>
#include "Symbol.hpp"

#define  READ_ALL       (1<<29)



class SequenceWindow {
public:
    unsigned long long  windowID;   // A unique number assigned to this window
    unsigned long long  seqID;      // A unique number representing the whole sequence
    unsigned long long  readID;     // A unique number representing this particular read
    
    std::string         sequenceName;
    
    unsigned long long  start; // where the window starts in the sequence
    unsigned long long  end; // Where the window ends in the sequence
    
    
    // We don't actually need seq (which will always be pointing to overlapSeq)
    // but it makes things a little bit simpler to keep track of stuff
    
    symbol* seqMem;     // The actual sequence memory
    symbol* overlapSeq; // where in seqMem the overlap begins
    symbol* newSeq;     // Where in seqMem the nonoverlap bit of the sequence begins
    symbol* seq;        // The beginning of the entire sequence
    
    unsigned long long  length;  // How many nucleotides are in this window
    unsigned long long  overlap; // The length of sequence shared with this.seqID - 1
    
    unsigned long long  symbolCounts[27];
    double              symbolFreqs[27];
    
    unsigned long long  maxLength;  // The largest length of sequence we can store
    unsigned long long  maxOverlap; // The largest amount of overlap we can store
    
    
    void        PrepareWindow(std::string seqName,
                       unsigned long long seqID,
                       unsigned long long start,
                       unsigned long long overlap);
    
    // ReadLine() returns how much of line was read
    // and will returns -1 if the line is a new sequence
    long long   ReadLine(std::string line,
                           long long place);
    void        CopyOverlap(symbol *b);
    void        CalculateSymbolFrequencies();
    
    
    //bool operator<(const SequenceWindow& rhs);
    
    SequenceWindow(unsigned long long maxLength,
                   unsigned long long maxOverlap);
    ~SequenceWindow();
    
    
};

class CompareSequenceWindows {
public:
    bool operator() (SequenceWindow *lhs, SequenceWindow *rhs);
};


#endif /* SequenceWindow_hpp */
