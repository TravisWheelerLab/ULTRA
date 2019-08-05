//
//  SequenceWindow.cpp
//  ultraP
//


#include "SequenceWindow.hpp"


/*
bool compareSequenceWindows(SequenceWindow *lhs,
                            SequenceWindow *rhs)
{
    return lhs->readID < rhs->readID;
}*/

void SequenceWindow::PrepareWindow(std::string seqName,
                                   unsigned long long sid,
                                   unsigned long long strt,
                                   unsigned long long overlapLength)
{
    
    overlap = overlapLength;
    
    sequenceName = seqName;
    
    start = strt;
    end = start;
    seqID = sid;
    
    length  = 0;
    newSeq = &seqMem[overlap];
    
    
    for (int i = 0; i < 27; ++i) {
        symbolFreqs[i] = 0;
    }
}

long long SequenceWindow::ReadLine(std::string line,
                                   long long place)
{
    //printf ("(mem: %llx overlap:%llx newseq: %llx, length: %llu, place %llu) Reading line: %s\n", (unsigned long long)seqMem, (unsigned long long)overlapSeq, (unsigned long long)newSeq, length, place,  line.c_str());
    unsigned long long windowRemaining = maxLength - length;
    unsigned long long charactersToRead = windowRemaining;
    
    if (charactersToRead > line.length() - place) {
        charactersToRead = line.length() - place;
    }
    
    if (line[0] == ';')
        return READ_ALL;
    
    else if (line[0] == '>')
        return -1;
    
    
    int charactersRead = 0;
    for (int i = 0; i < charactersToRead; ++i) {
        char c = line[place++];
        if (isalnum(c)) {
            charactersRead++;
            symbol s = SymbolForChar(c);
            newSeq[length++] = s;
            symbolCounts[s] += 1;
            
        }
    }
    
  //  printf("%lli vs %lli\n", maxLength, length);
    
    
    end += charactersRead;
    
    if (place >= line.length())
        return READ_ALL;
    
  /*  for (int i = 0; i < length; ++i) {
        printf("%c", CharForSymbol(newSeq[i]));
        if (i % 70 == 0 && i > 0)
            printf("\n");
    }
    printf("\n");*/
    
    return place;
}

void SequenceWindow::CopyOverlap(symbol *b)
{
    
    //printf("%i %i\n", readID, overlap);
    
    for (unsigned long long i = 0; i < overlap; ++i) {
        overlapSeq[i] = b[i];
        
    }
}


void SequenceWindow::CalculateSymbolFrequencies() {
    unsigned long long totalSymbols = 0;
    
    for (int i = 0; i < 27; ++i) {
        totalSymbols += symbolCounts[i];
    }
    
    // Check to see if there are any symbols at all..
    if (totalSymbols == 0)
        return;
    
    for (int i = 0; i < 27; ++i) {
        symbolFreqs[i] = (double)symbolCounts[i] / (double)totalSymbols;
    }
}

/*
bool SequenceWindow::operator<(const SequenceWindow& rhs) {
    if (this->seqID != rhs.seqID)
        return this->seqID < rhs.readID;
    return (this->readID < rhs.readID);
}*/

SequenceWindow::SequenceWindow (unsigned long long ml,
                                unsigned long long mo) {
    
    maxLength = ml;
    maxOverlap = mo;
    
    
    sequenceName.reserve(4096);
    sequenceName = "";
    
    
    seqMem      = (symbol* )malloc(sizeof(symbol) * (maxLength + maxOverlap + 1));
   // unsigned long long memloc = (unsigned long long)seqMem;
   // unsigned long long size = (unsigned long long)(sizeof(symbol) * (maxLength + maxOverlap + 1));
    //printf("Sequence Window (%llx) has created %llx to %llx\n", (unsigned long long)this, (unsigned long long)seqMem, size + memloc);
    overlapSeq  = &seqMem[0];
    newSeq      = &seqMem[mo];
    seq         = &seqMem[0];
    
    length      = 0;
    overlap     = 0;
    start       = 0;
    end         = 0;
    
    windowID    = -1;
    seqID       = -1;
    readID      = -1;
    
    for (int i = 0; i < 27; ++i) {
        symbolFreqs[i] = 0;
        symbolCounts[i] = 0;
    }
    
}

SequenceWindow::~SequenceWindow() {
    //printf("%llx has is being deleted\n", (unsigned long long)this);
    if (seqMem != NULL) {
        free(seqMem);
        seqMem = NULL;
    }
}


bool CompareSequenceWindows::operator() (SequenceWindow *lhs, SequenceWindow *rhs) {
   // printf("comparing %lli vs %lli\n", lhs->readID, rhs->readID);
    return lhs->readID > rhs->readID;
}
