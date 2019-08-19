//
//  repeat.cpp
//  ultraP
//

#include "repeat.hpp"




void RepeatRegion::CreateLogo(SequenceWindow *window, UMatrix *matrix)
{
    
    
    logoMemory = (int*)malloc(sizeof(int) * repeatPeriod * (NUM_SYMBOLS + 1));
    logo = (int**)malloc(sizeof(int *) * repeatPeriod);
    
    for (int i = 0; i < repeatPeriod; ++i) {
        int p = i * NUM_SYMBOLS;
        logo[i] = &logoMemory[p];
        for (int j = 0; j < NUM_SYMBOLS; ++j) {
            logo[i][j] = 0;
        }
        
    }
    
    
    
    unsigned long long backset = repeatPeriod;
    unsigned long long endpoint = repeatLength - repeatPeriod;
    
    unsigned long long fplace = repeatPeriod + windowStart;
    unsigned long long bplace = windowStart;
    
    //printf("%llu %llu %llu\n", backset, fplace, endpoint);
    if (endpoint > window->length + window->overlap)
        endpoint = window->length + window->overlap;
    
    
    
    //printf("\n%i vs %i\n", windowStart, endpoint);
    
    int insertionReset = -1;
    int deletionReset = -1;
    
    int p = 0;
    
    double lastUnitScore = 0;
    int modp = 0;
    
    
    for (unsigned long long i = 0; i < endpoint; ++i) {
        
        int frow = matrix->traceback[fplace];
        int brow = matrix->traceback[bplace];
     //   printf("%f\n", matrix->scoreColumns[fplace][frow]);
       // printf("%i(%i,%i,%i).", row, insertions, deletions, mismatches);
        
        cell fdesc = matrix->cellDescriptions[frow];
        cell bdesc = matrix->cellDescriptions[brow];
        
        if (i < repeatPeriod)
            bdesc = matrix->cellDescriptions[0];
        
        if (fdesc.type == CT_INSERTION) {
            insertions += fdesc.indelNumber;
            i += fdesc.indelNumber;
            fplace += fdesc.indelNumber;
            for (int j = 0; j < fdesc.indelNumber; ++j) {
                
            }
            
           // printf("Insertion.\n");
        }
        
        else if (fdesc.type == CT_DELETION) {
            deletions += fdesc.indelNumber;
            bplace += fdesc.indelNumber;
            p += fdesc.indelNumber;
            modp += fdesc.indelNumber + 1;
            //printf("Deletion.\n");
        }
        
        else {
            ++modp;
        }
        
        
        modp = modp % repeatPeriod;
        
        if (bdesc.type == CT_INSERTION) {
            bplace += bdesc.indelNumber;
        }
        
        else if (bdesc.type == CT_DELETION) {
            fplace += bdesc.indelNumber;
            i += bdesc.indelNumber;
        }
        
        
        

        
        symbol s = window->seq[fplace];
        symbol sb = window->seq[bplace];
      //  printf("%i (%i): %c %c\n", fplace, bplace, CharForSymbol(s), CharForSymbol(sb));
        p = p % repeatPeriod;
        
        if (s != sb) {
            mismatches += 1;
            logo[p][s] += 1;
        }
        
        else {
           // printf("\n**%i %i**\n", p, s);
            logo[p][s] += 4;
            
        }
        
        
        
        
        ++p;
        ++fplace;
        ++bplace;
        --insertionReset;
        --deletionReset;
        
        if (deletionReset == 0) {
            backset = backset + repeatPeriod;
        }
        
        if (insertionReset == 0) {
            backset = repeatPeriod;
        }
        
    }
  /*
    unsigned long long p1 = windowStart + repeatLength;
    unsigned long long p2 = windowStart + repeatLength - modp;
    int t1 = matrix->traceback[p1];
    int t2 = matrix->traceback[p2];
    
    lastUnitScore = matrix->scoreColumns[p2][t2] - matrix->scoreColumns[p1][t1];
    printf("%f - %f = %f with rs = %f\n", matrix->scoreColumns[p2][t2], matrix->scoreColumns[p1][t1], lastUnitScore, regionScore);
    printf("A: %llu %llu %f %f\n", modp, repeatLength, lastUnitScore, regionScore);
    repeatLength -= modp;
    regionScore -= lastUnitScore;
    printf("B: %llu %llu %f %f\n\n", modp, repeatLength, lastUnitScore, regionScore);*/
    
}

void RepeatRegion::CreateConsensusFromLogo()
{
    consensus = (symbol*)malloc(sizeof(symbol) * repeatPeriod);
    
    //printf("%i %i\n", repeatLength, repeatPeriod);
    for (int i = 0; i < repeatPeriod; ++i) {
        symbol s = 0;
        
        for (int j = 1; j < NUM_SYMBOLS; ++j) {
            if (logo[i][j] > logo[i][s])
                s = j;
         //   printf("%i ", logo[i][j]);
        }
      //  printf("\n");
        consensus[i] = s;
    }

}

void RepeatRegion::StoreSequence(SequenceWindow *window) {
    sequence = "";
    sequence.reserve(repeatLength + 1);
    
    for (int i = 0; i < repeatLength; ++i) {
        sequence.push_back(CharForSymbol(window->seq[i + windowStart]));
    }
    sequence[repeatLength] = '\0';
}

RepeatRegion::RepeatRegion() {
    sequenceName = "";
    sequence = "";
    sequenceID = 0;
    windowStart = 0;
    sequenceStart = 0;
    repeatLength = 0;
    repeatPeriod = 0;
    startScore = 0;
    endScore = 0;
    regionScore = 0;
    readID = 0;
    falseStart = 0;
    
    mismatches = 0;
    insertions = 0;
    deletions = 0;
    
    winOverlapSize = 0;
    
    logoMemory = NULL;
    logo = NULL;
    consensus = NULL;
    
    combinedRepeat = false;
}

RepeatRegion::~RepeatRegion() {
    if (logo != NULL) {
        free(logo);
        logo = NULL;
    }
    
    if (consensus != NULL) {
        free(consensus);
        consensus = NULL;
    }
    
    if (logoMemory != NULL) {
        free(logoMemory);
        logoMemory = NULL;
    }
}

std::string RepeatRegion::GetConsensus() {
    std::string con = "";
    con.reserve(repeatPeriod + 1);
    for (int i = 0; i < repeatPeriod; ++i) {
        con.push_back(CharForSymbol(consensus[i]));
    }
    con[repeatPeriod] = '\0';
    
    return con;
}

void RepeatRegion::StoreTraceback(UMatrix *matrix) {
    for (int i = 0; i < repeatPeriod; ++i) {
        traceback.push_back('.');
    }
    
    unsigned long long end = windowStart + repeatLength;
    for (unsigned long long i = windowStart + repeatPeriod; i < end; ++i) {
        cell c = matrix->cellDescriptions[matrix->traceback[i]];
        
        if (c.type == CT_INSERTION) {
            i += c.indelNumber - 1;
            for (int j = 0; j < c.indelNumber; ++j) {
                traceback.push_back('I');
            }
        }
        
        else if (c.type == CT_DELETION) {
            i += c.indelNumber - 1;
            for (int j = 0; j < c.indelNumber; ++j) {
                traceback.push_back('D');
            }
        }
        
        
        else  {
            traceback.push_back('.');
        }
    }
}

// This assumes that matrix has already calculated the traceback
RepeatRegion *GetNextRepeat(SequenceWindow *window, UMatrix *matrix, int *pos)
{
    if (*pos == 0)
        *pos = 1;
    
    bool foundRepeat = false;
    int i;
    
    unsigned long long seqLength = window->length + window->overlap;
    
    RepeatRegion *region = NULL;
    
    for (i = *pos; i < seqLength; ++i) {
        
       // printf("%i: %llx, %llx, %llx\n", i, (unsigned long long)matrix, (unsigned long long)matrix->traceback, (unsigned long long)&matrix->traceback[i]);
        if (matrix->traceback[i] > 0) {
            foundRepeat = true;
            break;
        }
        
        else if (matrix->traceback[i] < 0) {
       //     printf("%i %i\n", i, matrix->traceback[i]);
        }
    }
    
    
   
    
    if (foundRepeat) {
        int row = matrix->traceback[i];
       // printf("%i %i\n", matrix->traceback[i], matrix->traceback[i-1]);
        cell desc = matrix->cellDescriptions[row];
        
        region = new RepeatRegion();
        region->winOverlapSize = window->overlap;
        region->winTotalLength = window->length + window->overlap;
        region->sequenceName = window->sequenceName;
        region->sequenceID = window->seqID;
        region->readID = window->readID;
        region->repeatPeriod = desc.order;
        
        //printf("%i %i %i %i\n", window->start, i, window->overlap, region->repeatPeriod);
        region->windowStart = i - region->repeatPeriod - 1;
        region->sequenceStart = window->start + (i - window->overlap - region->repeatPeriod) - 1;

        region->startScore = matrix->scoreColumns[i][row];
        
        int length = region->repeatPeriod;
        
        for (; i < seqLength; ++i) {
            if (matrix->traceback[i] != 0) {
                length++;
            }
            
            else {
                ++i; // always make i end + 1
                break;
            }
        }
        
        region->repeatLength = length;
        
        row = matrix->traceback[i - 1];
        
        region->endScore = matrix->scoreColumns[i - 1][row];
        region->regionScore = region->endScore - region->startScore;
    }
    
    *pos = i;
    
    if (region != NULL) {
        region->readID = window->readID;
        region->CreateLogo(window, matrix);
        region->CreateConsensusFromLogo();
    }
    
    return region;
}
