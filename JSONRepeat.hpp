//
//  JSONRepeat.hpp
//  ultrax
//
//

#ifndef JSONRepeat_hpp
#define JSONRepeat_hpp

#include <stdio.h>

#include <string>
#include <vector>
#include "json11.hpp"

class JSONRepeat {
public:
    
    int passID;
    
    unsigned long start;
    unsigned long length;
    
    int period;
    
    double score;
    
    int substitutions;
    int insertions;
    int deletions;
    
    std::string consensus;
    std::string sequence;
    std::string traceback;
    
    std::string sequenceName;
    
    int overlap_code;
    
    std::vector<JSONRepeat *> subrepeats;
    
    void OutputRepeat(FILE *out, bool first=true);
    
    bool InterpretRepeat(json11::Json repeat);
};

#endif /* JSONRepeat_hpp */
