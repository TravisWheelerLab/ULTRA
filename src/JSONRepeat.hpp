//
//  JSONRepeat.hpp
//  ultrax
//
//  Created by Daniel Olson on 9/25/19.
//  Copyright © 2019 Daniel Olson. All rights reserved.
//

#ifndef JSONRepeat_hpp
#define JSONRepeat_hpp

#include <stdio.h>

#include <string>
#include <vector>
#include "json11.hpp"

class JSONRepeat {
public:
    
    int             passID;
    
    unsigned long   start;
    unsigned long   length;
    
    int             period;
    
    double          score;
    
    int             substitutions;
    int             insertions;
    int             deletions;
    
    std::string     consensus;
    std::string     sequence;
    std::string     traceback;
    std::string     scores;
    std::string     sequenceName;
    
    int             overlap_code;
    
    
    std::vector<JSONRepeat*> subrepeats;
    
    void OutputRepeat(FILE *out, float threshold, bool first);
    
    bool InterpretRepeat(json11::Json repeat);
    void FixSubrepeatOverlap();
    
};

#endif /* JSONRepeat_hpp */
