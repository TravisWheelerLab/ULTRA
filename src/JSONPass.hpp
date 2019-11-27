//
//  JSONPass.hpp
//  ultrax
//
//  Created by Daniel Olson on 10/2/19.
//  Copyright © 2019 Daniel Olson. All rights reserved.
//

#ifndef JSONPass_hpp
#define JSONPass_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include "json11.hpp"


class JSONPass {
public:
    
    int         passID;
    std::string version;
    
    std::vector<std::pair<std::string, std::string>>parameters;
    
    void OutputPass(FILE *out);
    bool InterpretPass(json11::Json pass);
    
};

#endif /* JSONPass_hpp */
