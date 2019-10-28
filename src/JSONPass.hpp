//
//  JSONPass.hpp
//  ultrax
//
//

#ifndef JSONPass_hpp
#define JSONPass_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include "json11.hpp"

class JSONPass {
public:
    int passID;
    std::string version;
    std::vector<std::pair<std::string, std::string>>parameters;
    
    bool InterpretPass(json11::Json pass);
    
};

#endif /* JSONPass_hpp */
