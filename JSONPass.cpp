//
//  JSONPass.cpp
//  ultrax
//
//

#include "JSONPass.hpp"

bool JSONPass::InterpretPass(json11::Json pass) {
    
    if (!pass.is_object()) {
        printf("Not an object\n");
        return false;
    }
    
    passID = pass["Pass ID"].int_value();
    version = pass["Version"].string_value();
    auto params = pass["Parameters"].object_items();
    std::map<std::string, json11::Json>::iterator it = params.begin();
    
    while (it != params.end()) {
        std::string key = it->first;
        std::string value = it->second.string_value();
        
        parameters.push_back(std::pair<std::string, std::string>(key, value));
        
        it++;
    }
    
    return true;
}
