//
//  JSONPass.cpp
//  ultrax
//
//  Created by Daniel Olson on 10/2/19.
//  Copyright © 2019 Daniel Olson. All rights reserved.
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

void JSONPass::OutputPass(FILE *out) {
    fprintf(out, "{");
    fprintf(out, "\"Pass ID\": %i,\n", passID);
    fprintf(out, "\"Version\": \"%s\",\n", version.c_str());
    fprintf(out, "\"Parameters\": {\n");
    
    for (int i = 0; i < parameters.size(); ++i) {
        if (i > 0)
            fprintf(out, ",\n");
        fprintf(out, "\"%s\": \"%s\"", parameters[i].first.c_str(), parameters[i].second.c_str());
    }
    
    fprintf(out, "\n}");
    
}
