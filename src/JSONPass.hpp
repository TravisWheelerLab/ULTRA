//
//  JSONPass.hpp
//  ultrax
//

#ifndef JSONPass_hpp
#define JSONPass_hpp

#include "../lib/json11.hpp"
#include <stdio.h>
#include <string>
#include <vector>

class JSONPass {
public:
  int passID;
  std::string version;

  std::vector<std::pair<std::string, std::string>> parameters;

  void OutputPass(FILE *out);
  bool InterpretPass(json11::Json pass);
};

#endif /* JSONPass_hpp */
