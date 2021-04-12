//
//  main.cpp
//  ultraP
//

#include "FASTAReader.hpp"
#include "settings.hpp"
#include "ultra.hpp"
#include <iostream>
#include <math.h>

#include "../lib/json11.hpp"
#include <fstream>
#include <streambuf>
#include <string>

#include "JSONReader.hpp"

int main(int argc, const char *argv[]) {

  Settings settings = Settings(argc, argv);
  // const char * argu[] ={"./ultra", "/Users/danielolson/Desktop/chr22.fa"};
  // Settings settings = Settings(2, argu);

  Ultra *ultra = new Ultra(&settings, 0);
  ultra->AnalyzeFile();
  ultra->OutputRepeats(true);

  if (ultra->settings->v_outFilePath != "")
    fclose(ultra->out);

  return 0;
}
