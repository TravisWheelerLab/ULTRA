//
//  main.cpp
//  ultraP
//


#include <iostream>
#include "FASTAReader.hpp"
#include "ultra.hpp"
#include "settings.hpp"
#include <math.h>

#include "json11.hpp"
#include <string>
#include <fstream>
#include <streambuf>

#include "JSONReader.hpp"

int main(int argc, const char * argv[]) {
    
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
