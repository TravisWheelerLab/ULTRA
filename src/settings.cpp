//
//  settings.cpp
//  ultraP
//
//  Created by Daniel Olson on 7/3/19.
//  Copyright © 2019 Daniel Olson. All rights reserved.
//

#include "settings.hpp"

bool setting_param::operator==(const t_set_param  &p) const {
    return p.name == name;
}

std::string Settings::StringUsage() {
    
    return "ultra <arguments> <input sequence path>\n";
}

std::string Settings::StringVersion() {
    
    return "Ultra 0.99.1\n";
}

std::string Settings::StringHelp() {
    std::string helpString = "------------------\n";
    helpString += StringVersion();
    helpString += StringUsage();
    helpString += "------------------\n";
    helpString += "-flag\tname\tdescription [default value / behavior]\n";
    helpString += "------------------\n";
    
    long long longestName = 0;
    for (int i = 0; i < settings.size(); ++i) {
        if (settings[i] == NULL)
            continue;
        
        if (settings[i]->name.size() > longestName)
            longestName = settings[i]->name.size();
    }
    
    for (int i = 0; i < settings.size(); ++i) {
        if (settings[i] == NULL) {
            helpString += "\n";
            continue;
        }
        
        std::string argument = "-" + settings[i]->flag;
        argument = argument + "\t";
        argument = argument + settings[i]->name;
        argument = argument + ":";
        
        long long end = longestName - settings[i]->name.size();
        argument += "  ";
        for (int i = 0; i < end; ++i) {
            argument += " ";
        }
        
        argument = argument + settings[i]->description;
        
        argument += ". [";
        std::string def = DefaultParam(*settings[i]);
        argument += def.c_str();
        argument += "]\n";
        
        helpString += argument;
    }
    
    return helpString;
}

std::string Settings::DefaultParam(setting_param param) {
    
    std::string defaultValue = "";
 //   printf("-%s = %s", param.flag.c_str(), param.name.c_str());
    if (param == scoreThreshold) {
        defaultValue = std::to_string(v_scoreThreshold);
    }
    
    else if (param == repeatUnits) {
        defaultValue = std::to_string(v_repeatThreshold);
    }
    
    else if (param == lengthThreshold) {
        defaultValue = std::to_string(v_lengthThreshold);
    }
    
    else if (param == ATRichness) {
        defaultValue = std::to_string((v_Apctg + v_Tpctg));
      /*  defaultValue += ", A = ";
        defaultValue += std::to_string(v_Apctg);
        defaultValue += ", T = ";
        defaultValue += std::to_string(v_Tpctg);
        defaultValue += ", C = ";
        defaultValue += std::to_string(v_Cpctg);
        defaultValue += ", G = ";
        defaultValue += std::to_string(v_Gpctg);*/

    }
    
    else if (param == ATCGDistribution) {
        defaultValue += "A=";
        defaultValue += std::to_string(v_Apctg);
        defaultValue += ",T=";
        defaultValue += std::to_string(v_Tpctg);
        defaultValue += ",C=";
        defaultValue += std::to_string(v_Cpctg);
        defaultValue += ",G=";
        defaultValue += std::to_string(v_Gpctg);
    }
    
    else if (param == matchProbability) {
        defaultValue = std::to_string(v_matchProbability);
    }
    
    else if (param == numberOfThreads) {
        defaultValue = std::to_string(v_numberOfThreads);
    }
    
    else if (param == maxPeriod) {
        defaultValue = std::to_string(v_maxPeriod);
    }
    
    else if (param == maxInsertions) {
        defaultValue = std::to_string(v_maxInsertion);
    }
    
    else if (param == maxDeletions) {
        defaultValue = std::to_string(v_maxDeletion);
    }
    
    else if (param == zeroToRepeat) {
        defaultValue = std::to_string(v_zeroToMatch);
    }
    
    else if (param == repeatToZero) {
        defaultValue = std::to_string(v_matchToZero);
    }
    
    else if (param == repeatDecay) {
        defaultValue = std::to_string(v_repeatPeriodDecay);
    }
    
    else if (param == repeatToInsertion) {
        defaultValue = std::to_string(v_matchToInsertion);
    }
    
    else if (param == repeatToDeletion) {
        defaultValue = std::to_string(v_matchToDeletion);
    }
    
    else if (param == consecutiveInsertion) {
        defaultValue = std::to_string(v_consecutiveInsertion);
    }
    
    else if (param == consecutiveDeletion) {
        defaultValue = std::to_string(v_consecutiveDeletion);
    }
    
    
    
    else if (param == outFilePath) {
        defaultValue = "STDOUT";
    }
    
    else if (param == hideRepeatSequence) {
        defaultValue = "Show Repeat Sequence";
    }
    
    else if (param == showWindowID) {
        defaultValue = "Hide Window ID";
    }
    
    else if (param == showTraceback) {
        defaultValue = "Hide traceback";
    }
    
    else if (param == readWholeFile) {
        defaultValue = "Do not read whole file";
    }
    
    else if (param == windowSize) {
        defaultValue = std::to_string(v_windowSize);
    }
    
    else if (param == overlapSize) {
        defaultValue = std::to_string(v_overlapSize);
    }
    
    else if (param == numberOfWindows) {
        defaultValue = std::to_string(v_numberOfWindows);
    }
    
    return defaultValue;
}

bool Settings::CheckFlags(bool printErrors) {
    
    for (int i = 0; i < settings.size(); ++i) {
        if (settings[i] == NULL)
            continue;
        
        for (int j = 0; j < settings.size(); ++j) {
            
            if (i == j || settings[j] == NULL)
                continue;
            
            
            if (settings[i]->flag == settings[j]->flag) {
                if (printErrors) {
                    printf("%i (%s) and %i (%s) share the same flag (%s)",
                           i,
                           settings[i]->name.c_str(),
                           j,
                           settings[j]->name.c_str(),
                           settings[i]->flag.c_str());
                }
                
                return false;
            }
        }
    }
    
    return true;
}

int Settings::InterpretArgument(setting_param   arg,
                                int             argc,
                                const char**    argv,
                                int             index)
{
    index++;
    int itemsLeft = argc - index;
    std::vector<std::string> arguments;
    arguments.reserve(arg.numberOfArguments);
    
    
    if (arg.numberOfArguments > itemsLeft) {
        printf("%s requires %i arguments. See 'ultra - h'.\n", arg.flag.c_str(), arg.numberOfArguments);
        exit(0);
    }
    
    else {
        for (int i = 0; i < arg.numberOfArguments; ++i) {
            int apos = index + i;
            
            if (apos >= argc) {
                printf("Not enough arguments provided for: '%s'\n.", arg.flag.c_str());
                exit(-1);
            }
            std::string arg = argv[apos];
            arguments.push_back(arg);
        }
        
        if (arg == scoreThreshold) {
            v_scoreThreshold = std::stod(arguments[0]);
        }
        
        else if (arg == repeatUnits) {
            v_repeatThreshold = std::stoi(arguments[0]);
        }
        
        else if (arg == lengthThreshold) {
            v_lengthThreshold = std::stoi(arguments[0]);
        }
        
        else if (arg == ATRichness) {
            double val = std::stod(arguments[0]);
            v_Apctg = val / 2.0;
            v_Tpctg = val / 2.0;
            
            val = 1.0 - val;
            
            v_Cpctg = val / 2.0;
            v_Gpctg = val / 2.0;
        }
        
        else if (arg == ATCGDistribution) {
            v_Apctg = std::stod(arguments[0]);
            v_Tpctg = std::stod(arguments[1]);
        
            v_Cpctg = std::stod(arguments[2]);
            v_Gpctg = std::stod(arguments[3]);
        }
        
        
        else if (arg == matchProbability) {
            v_matchProbability = std::stod(arguments[0]);
        }
        
        else if (arg == numberOfThreads) {
            v_numberOfThreads = std::stoi(arguments[0]);
        }
        
        else if (arg == maxPeriod) {
            v_maxPeriod = std::stoi(arguments[0]);
        }
        
        else if (arg == maxInsertions) {
            v_maxInsertion = std::stoi(arguments[0]);
        }
        
        else if (arg == maxDeletions) {
            printf("%s\n", arguments[0].c_str());
            v_maxDeletion = std::stoi(arguments[0]);
        }
        
        else if (arg == zeroToRepeat) {
            v_zeroToMatch = std::stod(arguments[0]);
        }
        
        else if (arg == repeatToZero) {
            v_matchToZero = std::stod(arguments[0]);
        }
        
        else if (arg == repeatDecay) {
            v_repeatPeriodDecay = std::stod(arguments[0]);
        }
        
        else if (arg == repeatToInsertion) {
            v_matchToInsertion = std::stod(arguments[0]);
        }
        
        else if (arg == repeatToDeletion) {
            v_matchToDeletion = std::stod(arguments[0]);
        }
        
        else if (arg == consecutiveInsertion) {
            v_consecutiveInsertion = std::stod(arguments[0]);
        }
        
        else if (arg == consecutiveDeletion) {
            v_consecutiveDeletion = std::stod(arguments[0]);
        }
        
        else if (arg == outFilePath) {
            v_outFilePath = arguments[0];
        }
        
        else if (arg == hideRepeatSequence) {
            v_outputRepeatSequence = false;
        }
        
        else if (arg == showTraceback) {
            v_showTraceback = true;
        }
        
        else if (arg == showWindowID) {
            v_showWindowID = true;
        }
        
        else if (arg == readWholeFile) {
            v_readWholeFile = true;
        }
        
        else if (arg == windowSize) {
            v_windowSize = std::stoi(arguments[0]);
        }
        
        else if (arg == overlapSize) {
            v_overlapSize = std::stoi(arguments[0]);
        }
        
        else if (arg == numberOfWindows) {
            v_numberOfWindows = std::stoi(arguments[0]);
        }
        
        // Do some checks for the help flag and what have you...
    }
    
    return (int)arguments.size();
}


Settings::Settings(int argc, const char * argv[]) {
    
    
    settings.push_back(&outFilePath);
    settings.push_back(&scoreThreshold);
    settings.push_back(&lengthThreshold);
    settings.push_back(&repeatUnits);
    settings.push_back(&ATRichness);
    settings.push_back(&ATCGDistribution);
    settings.push_back(&matchProbability);
    settings.push_back(NULL);
    
    settings.push_back(&numberOfThreads);
    settings.push_back(&maxPeriod);
    settings.push_back(&maxInsertions);
    settings.push_back(&maxDeletions);
    
    settings.push_back(NULL);
    
    settings.push_back(&zeroToRepeat);
    settings.push_back(&repeatToZero);
    settings.push_back(&repeatDecay);
    settings.push_back(&repeatToInsertion);
    settings.push_back(&repeatToDeletion);
    settings.push_back(&consecutiveInsertion);
    settings.push_back(&consecutiveDeletion);
    
    settings.push_back(NULL);
    
    settings.push_back(&hideRepeatSequence);
    settings.push_back(&showTraceback);
    settings.push_back(&showWindowID);
    
    settings.push_back(&readWholeFile);
    settings.push_back(&windowSize);
    settings.push_back(&overlapSize);
    settings.push_back(&numberOfWindows);
    
    if (!CheckFlags()) {
        printf("Flags were found to be inconsistent.\n \
               Check settings.hpp and see if there are duplicate flags being reused.\n");
        exit(0);
    }
    
    if (argc < 2) {
        printf("No arguments recieved. \n");
        printf("%s\n", StringHelp().c_str());
        exit(0);
    }
    
    
    if (argv[argc - 1][0] != '-') {
        v_filePath = argv[argc - 1];
        argc = argc - 1;
    }
    
    
    
    for (int i = 1; i < argc; ++i) {
        // We found a command
        if (argv[i][0] == '-') {
            std:: string argument = argv[i];
            argument = argument.substr(1);
            
            setting_param *param = NULL;
            for (int x = 0; x < settings.size(); ++x) {
                if (settings[x] == NULL)
                    continue;
                
                if (argument == settings[x]->flag) {
                    param = settings[x];
                }
            }
            
            if (param == NULL) {
                if (strcmp(argv[i], "-h") == 0) {
                    printf("%s\n", StringHelp().c_str());
                }
                
                else {
                    printf("Unknown flag '%s'.\n", argv[i]);
                }
                exit(0);
            }
            
            else {
               i += InterpretArgument(*param, argc, argv, i);
            }
            
            
            
        }
        
        else {
            printf("Unrecognized argument %s.\n Commands must start with '-'\n", argv[i]);
            exit(0);
        }
    }
}

Settings::Settings() {
    
}
