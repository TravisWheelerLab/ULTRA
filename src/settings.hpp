//
//  settings.hpp
//  ultraP
//
//  Created by Daniel Olson on 7/3/19.
//  Copyright © 2019 Daniel Olson. All rights reserved.
//

#ifndef settings_hpp
#define settings_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <string.h>
#define ULTRA_VERSION_STRING       "1.0b"


typedef struct t_set_param{
    std::string name;
    std::string description;
    std::string flag;
    int numberOfArguments;
    
    bool operator==(const t_set_param  &) const;
} setting_param;



class Settings {
public:
    
    double v_scoreThreshold = 0.0;
    
    double v_Apctg = 0.3;
    double v_Tpctg = 0.3;
    double v_Cpctg = 0.2;
    double v_Gpctg = 0.2;
    
    double v_matchProbability       = 0.75;
    
    double v_zeroToMatch            = 0.01;
    double v_matchToZero            = 0.05;
    double v_repeatPeriodDecay      = 0.9;
    
    double v_matchToInsertion       = 0.07;
    double v_matchToDeletion        = 0.05;
    
    double v_consecutiveInsertion   = 0.1;
    double v_consecutiveDeletion    = 0.1;
    
    int v_numberOfThreads   = 1;
    int v_maxPeriod         = 15;
    int v_maxInsertion      = 8;
    int v_maxDeletion       = 7;
    
    std::string v_filePath      = "";
    std::string v_outFilePath   = "";

    bool v_readWholeFile        = false;
    int  v_numberOfWindows      = 1024;
    int  v_windowSize           = 8192;
    int  v_overlapSize          = 100;
    
    bool v_showTraceback        = false;
    
    
    bool v_showWindowID         = false;
    
    bool v_outputRepeatSequence = true;
    
    
    //***************************
    //*****PARAMS/INTERFACE******
    //***************************
    
    // Model parameters
    
    
    
    setting_param scoreThreshold = {
        "Score Threshold",
        "Minimum score necessary for a repeat to be recorded",
        "s",
        1
    };
    
    setting_param ATRichness = {
        "AT Richness",
        "Frequency of As and Ts in the input sequence",
        "at",
        1
    };
    
    setting_param ATCGDistribution = {
        "ACTG Distribution",
        "Frequency of As Ts Cs and Gs in the input sequence",
        "atcg",
        4
    };
    
    setting_param matchProbability = {
        "Match Probability",
        "Probability of two nucleotides in consecutive repeat units being the same",
        "m",
        1
    };
    
    
    setting_param numberOfThreads = {
        "Number of Threads",
        "Number of threads used by ultra to analyze the input sequence",
        "n",
        1
    };
    
    setting_param maxPeriod = {
        "Maximum Period",
        "Largest repeat period detectable by Ultra",
        "p",
        1
    };
    
    setting_param maxInsertions = {
        "Maximum Consecutive Insertions",
        "Maximum number of insertions that can occur in tandem",
        "mi",
        1
    };
    
    setting_param maxDeletions = {
        "Maximum Consecutive Deletions",
        "Maximum number of deletions that can occur in tandem",
        "md",
        1
    };
    
    setting_param zeroToRepeat = {
        "Repeat Probability",
        "Probability of transitioning from the nonrepetitive state to a repetitive state",
        "nr",
        1
    };
    
    setting_param repeatToZero = {
        "Repeat Stop Probability",
        "Probability of transitioning from a repetitive state to the nonrepetitive state",
        "rn",
        1
    };
    
    setting_param repeatDecay = {
        "Repeat probability decay",
        "**reminder to write a good description later",
        "pd",
        1
    };
    
    setting_param repeatToInsertion = {
        "Insertion Probability",
        "Probability of transitioning from a repetitive state to an insertion state",
        "ri",
        1
    };
    
    setting_param repeatToDeletion = {
        "Deletion Probability",
        "Probability of transitioning from a repetitive state to a deletion state",
        "rd",
        1
    };
    
    setting_param consecutiveInsertion = {
        "Consecutive Insertion Probability",
        "The probability of transitioning from an insertion state to another insertion state",
        "ii",
        1
    };
    
    setting_param consecutiveDeletion = {
        "Consecutive Deletion Probability",
        "The probability of transitioning from a deletion state to another deletion state",
        "dd",
        1
    };
    
    setting_param outFilePath = {
        "Output File Path",
        "Path to output csv",
        "f",
        1
    };
    
    setting_param hideRepeatSequence = {
        "Hide Repeat Sequence",
        "Don't show the repetitive sequences in the results CSV",
        "hs",
        0
    };
    
    setting_param showWindowID = {
        "Show Window ID",
        "Display the windowID corresponding to a repeat in the results CSV",
        "wid",
        0
    };
    
    setting_param showTraceback = {
        "Show traceback",
        "Output the Viterbi traceback in the results CSV",
        "st",
        0
    };
    
    
    setting_param readWholeFile = {
        "Completely Read File",
        "Read the entire input file during initialization.",
        "R",
        0
    };
    
    setting_param windowSize = {
        "Window Size",
        "The number of nucleotides per sequence window",
        "ws",
        1,
    };
    
    setting_param overlapSize = {
        "Window Overlap",
        "The number of nucleotides overlaped between two consecutive windows",
        "os",
        1
    };
    
    setting_param numberOfWindows = {
        "Number of Windows",
        "Maximum number of windows stored in memory at once",
        "wn",
        1
    };
    
    
    std::vector<setting_param *> settings;
    
    bool CheckFlags(bool printErrors=true); // Makes sure flags aren't reused
    int InterpretArgument(setting_param arg, int argc, const char * argv[], int index); // returns the new index
    std::string StringUsage();
    std::string StringVersion();
    std::string StringHelp();
    
    std::string DefaultParam(setting_param param);
    
    Settings();
    Settings(int argc, const char * argv[]);
    
};
#endif /* settings_hpp */
