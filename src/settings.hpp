//
//  settings.hpp
//  ultraP
//

#ifndef settings_hpp
#define settings_hpp

#include "FileReader.hpp"
#include <cstring>
#include <stdio.h>
#include <string>
#include <vector>

#define ULTRA_VERSION_STRING "0.99.17"

typedef file_type output_type;

typedef struct t_set_param {
  std::string name;
  std::string description;
  std::string flag;
  int numberOfArguments;

  bool operator==(const t_set_param &) const;
} setting_param;

class Settings {
public:
  int o_argc;
  const char **o_argv;
  double v_scoreThreshold = -10000.0;
  int v_repeatThreshold = 3;
  int v_lengthThreshold = 10;

  double v_Apctg = 0.3;
  double v_Tpctg = 0.3;
  double v_Cpctg = 0.2;
  double v_Gpctg = 0.2;

  double v_matchProbability = 0.80;

  double v_zeroToMatch = 0.01;
  double v_matchToZero = 0.05;
  double v_repeatPeriodDecay = 0.85;

  double v_matchToInsertion = 0.02;
  double v_matchToDeletion = 0.02;

  double v_consecutiveInsertion = 0.02;
  double v_consecutiveDeletion = 0.02;

  bool v_calculateLogPVal = true;

  double v_exponLocM = -0.13;
  double v_exponLocB = 2.82;

  double v_exponScaleM = 0.081;
  double v_exponScaleB = 1.28;

  int v_numberOfThreads = 1;
  int v_maxPeriod = 15;
  int v_maxInsertion = 8;
  int v_maxDeletion = 7;

  unsigned long v_randomWindows = 0;

  std::string v_filePath = "";
  std::string v_outFilePath = "";
  bool v_JSONInput = false;
  std::vector<int> v_JSONPasses;
  int v_passID = -1;

  bool v_readWholeFile = false;
  int v_numberOfWindows = 1024;
  int v_windowSize = 8192;
  int v_overlapSize = 100;

  bool v_showLogoNumbers = false;
  bool v_showScores = false;
  bool v_showTraceback = false;

  bool v_showWindowID = false;

  bool v_outputRepeatSequence = true;
  bool v_correctOverlap = false;
  bool v_debugOverlapCorrection = false;

  int v_maxSplitPeriod = 20;
  int v_minSplitWindow = 20;
  int v_splitDepth = 5;
  float v_splitThreshold = 3.6;


  // TODO
  //  Plug in the settings for -json
  //  with BED as default
  output_type v_outputFormat = BED;

  //***************************
  //*****PARAMS/INTERFACE******
  //***************************

  // Model parameters

  setting_param maxSplitPeriod = {
      "Repeat Split Max Period",
      "Maximum repeat period that will be considered for splitting", "sr", 1};

  setting_param splitThreshold = {
      "Repeat Split Threshold",
      "Threshold value used during splitting (larger is more conservative)", "sc",
      1};

  setting_param splitDepth = {
      "Split Depth", "Number of repeat units to consider when splitting repeat",
      "sd", 1};

  setting_param minSplitWindow = {
      "Min Repeat Splitting Window", "Minimum value a repeat splitting window can be",
      "sw", 1};

  setting_param scoreThreshold = {
      "Score Threshold", "Minimum score necessary for a repeat to be recorded",
      "s", 1};

  setting_param repeatUnits = {
      "Repeat Unit Threshold",
      "Minimum number of repeat units necessary for a repeat to be recorded",
      "mu", 1};

  setting_param lengthThreshold = {
      "Region Length Trheshold",
      "Minimum total repeat length necessary for a repeat to be recorded", "ml",
      1};

  setting_param ATRichness = {
      "AT Richness", "Frequency of As and Ts in the input sequence", "at", 1};

  setting_param ATCGDistribution = {
      "ATCG Distribution", "Frequency of As Ts Cs and Gs in the input sequence",
      "atcg", 4};

  setting_param matchProbability = {"Match Probability",
                                    "Probability of two nucleotides in "
                                    "consecutive repeat units being the same",
                                    "m", 1};

  setting_param numberOfThreads = {
      "Number of Threads",
      "Number of threads used by ultra to analyze the input sequence", "n", 1};

  setting_param maxPeriod = {
      "Maximum Period", "Largest repeat period detectable by Ultra", "p", 1};

  setting_param maxInsertions = {
      "Maximum Consecutive Insertions",
      "Maximum number of insertions that can occur in tandem", "mi", 1};

  setting_param maxDeletions = {
      "Maximum Consecutive Deletions",
      "Maximum number of deletions that can occur in tandem", "md", 1};

  setting_param zeroToRepeat = {"Repeat Probability",
                                "Probability of transitioning from the "
                                "nonrepetitive state to a repetitive state",
                                "nr", 1};

  setting_param repeatToZero = {"Repeat Stop Probability",
                                "Probability of transitioning from a "
                                "repetitive state to the nonrepetitive state",
                                "rn", 1};

  setting_param repeatDecay = {"Repeat probability decay",
                               "**reminder to write a good description later",
                               "pd", 1};

  setting_param repeatToInsertion = {"Insertion Probability",
                                     "Probability of transitioning from a "
                                     "repetitive state to an insertion state",
                                     "ri", 1};

  setting_param repeatToDeletion = {"Deletion Probability",
                                    "Probability of transitioning from a "
                                    "repetitive state to a deletion state",
                                    "rd", 1};

  setting_param consecutiveInsertion = {
      "Consecutive Insertion Probability",
      "The probability of transitioning from an insertion state to another "
      "insertion state",
      "ii", 1};

  setting_param consecutiveDeletion = {
      "Consecutive Deletion Probability",
      "The probability of transitioning from a deletion state to another "
      "deletion state",
      "dd", 1};

  setting_param outFilePath = {"Output File Path", "Path to output JSON data",
                               "f", 1};

  setting_param hideRepeatSequence = {
      "Hide Repeat Sequence",
      "Don't show the repetitive sequences in the results JSON", "hs", 0};

  setting_param showWindowID = {
      "Show Window ID",
      "Display the windowID corresponding to a repeat in the results JSON",
      "wid", 0};

  setting_param showScores = {
      "Show score Deltas",
      "Output the score change per residue",
      "ss",
  };

  setting_param showTraceback = {
      "Show traceback", "Output the Viterbi traceback in the results JSON",
      "st", 0};

  setting_param showLogoNumbers = {
      "Show logo numbers",
      "Output the corresponding logo annotation for a given repeat", "sl", 0};

  setting_param JSONInput = {"Read JSON file",
                             "Process all passes in JSON file", "json", 0};

  setting_param JPasses = {"Process passes in JSON file",
                           "Process selected passes in JSON file", "jpass", 1};

  setting_param JSONPassID = {"Pass ID", "Assigns a custom pass ID", "pid", 1};

  setting_param readWholeFile = {
      "Completely Read File",
      "Read the entire input file during initialization.", "R", 0};

  // lots of bugs, we're not using this for now
  setting_param randomSeq = {
      "Randomly generate n windows",
      "Run Ultra on a randomly generated windows instead of a FASTA file",
      "RAN", 1};

  setting_param windowSize = {
      "Window Size",
      "The number of nucleotides per sequence window",
      "ws",
      1,
  };

  setting_param overlapSize = {
      "Window Overlap",
      "The number of nucleotides overlaped between two consecutive windows",
      "os", 1};

  setting_param numberOfWindows = {
      "Number of Windows", "Maximum number of windows stored in memory at once",
      "wn", 1};

  setting_param debugOverlapCorrection = {
      "Debug overlap correction",
      "Report overlap correction in repeat information", "doc", 0};

  setting_param showVersion = {"Ultra Version", "Shows Ultra's version", "v",
                               0};

  setting_param useBED = {"BED output", "Output BED", "bed", 0};

  std::vector<setting_param *> settings;

  bool CheckFlags(bool printErrors = true); // Makes sure flags aren't reused
  int InterpretArgument(setting_param arg, int argc, const char *argv[],
                        int index); // returns the new index
  std::string StringUsage();
  std::string StringVersion();
  std::string StringHelp();
  std::string JSONStringForArgument(setting_param arg);
  std::string JSONString();

  std::string DefaultParam(setting_param param);

  Settings();
  Settings(int argc, const char *argv[]);
};
#endif /* settings_hpp */
