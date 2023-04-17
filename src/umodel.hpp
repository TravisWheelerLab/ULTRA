/*
 umodel.hpp


 */
#ifndef umodel_hpp
#define umodel_hpp

#include "SequenceWindow.hpp"
#include "umatrix.hpp"
class UModel {
public:
  // Class variables
  UMatrix *matrix;

  // (tp = transition probabilities)
  float tp_zeroToMatch = 0.01;
  float tp_matchToZero = 0.02;
  float tp_matchToInsertion = 0.07;
  float tp_matchToDeletion = 0.05;
  float tp_consecutiveInsertion = 0.1;
  float tp_consecutiveDeletion = 0.1;

  float periodDecay = 0.9;
  float periodDecayOffset = 0.0;

  // If true then pt_zeroToMatch = tp_zeroToMatch / matrix->maxPeriod
  bool adjustForMaximumPeriod = false;

  // Eventually this will need to be expanded for proteins
  float matchProbabilities[5][5]{{1.0, .3, .2, .3, .2},           // [N][*]
                                 {1.00, 0.61, 0.13, 0.13, 0.13},  // [T][*]
                                 {1.00, 0.13, 0.61, 0.13, 0.13},  // [C][*]
                                 {1.00, 0.13, 0.13, 0.61, 0.13},  // [A][*]
                                 {1.00, 0.13, 0.13, 0.13, 0.61}}; // [G][*]

  //                                     N     T     C     A     G
  float backgroundProbabilties[5] = {1.00, 0.3, 0.2, 0.3, 0.2};

  float windowedWeight = 0.75;

  //  float matchScores[256][5][5];
  // float backgroundScores[256][5];

  float logmp[5][5];
  float lbp[5][5];
  float mscore[5][5];
  float bscore[5];
  float tscore[4][4];

  bool immediateTransitionToRepeat = false;
  bool adjustMatchMismatchMatrix = true;

  // Methods

  // void    CalculateEmissionScores (unsigned char symbolFreq);

  void SetMatchProbabilities(float p);
  void SetATCGProbabilities(float A, float T, float C, float G);

  void CalculateScores(); // Needs to be called before general use

  float PreviousEmissionScore(SequenceWindow *seq, int index, int order, int d);

  float EmissionScore(SequenceWindow *seq, int index, int order);

  void CalculateCurrentColumn(SequenceWindow *seq, int nucIndex,
                              bool *canBeRepetitive = NULL);

  // Clas management

  UModel(UMatrix *matrix);

  UModel(int maxPeriod, int maxInsertions, int maxDeletions, int matrixLength);

  ~UModel();
};

#endif /* umodel_hpp */
