/*
 umodel.hpp
 
 
 */
#ifndef umodel_hpp
#define umodel_hpp

#include "umatrix.hpp"
#include "SequenceWindow.hpp"
class UModel {
public:
    
    // Class variables
    UMatrix *matrix;
    
    // (tp = transition probabilities)
    double   tp_zeroToMatch              = 0.01;
    double   tp_matchToZero              = 0.02;
    double   tp_matchToInsertion         = 0.07;
    double   tp_matchToDeletion          = 0.05;
    double   tp_consecutiveInsertion     = 0.1;
    double   tp_consecutiveDeletion      = 0.1;
    
    double   periodDecay                 = 0.9;
    double   periodDecayOffset           = 0.0;
    
    // If true then pt_zeroToMatch = tp_zeroToMatch / matrix->maxPeriod
    bool    adjustForMaximumPeriod      = false;
    
    // Eventually this will need to be expanded for proteins
    double   matchProbabilities[5][5] {
        {1.0, .3, .2, .3, .2},          // [N][*]
        {1.00, 0.61, 0.13, 0.13, 0.13},     // [T][*]
        {1.00, 0.13, 0.61, 0.13, 0.13},     // [C][*]
        {1.00, 0.13, 0.13, 0.61, 0.13},     // [A][*]
        {1.00, 0.13, 0.13, 0.13, 0.61}};    // [G][*]
    
    //                                     N     T     C     A     G
    double   backgroundProbabilties[5]   = {1.00, 0.3, 0.2, 0.3, 0.2};
    
    double windowedWeight = 0.75;
    
    //  double matchScores[256][5][5];
    //double backgroundScores[256][5];
    
    double logmp[5][5];
    double lbp[5][5];
    double mscore[5][5];
    double bscore[5];
    double tscore[4][4];
    
    bool immediateTransitionToRepeat = false;
    bool adjustMatchMismatchMatrix = true;
    
    
    // Methods
    
    //void    CalculateEmissionScores (unsigned char symbolFreq);
    
    void SetMatchProbabilities(double p);
    
    void    CalculateScores();    // Needs to be called before general use
    
    double   PreviousEmissionScore   (SequenceWindow   *seq,
                                      int                index,
                                      int                order,
                                      int                d);
    
    double   EmissionScore           (SequenceWindow   *seq,
                                      int                index,
                                      int                order);
    
    void    CalculateCurrentColumn  (SequenceWindow   *seq,
                                     int               nucIndex,
                                     bool             *canBeRepetitive = NULL);
    
    // Clas management
    
    UModel  (UMatrix    *matrix);
    
    UModel  (int        maxPeriod,
             int        maxInsertions,
             int        maxDeletions,
             int        matrixLength);
    
    ~UModel();
};


#endif /* umodel_hpp */
