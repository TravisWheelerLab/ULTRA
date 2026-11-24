/*
 umodel.hpp


 */
#ifndef umodel_hpp
#define umodel_hpp

#include "SequenceWindow.hpp"
#include "umatrix.hpp"
#include <stdint.h>
#include <vector>

// Cells per block in the steady column's indel loops.
#define STEADY_LANES 4

// Floats per window position that the blocked kernel keeps (see
// PositionEmissions).
#define POSITION_EMISSION_ARRAYS 10ULL

// Whether a model with these indel counts fills its steady columns with the
// blocked kernel, which also keeps POSITION_EMISSION_ARRAYS floats per window
// position. With fewer than two full blocks per order the fused scalar kernel
// is faster.
inline bool UsesBlockedKernel(unsigned long long maxInsertions,
                              unsigned long long maxDeletions) {
  return maxInsertions / STEADY_LANES + maxDeletions / STEADY_LANES >= 2;
}

// Emission scores by sequence position for the steady-state column. Every
// emission the column scores is diff[seq[a]][seq[b]] (diff[x][y] =
// mscore[x][y] - bscore[y]) with one of the two symbols fixed for a whole run
// of the other, so with ten arrays, Against(y)[a] = diff[seq[a]][y] and
// From(x)[a] = diff[x][seq[a]], each run is a contiguous slice.
//
// A position's values depend only on its symbol and on diff, so the cache
// keeps the symbol each position was filled from. Before every column, Fill
// checks the positions the column will read against the sequence as it is
// now and refills from the first one that differs, then fills up to the
// current position (never ahead of it). The cached values are therefore
// always those the general path would compute from the live sequence. Every
// value is refilled when any bit of mscore or bscore changes.
class PositionEmissions {
public:
  size_t filled = 0; // Positions 0..filled hold values

  const float *Against(symbol y) const { return &values[(size_t)y * stride]; }
  const float *From(symbol x) const {
    return &values[(5 + (size_t)x) * stride];
  }

  // Makes positions upto - back .. upto of seq (seqLen symbols) hold the
  // values for their current symbols and these mscore and bscore; back is
  // how far back the column reads.
  void Fill(const symbol *seq, size_t seqLen, size_t upto, size_t back,
            const float mscore[5][5], const float bscore[5]);

private:
  uint32_t mscoreBits[5][5] = {};
  uint32_t bscoreBits[5] = {};
  float diff[5][5] = {};
  std::vector<symbol> symbols;
  size_t stride = 0;
  std::vector<float> values;
};
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
  void CalculateCurrentColumnWithoutEmission(bool *canBeRepetitive = NULL);

  // Fill steady-state columns with the fast kernels (the result is the same
  // either way; false is for testing)
  bool fastColumns = true;

private:
  PositionEmissions positionEmissions;
  // Per-order scratch for the steady-state column
  std::vector<float> orderBest;
  std::vector<int> orderBestFrom;
  std::vector<float> orderEs;

  void CalculateGeneralRows(SequenceWindow *seq, int nucIndex, int startRow);
  void SteadyColumn(const symbol *seq, size_t i, size_t nOrders);
  void CalculateSteadyColumnFused(const symbol *seq, size_t i,
                                  size_t nOrders);
  void CalculateSteadyColumn(const symbol *seq, size_t i, size_t nOrders);
  bool UseBlockedKernel() const {
    return UsesBlockedKernel(matrix->maxInsertions, matrix->maxDeletions);
  }
  size_t SteadyOrders(size_t seqLen, int nucIndex) const;
  bool ColumnIsSteady(size_t seqLen, int nucIndex) const;

public:

  // Clas management
  UModel(int maxPeriod, int maxInsertions, int maxDeletions, int matrixLength);

  ~UModel();
};

#endif /* umodel_hpp */
