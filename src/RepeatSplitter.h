//
// Created by Olson, Daniel (NIH/NIAID) [E] on 4/25/22.
//

#ifndef REPEATSPLITTER_REPEATSPLITTER_H
#define REPEATSPLITTER_REPEATSPLITTER_H
#include <vector>
#include <string>
#include "repeat.hpp"
#include "Symbol.hpp"
#include <math.h>

typedef struct s_splitWindow {
  float Np;
  float Nq;

  float logNp;
  float logNq;


  // Indexing is: c[symbol][position]
  float **cp = nullptr;
  float **cq = nullptr;
  float **logo = nullptr;
  float *logo_totals = nullptr;

  float *p_mat = nullptr;
  float *q_mat = nullptr;
  float *logo_mat = nullptr;

  float plogp;
  float plogq;
  float qlogq;
  float qlogp;

  int alphabetSize;
  int maxPeriod;
  int period;

  int place = 0;
  int length = 0;


  void AllocateSplitWindow(int alphabetSize, int maxPeriod);
  void StartSplitWindow(int period, float startingScore);
  void DeallocSplitWindow();

  float slow_KLD_pq();
  float slow_KLD_qp();

  float KLD_pq();
  float KLD_qp();

  float AdjustMatrix(bool is_p, int i, int symbol, float delta);

  void FillWindow(RepeatRegion *r, int windowSize);
  void MoveWindowForward(RepeatRegion *r);


  std::vector<float>*RegionScores(RepeatRegion *r,
                                   int windowSize);

  std::vector<int> *SplitsForRegion(RepeatRegion *r,
                                    int window_size,
                                    float threshold);

  void FillLogoForRegion(RepeatRegion *r, int start, int end);
  std::string ConsensusForRegion(RepeatRegion *r,
                                 int start,
                                 int end,
                                 float consensus_threshold);
  std::vector<std::string> *ConsensiForSplit(RepeatRegion *r,
                                             std::vector<int> *splits,
                                             float consensus_threshold);


} SplitWindow;

float ConsensusSimilarity(std::string *c1, std::string *c2, int offset);
bool ShouldJoinConsensus(std::string *c1, std::string *c2, float join_threshold);
void ValidateSplits(std::vector<std::string> *consensi,
                    std::vector<int> *splits,
                    float join_threshold);
#endif //REPEATSPLITTER_REPEATSPLITTER_H
