//
// Created by Daniel Olson on 10/20/21.
//

#ifndef ULTRA_REPEATSPLITTER_H
#define ULTRA_REPEATSPLITTER_H

#include <stdio.h>
#include <string>

#include "SequenceWindow.hpp"
#include "Symbol.hpp"
#include "umatrix.hpp"
#include <cmath>


class RepeatRegion;

class PairWindow {
public:
  int start;
  int end;

  float *logo;

  float sum;
  float plogp;
  float plogq;

  int period;



  PairWindow *pair;

  PairWindow(int period, float pseudo_count = 1.0);

  void fillWithPseudoCount(float pseudo_count);
  void calculateTerms();

  float KLD();
  float slowKLD();
  float splitValue();

  // Updates both us and our partner
  void editCell(int i, int j, float delta);
  // Only updates us
  void selfishEditCell(int i, int j, float delta);
  void dumbEditCell(int i , int j, float delta);
  void fillWindow(RepeatRegion *r, int start, int end);
  void dumbMoveWindowForward(RepeatRegion *r);
  void moveWindowForward(RepeatRegion *r);

  //TODO - Need to include some book keeping
  // for moving the windows and keeping track
  // of where in a string we are as well as what our window
  // size is-
};

std::vector<int> *SplitRepeat(RepeatRegion *r,
                              float threshold,
                              int windowUnits,
                              int minSize,
                              int minLagtime);

std::vector<std::string> *ConsensusForSplits(RepeatRegion *r,
                                             std::vector<int> *splits,
                                             float consensusThreshold=0.6);

// This will set all invalid splits to -1
void FilterSplits(std::vector<int> *splits,
                  std::vector<std::string> *consensus,
                  float threshold=0.1,
                  float wildstar_weight=0.25);





#endif // ULTRA_REPEATSPLITTER_H
