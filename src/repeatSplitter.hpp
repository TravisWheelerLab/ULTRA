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

typedef struct {
  int pos;
  std::string leftConsensus;
  std::string rightConsensus;
} RepeatSplit;

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

  std::string consensus();

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

std::vector<RepeatSplit> *SplitRepeat(RepeatRegion *r,
                              float threshold,
                              int windowUnits,
                              int minSize,
                              int minLagtime);

class RepeatSplitter {
public:
  PairWindow *left;
  PairWindow *right;

  RepeatRegion *repeat;


};

class SplitWindow {
public:
  int start1, start2;
  int end1, end2;

  int minWindowSize;
  int maxWindowSize;

  int minWindowUnits;
  int maxWindowUnits;

  float *window1;
  float *window2;

  float *window1Counts;
  float *window2Counts;


  RepeatRegion *repeat;

  SplitWindow(RepeatRegion *r, int minUnits, int maxUnits, float backgroundCount);
  ~SplitWindow();
  double JSDivergence();
  double FillWindows(int start);
  double ConvolveWindow(bool &notDone);

};



#endif // ULTRA_REPEATSPLITTER_H
