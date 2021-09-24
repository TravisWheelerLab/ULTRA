//
//  repeat.hpp
//  ultraP
//

#ifndef repeat_hpp
#define repeat_hpp

#include <stdio.h>
#include <string>

#include "SequenceWindow.hpp"
#include "Symbol.hpp"
#include "umatrix.hpp"

#define OC_NONE 0
#define OC_OVERLAP_LEFT 1
#define OC_OVERLAP_RIGHT 2
#define OC_OVERLAP_SPLIT 3
#define OC_PERIOD_LEFT 4
#define OC_PERIOD_RIGHT 5
#define OC_PERIOD_SPLIT 6

class RepeatRegion {
public:
  std::string sequenceName;

  std::string traceback;
  std::string sequence;

  unsigned long sequenceID;
  unsigned long readID;

  unsigned long winOverlapSize;
  unsigned long winTotalLength;

  bool combinedRepeat;

  int falseStart;
  unsigned long windowStart;
  unsigned long sequenceStart;

  unsigned long repeatLength;

  int repeatPeriod;

  double startScore;
  double endScore;
  double regionScore;
  double *scores;

  int *logoMemory;
  int **logo;
  symbol *consensus;

  int *logoNumbers;

  int *lookBack;
  int *lookForward;
  char *forwardCounts;
  char *backCounts;

  int mismatches;
  int insertions;
  int deletions;

  int overlapCorrection;

  double logPVal;

  void CreateLogo(SequenceWindow *window, UMatrix *matrix);
  void CreateLogoWithoutMatrix(); // Requires traceback + sequence + lookBack

  void CreateConsensusFromLogo();
  std::string GetConsensus();

  void StoreSequence(SequenceWindow *window);
  void StoreTraceback(UMatrix *matrix);
  void StoreScores(UMatrix *matrix);

  void LookBackDistance();
  void LookForwardDistance();

  void CountFromDistance(int *look, char *counts, const int dir,
                         const int depth);

  void CountForwardBackwardConsensus(const int depth);
  std::vector<RepeatRegion *> *SplitRepeats(const int depth, const int cutoff);
  void GetLogoNumbers();

  RepeatRegion *SubRepeat(int start, int length);

  RepeatRegion();
  ~RepeatRegion();
};

RepeatRegion *GetNextRepeat(SequenceWindow *window, UMatrix *matrix, int *i);

#endif /* repeat_hpp */
