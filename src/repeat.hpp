//
//  repeat.hpp
//  ultraP
//

#ifndef repeat_hpp
#define repeat_hpp

#include <stdio.h>
#include <string>
#include <vector>

#include "SequenceWindow.hpp"
#include "Symbol.hpp"
#include "umatrix.hpp"

#define OC_NONE 0
#define OC_TRUE 1

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

  float startScore;
  float endScore;
  float regionScore;
  float *scores;

  int *logoMemory;
  int **logo;
  symbol *consensus;

  int *logoNumbers;

  int *lookBack;
  int *lookForward;
  char *forwardCounts;
  char *backCounts;

  std::vector<int> *splits;
  std::vector<std::string> *consensi;

  int mismatches;
  int insertions;
  int deletions;

  int overlapCorrection;

  float logPVal;

  void CreateLogo(SequenceWindow *window, UMatrix *matrix);
  void CreateLogoWithoutMatrix(); // Requires traceback + sequence + lookBack

  void CreateConsensusFromLogo();
  std::string GetConsensus();
  std::string string_consensus = "";

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

  std::string PermutationForString(const std::string &str, int offset);
  int CompareStrPerm(const std::string &str, int perm1, int perm2);
  void SortConsensus(int index);
  void SortConsensi();

  RepeatRegion();
  ~RepeatRegion();
};

RepeatRegion *GetNextRepeat(SequenceWindow *window, UMatrix *matrix, int *i);

bool repeats_overlap(RepeatRegion *r1, RepeatRegion *r2,
                     bool require_same_period = true);

// This does not dealloc r1/r2.
// this also assuems that r2
RepeatRegion *joint_repeat_region(RepeatRegion *r1, RepeatRegion *r2);

#endif /* repeat_hpp */
