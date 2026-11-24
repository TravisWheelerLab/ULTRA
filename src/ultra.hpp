//
//  ultra.hpp
//  ultraP
//

#ifndef ultra_hpp
#define ultra_hpp

#include <algorithm>
#include <pthread.h>
#include <stdio.h>
#include <unordered_map>
#include <vector>

#include "BEDFileWriter.hpp"
#include "FASTAReader.hpp"
#include "FileReader.hpp"
#include "JSONFileWriter.hpp"
#include "RepeatFileWriter.hpp"
#include "RepeatSplitter.hpp"
#include "TabFileWriter.hpp"
#include "cli.hpp"
#include "mask.hpp"
#include "repeat.hpp"
#include "umatrix.hpp"
#include "umodel.hpp"

class Ultra;

// The model and splitter one worker thread analyses windows with.
typedef struct s_uthread {
  int id;
  Ultra *ultra;
  UModel *model;
  SplitWindow *splitter;
} uthread;

class Ultra {
public:
  Settings *settings = NULL;

  FileReader *reader = NULL;
  std::vector<RepeatFileWriter *> writers;
  std::vector<FILE *> outs;

  int numberOfThreads = 1;
  int minReaderSize = 100;

  FILE *settings_out;

  float scoreThreshold = 0.0;

  bool outputRepeatSequence = true;
  bool outputReadID = false;
  bool multithreading = false;
  bool shuffleSequence = false;

  bool storeSequence = false;
  bool storeTraceback = false;
  bool storeScores = false;
  bool storeProfileNumbers = false;

  bool canOutput = false;
  bool firstRepeat = true;
  unsigned long repeatBuffer = 2000;

  unsigned long long total_coverage = 0;
  unsigned long long last_seq_id = 0;
  unsigned long long last_rep_end = 0;

  bool storeTraceAndSequence = false;

  std::unordered_map<unsigned long long, std::vector<mregion> *> masks_for_seq{};
  std::vector<RepeatRegion *> outRepeats{};
  std::vector<UModel *> models{};

  int count = 0;
  int passID = 0;
  std::vector<uthread *> threads{};
  int correctedReadIDs = -1;

  void StoreMaskForRegion(RepeatRegion *r);

  void AnalyzeFile();
  void AnalyzeFileWithThread();
  void AnalyzeFileWithWorkers();
  void FindRepeatsInWindow(SequenceWindow *sequence, uthread *thrd,
                           std::vector<RepeatRegion *> &repeats);
  void BeginWindow(SequenceWindow *sequence);
  void FinishWindow(std::vector<RepeatRegion *> &repeats);

  // Read ids of the window being applied and of the oldest window that may
  // still have repeats to come
  unsigned long long activeReadID = 0;
  unsigned long long smallestReadID = 0;
  void AnalyzeSequenceWindow(SequenceWindow *sequence, uthread *thrd);
  void OutputRepeats(bool flush = false);
  void OutputRepeat(RepeatRegion *r, bool isSubRep = false);

  std::vector<RepeatRegion *> *FindRepeatsInString(const std::string &s);

  void OutputULTRASettings();
  void InitializeWriter();

  void CorrectOverlap(int maxReadID);

  double Log2PvalForScore(float score, float period) const;
  double PvalForScore(float score) const;

  SequenceWindow *GetSequenceWindow(SequenceWindow *seq);

  unsigned long long SmallestReadID();

  void SortRepeatRegions();
  unsigned long long Coverage();


  Ultra(Settings *settings);
  ~Ultra();
};

class CompareRepeatOrder {
public:
  bool operator()(RepeatRegion *lhs, RepeatRegion *rhs);
};


#endif /* ultra_hpp */
