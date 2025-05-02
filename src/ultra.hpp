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
#include "TabFileWriter.hpp"
#include "FASTAReader.hpp"
#include "FileReader.hpp"
#include "JSONFileWriter.hpp"
#include "RepeatFileWriter.hpp"
#include "RepeatSplitter.hpp"
#include "cli.hpp"
#include "mask.hpp"
#include "repeat.hpp"
#include "umatrix.hpp"
#include "umodel.hpp"

class Ultra;

typedef struct s_uthread {
  int id;
  pthread_t p_thread;
  Ultra *ultra;
  UModel *model;
  std::vector<RepeatRegion *> repeats{};
  int activeReadID = 0;
  int smallestReadID = 0;
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
  int primaryThread = 0;

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
  int last_read_id = 0;
  unsigned long long last_rep_end = 0;

  bool storeTraceAndSequence = false;


  std::unordered_map<unsigned long long, std::vector<mregion> *>
      masks_for_seq{};
  std::vector<RepeatRegion *> outRepeats{};
  std::vector<UModel *> models{};

  int count = 0;
  int passID = 0;
  std::vector<uthread *> threads{};
  int correctedReadIDs = -1;

  void StoreMaskForRegion(RepeatRegion *r);

  void AnalyzeFile();
  void AnalyzeFileWithThread(void *tid);
  void AnalyzeSequenceWindow(SequenceWindow *sequence, uthread *thrd);
  void OutputRepeats(bool flush = false);
  void OutputRepeat(RepeatRegion *r, bool isSubRep = false);

  std::vector<RepeatRegion *>* GetRepeatsForSequence(const std::string &s);

  void OutputULTRASettings();
  void InitializeWriter();

  void CorrectOverlap(int maxReadID);

  double Log2PvalForScore(float score, float period) const;
  double PvalForScore(float score) const;

  SequenceWindow *GetSequenceWindow(SequenceWindow *seq, uthread *uth);

  int SmallestReadID();

  void SortRepeatRegions();
  unsigned long long Coverage();

  pthread_mutex_t outerLock;
  pthread_mutex_t innerLock;
  pthread_mutex_t repeatLock;

  Ultra(Settings *settings);
  ~Ultra();
};

class CompareRepeatOrder {
public:
  bool operator()(RepeatRegion *lhs, RepeatRegion *rhs);
};

void *UltraThreadLaunch(void *dat);

#endif /* ultra_hpp */
