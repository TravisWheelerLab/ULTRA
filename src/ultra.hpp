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
#include "RepeatSplitter.h"
#include "mask.h"
#include "repeat.hpp"
#include "cli.hpp"
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
  SplitWindow *splitter;
} uthread;

class Ultra {
public:
  Settings *settings = NULL;

  FileReader *reader = NULL;
  RepeatFileWriter *writer = NULL;

  int numberOfThreads = 1;
  int minReaderSize = 100;
  int primaryThread = 0;

  FILE *out;
  FILE *settings_out;

  float scoreThreshold = 0.0;

  bool outputRepeatSequence = true;
  bool outputReadID = false;
  bool multithreading = false;
  bool AnalyzingJSON = false;

  bool storeSequence = false;
  bool storeTraceback = false;
  bool storeScores = false;
  bool storeProfileNumbers = false;

  bool canOutput = false;
  bool firstRepeat = true;
  unsigned long repeatBuffer = 2000;

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
  void AnalyzeFileWithThread(void *tid);
  void AnalyzeSequenceWindow(SequenceWindow *sequence, uthread *thrd);
  void OutputRepeats(bool flush = false);
  void OutputRepeat(RepeatRegion *r, bool isSubRep = false);

  void OutputULTRASettings();
  void InitializeWriter();

  void CorrectOverlap(int maxReadID);

  float Log2PvalForScore(float score, float period) const;

  SequenceWindow *GetSequenceWindow(SequenceWindow *seq);

  int SmallestReadID();

  void SortRepeatRegions();
  void OutputMaskedFASTA(std::string in_file_path, FILE *out_file);
  pthread_mutex_t outerLock;
  pthread_mutex_t innerLock;
  pthread_mutex_t repeatLock;

  Ultra(Settings *settings, int numberOfThreads);
};

class CompareRepeatOrder {
public:
  bool operator()(RepeatRegion *lhs, RepeatRegion *rhs);
};

void *UltraThreadLaunch(void *dat);

#endif /* ultra_hpp */
