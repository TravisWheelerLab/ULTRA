//
//  umodeldriver.hpp
//  ULTRA
//

#ifndef umodeldriver_hpp
#define umodeldriver_hpp

#include "umodel.hpp"
#include <vector>

typedef struct t_str_region {
  // Where in the sequence this occurs
  int sequenceStart = 0;
  int sequenceEnd = 0;

  int startingPeriod = 0;

  double startScore = 0;
  double maxScore = -10000;
  double dScore = 0;

  // Traceback/score info
  int arStart = 0;
  int arEnd = 0;
  int arLength = 0;

  char *traceback = NULL;
  double *scores = NULL;
  std::string tracer = "";
  int windowIndex = 0;
  int windowStart = 0;
  int windowEnd = 0;

  double lastAdjustment;
  double totalAdjustment;

  bool continues = false;

} str_region;

class UModelDriver {
public:
  UModel *model;
  UMatrix *matrix;

  std::vector<str_region> *repeats = NULL;

  int windowCount = 0;

  int windowSize = 9500;
  double threshold = 10;

  int windowStart;
  int windowEnd;

  int *currentTraceback = NULL;

  double adjustmentThreshold = 10000;

  void AnalyzeSequence(SequenceWindow *sequence, int start, int end);
  void AnalyzeWindow(int length, int start, bool flush = false);

  UModelDriver(UModel *m, int ws);
};

#endif /* umodeldriver_hpp */
