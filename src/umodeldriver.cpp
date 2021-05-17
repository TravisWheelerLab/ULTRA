//
//  umodeldriver.cpp
//  ULTRA
//

#include "umodeldriver.hpp"
#include <cmath>
#include <vector>

void UModelDriver::AnalyzeSequence(SequenceWindow *sequence, int start,
                                   int end) {

  int tstart = start;
  int tend = end;

  if (start >= end)
    return;
  int numberOfPasses = ((end - start) / windowSize);

  if (numberOfPasses * windowSize < (end - start))
    numberOfPasses += 1;

  for (int pass = 0; pass < numberOfPasses; ++pass) {
    ++windowCount;
    if (pass == 0) {
      tend = tstart + matrix->length;

      if (tend > end)
        tend = end;
      int p = 0;

      for (int i = tstart; i < tend; ++i, ++p) {
        model->CalculateCurrentColumn(sequence, i);
      }
    }

    else {
      if (tstart >= end)
        return;
      int p = matrix->length - windowSize;

      for (int i = tstart; i < tend; ++i, ++p) {
        model->CalculateCurrentColumn(sequence, i);
      }
    }

    bool shouldFlush = false;
    if (pass == numberOfPasses - 1)
      shouldFlush = true;

    int length = matrix->length;
    if (tend - start < length)
      length = tend - start;

    int nucStart = tend - matrix->length;
    if (nucStart < 0) {
      nucStart = 0;
    }

    AnalyzeWindow(length, nucStart, shouldFlush);

    if (abs(matrix->previousScoreColumn[0]) > adjustmentThreshold) {
      matrix->AdjustScoreToZero((matrix->length - windowSize) + 2);
    }

    tstart = tend;
    tend = tstart + windowSize;

    if (tend > end)
      tend = end;
  }
}

void UModelDriver::AnalyzeWindow(int length, int start, bool flush) {
  if (repeats == NULL) {
    repeats = new std::vector<str_region>;
  }

  length = length - 1;

  matrix->currentTracebackColumn[0] = 0;
  matrix->ForwardTraceback(currentTraceback, length, 0);
  cell *desc = matrix->cellDescriptions;

  int end = length;
  if (!flush)
    end = end - (matrix->length - windowSize);

  str_region rep;

  int p = length - 1;
  for (int i = 0; i < end; ++i, --p) {
    int c = currentTraceback[i];

    if (c != 0) {
      rep.startScore = matrix->PreviousScore(c, p);
      rep.maxScore = NEG_INF;
      rep.windowStart = i;
      rep.startingPeriod = desc[c].order;
      rep.sequenceStart = i + start - rep.startingPeriod;

      for (i = i + 1, p = p - 1; i < end; ++i, --p) {
        c = currentTraceback[i];

        double s = matrix->PreviousScore(c, p);
        if (s > rep.maxScore)
          rep.maxScore = s;

        if (c == 0) {
          break;
        }
      }

      rep.windowEnd = i;
      rep.arEnd = i;
      rep.sequenceEnd = start + i;
      rep.dScore = rep.maxScore - rep.startScore;
      rep.windowIndex = windowCount;
      rep.lastAdjustment = matrix->lastScoreAdjustment;
      rep.totalAdjustment = matrix->totalScoreAdjustment;

      str_region prev;
      if (repeats->size() > 0)
        prev = repeats->back();

      if (rep.dScore >= threshold && rep.sequenceEnd > prev.sequenceEnd) {

        if (repeats->size() > 0) {
          if (prev.sequenceEnd >= rep.sequenceStart) {

            prev.continues = true;
            repeats->pop_back();
            repeats->push_back(prev);
          }
        }

        repeats->push_back(rep);
      }
    }
  }
}

UModelDriver::UModelDriver(UModel *m, int ws) {
  repeats = NULL;

  model = m;
  matrix = m->matrix;

  windowCount = 0;

  windowSize = ws;
  windowStart = 0;
  windowEnd = 0;

  currentTraceback = (int *)malloc(sizeof(int) * (matrix->length + 2));
}
