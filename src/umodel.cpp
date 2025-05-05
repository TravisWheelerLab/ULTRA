/*
 umodel.cpp

 */

#include "umodel.hpp"
#include <cmath>

/*** UMODEL SCORE ADJUSTMENT METHODS ***/

void UModel::SetMatchProbabilities(float p) {
  float ip = 1.0 - p;
  ip = ip / 3.0;

  for (int i = 1; i < 5; ++i) {
    for (int j = 1; j < 5; ++j) {
      if (i == j)
        matchProbabilities[i][j] = p;

      else
        matchProbabilities[i][j] = ip;
    }
  }
}

void UModel::SetATCGProbabilities(float A, float T, float C, float G) {

  backgroundProbabilties[N_A] = A;
  backgroundProbabilties[N_T] = T;
  backgroundProbabilties[N_C] = C;
  backgroundProbabilties[N_G] = G;
}

void UModel::CalculateScores() {

  // tscore[CT_BACKGROUND][*]
  tscore[CT_BACKGROUND][CT_BACKGROUND] = log2(1.0 - tp_zeroToMatch);

  float sum = 0.0;
  for (int i = 0; i < matrix->maxPeriod; ++i) {
    sum += pow(periodDecay, i);
  }

  sum = 1.0 / sum;
  periodDecayOffset = log2(sum);
  periodDecay = log2(periodDecay);

  if (adjustForMaximumPeriod) {
    tscore[CT_BACKGROUND][CT_MATCH] =
        log2(tp_zeroToMatch / (float)matrix->maxPeriod);
  }

  else {
    tscore[CT_BACKGROUND][CT_MATCH] = log2(tp_zeroToMatch);
  }

  tscore[CT_BACKGROUND][CT_INSERTION] = NEG_INF;
  tscore[CT_BACKGROUND][CT_DELETION] = NEG_INF;

  // tscore[CT_MATCH][*]
  tscore[CT_MATCH][CT_BACKGROUND] = log2(tp_matchToZero);

  tscore[CT_MATCH][CT_MATCH] =
      log2(1.0 - (tp_matchToZero + tp_matchToInsertion + tp_matchToDeletion));
  tscore[CT_MATCH][CT_INSERTION] = log2(tp_matchToInsertion);
  tscore[CT_MATCH][CT_DELETION] = log2(tp_matchToDeletion);

  // tscore[CT_INSERTION][*]
  tscore[CT_INSERTION][CT_BACKGROUND] = NEG_INF;
  tscore[CT_INSERTION][CT_MATCH] = 0;
  tscore[CT_INSERTION][CT_INSERTION] = log2(tp_consecutiveInsertion);
  tscore[CT_INSERTION][CT_DELETION] = NEG_INF;

  // tscore[CT_DELETION][*]
  tscore[CT_DELETION][CT_BACKGROUND] = NEG_INF;
  tscore[CT_DELETION][CT_MATCH] = 0;
  tscore[CT_DELETION][CT_INSERTION] = NEG_INF;
  tscore[CT_DELETION][CT_DELETION] = log2(tp_consecutiveDeletion);

  for (int i = 1; i < 5; ++i) {
    matchProbabilities[0][i] = backgroundProbabilties[i];
  }

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      logmp[i][j] = log2(matchProbabilities[i][j]);
      mscore[i][j] = logmp[i][j];
      lbp[i][j] = log2(backgroundProbabilties[i] * backgroundProbabilties[j]);
    }

    bscore[i] = log2(backgroundProbabilties[i]);
  }
}

// This does not check if index - d > 0
float UModel::PreviousEmissionScore(SequenceWindow *seq, int index, int order,
                                    int d) {

  int p = index - d;
  float b = bscore[seq->seq[p]];

  float v = mscore[seq->seq[p - order]][seq->seq[p]];

  return v - b;
}

float UModel::EmissionScore(SequenceWindow *seq, int index, int order) {

  symbol s1 = seq->seq[index];
  symbol s2 = seq->seq[index - order];

  return mscore[s2][s1] - bscore[s1];
}

void UModel::CalculateCurrentColumn(SequenceWindow *sequence, int nucIndex,
                                    bool *canBeRepetitive) {
  bool cbr = false; // can be repetitive - true if it is possible for next
                    // character to be repetitve
  float bestRepeatScore = NEG_INF;

  cell *desc = matrix->cellDescriptions;
  float *p = matrix->previousScoreColumn;
  float *c = matrix->currentScoreColumn;

  int *ct = matrix->currentTracebackColumn;

  int matind = matrix->previousColumnIndex;

  // CalculateEmissionScores(sequence->symbolFreqs[nucIndex]);

  for (int row = 0; row < matrix->cellsPerColumn; ++row) {
    int bestToZero = 0;
    switch (desc[row].type) {

      // ****ZEROTH ORDER CELL****
      // Can be transfered to from: zeroth order, nth order
    case CT_BACKGROUND: {
      // Calculate none->none
      c[row] = p[row] + tscore[CT_BACKGROUND][CT_BACKGROUND];
      ct[0] = row;

      break;
    }

      // ****NTH ORDER CELL****
      // Can be transfered to from: zeroth order, nth order, insertion,
      // deletion This will also check to see if it C_NONE inherits from
      // C_MATCH This does not check insertion->match or deletion->match -
      // those are done in
      //      insertion and deletion cells
    case CT_MATCH: {

      // Check to see if match->none is a better score than none->none
      int parentIndex = desc[row].parentIndex;
      int order = desc[row].order;
      float score = p[row] + tscore[CT_MATCH][CT_BACKGROUND];
      score += (order)*periodDecay;
      score += periodDecayOffset;

      if (score > c[parentIndex]) {
        c[parentIndex] = score;
        ct[0] = row;
        cbr = true;
        bestToZero = order;
      }

      // Before we calculate normal match scores we need to check to see
      // if we have had enough characters pass
      if (matind < order) {
        c[row] = NEG_INF;
        ct[desc[row].order] = 0;
        break;
      }

      // Calculate none->match

      // Testing out having transitions to state n come from 0th order n
      // characters previous

      float zscore = 0.0;

      if (immediateTransitionToRepeat) {
        zscore = p[0];
      }

      else {
        zscore = matrix->PreviousScore(0, order - 1);
        zscore += tscore[CT_BACKGROUND][CT_BACKGROUND] * order;
      }

      float es = EmissionScore(sequence, nucIndex, order);

      score = zscore + tscore[CT_BACKGROUND][CT_MATCH] + es;

      c[row] = score;
      ct[desc[row].order] = 0;

      // Calculate match->match
      score = p[row] + tscore[CT_MATCH][CT_MATCH] + es;
      if (score > c[row]) {
        c[row] = score;
        ct[desc[row].order] = row;
      }

      // Check to see if this is the best chance of being a m->0
      // transition
      if (c[row] > bestRepeatScore)
        bestRepeatScore = c[row];

      break;
    }

    case CT_INSERTION: {

      int order = desc[row].order;
      int parentIndex = desc[row].parentIndex;
      int peridocity = order + desc[row].indelNumber;
      int minIndex = peridocity + order + 1; // 2 * order + indelnum

      //  int indelNum = desc[row].indelNumber;

      // We can proceed normally`
      if (matind > minIndex) {

        // Check to see if the match state should transition from here
        float score = p[row] + tscore[CT_INSERTION][CT_MATCH] +
                      EmissionScore(sequence, nucIndex, order);
        if (score > c[parentIndex]) {
          c[parentIndex] = score;
          ct[desc[row].order] = row;
        }

        // Calculate new score based off of previous score
        score = p[row];

        // Update changes in match state score
        score -= matrix->PreviousScore(parentIndex, peridocity + 1);
        score += matrix->PreviousScore(parentIndex, peridocity);

        // Update J values

        score -= PreviousEmissionScore(sequence, nucIndex, peridocity, order);
        score += EmissionScore(sequence, nucIndex, peridocity);

        // Update insertion values
        // There is effectively no score for insertions

        c[row] = score;

      }

      // We have to do a manual calculation of the score
      else if (matind == minIndex) {
        float score = matrix->PreviousScore(
            parentIndex,
            peridocity); // The match score this insertion comes from
        score += 0.0;    // The emission score of the actual insertion = 0
                         // = bscore[n] - bscore[n]

        // Transition cost
        score += tscore[CT_MATCH][CT_INSERTION];
        score += tscore[CT_INSERTION][CT_INSERTION] *
                 (float)(desc[row].indelNumber - 1);
        score += log2(1.0 - tp_matchToZero) * (float)order;

        // We don not need to consider the insertion cost at this point
        // because the insertion state is 0th order, and will have
        // emission score = 1 = bscore/bscore = e^0

        // Calculate emission scores from j states
        for (int j = 0; j < desc[row].order; ++j) {

          score += EmissionScore(sequence, nucIndex - j, peridocity);
        }

        c[row] = score;

      }

      else {
        c[row] = NEG_INF;
      }

    } break;

    case CT_DELETION: {

      int parentIndex = desc[row].parentIndex;
      int order = desc[row].order;
      int indelNum = desc[row].indelNumber;
      int minIndex = 3 * order;

      int delOrder = order - indelNum;
      int EOrder = 2 * order - indelNum;

      // Proceed normally
      if (matind > minIndex) {

        // Check to see if the match state can transition from the del
        // state
        float score = p[row] + tscore[CT_DELETION][CT_MATCH] +
                      EmissionScore(sequence, nucIndex, order);
        if (score > c[parentIndex]) {
          c[parentIndex] = score;
          ct[desc[row].order] = row;

          // printf("del[%i for %i:%i,%i,%i,%i] :%f vs %f to p %f (%f,
          // %f= (%f with %f)\n", row, nucIndex, minIndex,
          // parentIndex, order, indelNum, p[row], score,
          // c[parentIndex], d, f, matrix->PreviousScore(parentIndex,
          // order), matrix->PreviousScore(parentIndex, order + 1));
        }

        score = p[row];
        // if (nucIndex > 3000000)
        // printf("%f\n", score);
        // Update changes in match state score
        score -= matrix->PreviousScore(parentIndex, order + 1);
        score += matrix->PreviousScore(parentIndex, order);

        // Update del values
        score -= PreviousEmissionScore(sequence, nucIndex, delOrder, order);
        score += PreviousEmissionScore(sequence, nucIndex, delOrder, indelNum);

        // Update E values
        score -= PreviousEmissionScore(sequence, nucIndex, EOrder, indelNum);
        score += EmissionScore(sequence, nucIndex, EOrder);

        c[row] = score; // score;

      }

      else if (matind == minIndex) {

        // Get the match state score
        float score = matrix->PreviousScore(parentIndex, order);
        // float oscore = p[row];

        // Calculate transition score
        score += tscore[CT_MATCH][CT_DELETION];
        score += tscore[CT_DELETION][CT_DELETION] * (float)(indelNum - 1);
        score += log2(1.0 - tp_matchToZero) * (float)order;

        // Calculate deletion state emission scores
        for (int j = 0; j < order - indelNum; ++j) {
          int place = nucIndex - order + 1 + j;
          score += EmissionScore(sequence, place, delOrder);
        }

        // Calculate E scores
        for (int j = 0; j < indelNum; ++j) {
          int place = nucIndex - j;
          score += EmissionScore(sequence, place, EOrder);
        }

        c[row] = score;

      }

      else {
        c[row] = NEG_INF;
      }

      break;
    }

    default:
      break;
    }

    //  if (bestToZero > 0)
    //    printf("%i\n", bestToZero);
  }

  matrix->MoveMatrixForward();
  // Right now we don't do anything with CanBeRepetitive
}

void UModel::CalculateCurrentColumnWithoutEmission(bool *canBeRepetitive) {
  bool cbr = false; // can be repetitive - true if it is possible for next
                    // character to be repetitve
  float bestRepeatScore = NEG_INF;

  cell *desc = matrix->cellDescriptions;
  float *p = matrix->previousScoreColumn;
  float *c = matrix->currentScoreColumn;

  int *ct = matrix->currentTracebackColumn;

  int matind = matrix->previousColumnIndex;

  // CalculateEmissionScores(sequence->symbolFreqs[nucIndex]);

  for (int row = 0; row < matrix->cellsPerColumn; ++row) {
    int bestToZero = 0;
    switch (desc[row].type) {

      // ****ZEROTH ORDER CELL****
      // Can be transfered to from: zeroth order, nth order
    case CT_BACKGROUND: {
      // Calculate none->none
      c[row] = p[row] + tscore[CT_BACKGROUND][CT_BACKGROUND];
      ct[0] = row;

      break;
    }

      // ****NTH ORDER CELL****
      // Can be transfered to from: zeroth order, nth order, insertion,
      // deletion This will also check to see if it C_NONE inherits from
      // C_MATCH This does not check insertion->match or deletion->match -
      // those are done in
      //      insertion and deletion cells
    case CT_MATCH: {

      // Check to see if match->none is a better score than none->none
      int parentIndex = desc[row].parentIndex;
      int order = desc[row].order;
      float score = p[row] + tscore[CT_MATCH][CT_BACKGROUND];
      score += (order)*periodDecay;
      score += periodDecayOffset;

      if (score > c[parentIndex]) {
        c[parentIndex] = score;
        ct[0] = row;
        cbr = true;
        bestToZero = order;
      }

      // Before we calculate normal match scores we need to check to see
      // if we have had enough characters pass
      if (matind < order) {
        c[row] = NEG_INF;
        ct[desc[row].order] = 0;
        break;
      }

      // Calculate none->match

      // Testing out having transitions to state n come from 0th order n
      // characters previous

      float zscore = 0.0;

      if (immediateTransitionToRepeat) {
        zscore = p[0];
      }

      else {
        zscore = matrix->PreviousScore(0, order - 1);
        zscore += tscore[CT_BACKGROUND][CT_BACKGROUND] * order;
      }

      float es = 0;

      score = zscore + tscore[CT_BACKGROUND][CT_MATCH] + es;

      c[row] = score;
      ct[desc[row].order] = 0;

      // Calculate match->match
      score = p[row] + tscore[CT_MATCH][CT_MATCH] + es;
      if (score > c[row]) {
        c[row] = score;
        ct[desc[row].order] = row;
      }

      // Check to see if this is the best chance of being a m->0
      // transition
      if (c[row] > bestRepeatScore)
        bestRepeatScore = c[row];

      break;
    }

    case CT_INSERTION: {

      int order = desc[row].order;
      int parentIndex = desc[row].parentIndex;
      int peridocity = order + desc[row].indelNumber;
      int minIndex = peridocity + order + 1; // 2 * order + indelnum

      //  int indelNum = desc[row].indelNumber;

      // We can proceed normally`
      if (matind > minIndex) {

        // Check to see if the match state should transition from here
        float score = p[row] + tscore[CT_INSERTION][CT_MATCH];
        if (score > c[parentIndex]) {
          c[parentIndex] = score;
          ct[desc[row].order] = row;
        }

        // Calculate new score based off of previous score
        score = p[row];

        // Update changes in match state score
        score -= matrix->PreviousScore(parentIndex, peridocity + 1);
        score += matrix->PreviousScore(parentIndex, peridocity);

        // Update insertion values
        // There is effectively no score for insertions

        c[row] = score;

      }

      // We have to do a manual calculation of the score
      else if (matind == minIndex) {
        float score = matrix->PreviousScore(
            parentIndex,
            peridocity); // The match score this insertion comes from
        score += 0.0;    // The emission score of the actual insertion = 0
                         // = bscore[n] - bscore[n]

        // Transition cost
        score += tscore[CT_MATCH][CT_INSERTION];
        score += tscore[CT_INSERTION][CT_INSERTION] *
                 (float)(desc[row].indelNumber - 1);
        score += log2(1.0 - tp_matchToZero) * (float)order;

        // We don not need to consider the insertion cost at this point
        // because the insertion state is 0th order, and will have
        // emission score = 1 = bscore/bscore = e^0

        c[row] = score;

      }

      else {
        c[row] = NEG_INF;
      }

    } break;

    case CT_DELETION: {

      int parentIndex = desc[row].parentIndex;
      int order = desc[row].order;
      int indelNum = desc[row].indelNumber;
      int minIndex = 3 * order;

      int delOrder = order - indelNum;
      int EOrder = 2 * order - indelNum;

      // Proceed normally
      if (matind > minIndex) {

        // Check to see if the match state can transition from the del
        // state
        float score = p[row] + tscore[CT_DELETION][CT_MATCH];
        if (score > c[parentIndex]) {
          c[parentIndex] = score;
          ct[desc[row].order] = row;

          // printf("del[%i for %i:%i,%i,%i,%i] :%f vs %f to p %f (%f,
          // %f= (%f with %f)\n", row, nucIndex, minIndex,
          // parentIndex, order, indelNum, p[row], score,
          // c[parentIndex], d, f, matrix->PreviousScore(parentIndex,
          // order), matrix->PreviousScore(parentIndex, order + 1));
        }

        score = p[row];
        // if (nucIndex > 3000000)
        // printf("%f\n", score);
        // Update changes in match state score
        score -= matrix->PreviousScore(parentIndex, order + 1);
        score += matrix->PreviousScore(parentIndex, order);

        c[row] = score; // score;

      }

      else if (matind == minIndex) {

        // Get the match state score
        float score = matrix->PreviousScore(parentIndex, order);
        // float oscore = p[row];

        // Calculate transition score
        score += tscore[CT_MATCH][CT_DELETION];
        score += tscore[CT_DELETION][CT_DELETION] * (float)(indelNum - 1);
        score += log2(1.0 - tp_matchToZero) * (float)order;

        c[row] = score;

      }

      else {
        c[row] = NEG_INF;
      }

      break;
    }

    default:
      break;
    }
  }

  matrix->MoveMatrixForward();
  // Right now we don't do anything with CanBeRepetitive
}

/*** UMODEL CLASS MANAGEMENT ***/

UModel::UModel(int maxPeriod, int maxInsertions, int maxDeletions,
               int matrixLength) {

  this->matrix =
      new UMatrix(maxPeriod, maxInsertions, maxDeletions, matrixLength);
  CalculateScores();
  // Initialize score distributions here
}

UModel::~UModel() {
  if (matrix) {
    delete matrix;
    matrix = nullptr;
  }
}
