/*
 umodel.cpp

 */

#include "umodel.hpp"
#include <algorithm>
#include <cstring>
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

void PositionEmissions::Fill(const symbol *seq, size_t seqLen, size_t upto,
                             size_t back, const float mscore[5][5],
                             const float bscore[5]) {
  // Any bit of mscore or bscore changed (this runs for every column)
  uint32_t changed = 0;
  for (int x = 0; x < 5; ++x) {
    for (int y = 0; y < 5; ++y) {
      uint32_t b;
      memcpy(&b, &mscore[x][y], sizeof(b));
      changed |= b ^ mscoreBits[x][y];
    }
    uint32_t b;
    memcpy(&b, &bscore[x], sizeof(b));
    changed |= b ^ bscoreBits[x];
  }
  if (changed != 0 || values.empty()) {
    for (int x = 0; x < 5; ++x) {
      for (int y = 0; y < 5; ++y) {
        memcpy(&mscoreBits[x][y], &mscore[x][y], sizeof(uint32_t));
        // diff[x][y] = mscore[x][y] - bscore[y]: the emission score of
        // symbol y against x, the same float as the general path
        diff[x][y] = mscore[x][y] - bscore[y];
      }
      memcpy(&bscoreBits[x], &bscore[x], sizeof(uint32_t));
    }
    filled = 0;
  }
  if (seqLen > stride) {
    stride = seqLen;
    values.assign(10 * stride, 0.0f);
    symbols.assign(stride, 0);
    filled = 0;
  }
  // Positions the column reads that are already filled: keep those up to the
  // first whose symbol has changed
  size_t lo = upto >= back ? upto - back : 0;
  size_t checkedTo = std::min(filled, upto + 1);
  for (size_t a = lo; a < checkedTo; ++a) {
    if (symbols[a] != seq[a]) {
      filled = a;
      break;
    }
  }
  if (filled > upto)
    return;
  size_t from = filled;
  size_t to = upto + 1;
  for (size_t a = from; a < to; ++a)
    symbols[a] = seq[a];
  for (int y = 0; y < 5; ++y) {
    float *against = &values[(size_t)y * stride];
    for (size_t a = from; a < to; ++a)
      against[a] = diff[seq[a]][y];
  }
  for (int x = 0; x < 5; ++x) {
    float *by = &values[(size_t)(5 + x) * stride];
    for (size_t a = from; a < to; ++a)
      by[a] = diff[x][seq[a]];
  }
  filled = to;
}

namespace {

// Updates best / bestFrom with the exits of consecutive rows firstRow.. as the
// general path does, one at a time, keeping the first strictly greater. The
// block's maximum only decides whether any exit needs that check, so the
// result is the same.
inline void TakeExits(const float *exit, int n, float &best, int &bestFrom,
                      size_t firstRow) {
  float m = exit[0];
  for (int l = 1; l < n; ++l)
    m = std::max(m, exit[l]);
  if (m > best) {
    for (int l = 0; l < n; ++l) {
      if (exit[l] > best) {
        best = exit[l];
        bestFrom = (int)(firstRow + l);
      }
    }
  }
}

// One order's insertion and deletion cells of a steady column. Every emission
// of these cells is a contiguous run of one of the position arrays, and every
// look-back a contiguous run of the order's history, so the indel loops run on
// blocks of STEADY_LANES cells, which the compiler turns into vector code.
// Within a cell the operations and their order are those of the general path.
struct SteadyOrder {
  const symbol *seq;
  size_t i;
  size_t o;
  size_t r; // The order's match row
  size_t cur;
  size_t nIns;
  size_t nDel;
  float es; // The order's match emission, EmissionScore(i, o)
  float tInMa;
  float tDeMa;
  const float *hO; // This order's history row
  const PositionEmissions *pe;

  // Writes the indel rows of the current column c from the previous column
  // p, and takes their exits into bestIn / bestFromIn. Kept out of line so
  // that the loop over orders, which most orders pass through with few or no
  // indel cells, stays small.
#if defined(__GNUC__)
  __attribute__((noinline))
#endif
  void
  Fill(const float *p, float *c, float &bestIn, int &bestFromIn) const {
    const int L = STEADY_LANES;
    const float *emitI = pe->Against(seq[i]);
    float best = bestIn;
    int bestFrom = bestFromIn;
    // PreviousEmissionScore(i, q, d) = diff[seq[i - d - q]][seq[i - d]]
    const float *atO = pe->Against(seq[i - o]);
    float parentO = hO[cur - o];
    float parentO1 = hO[cur - o - 1];

    // Fewer cells than a block, and the cells after the last block, are
    // filled one at a time straight from the arrays.
    if (nIns < (size_t)L) {
      for (size_t k = 1; k <= nIns; ++k) {
        size_t row = r + k;
        float pr = p[row];
        float exit = pr + tInMa + es;
        if (exit > best) {
          best = exit;
          bestFrom = (int)row;
        }
        float v = pr;
        v -= hO[cur - o - k - 1];
        v += hO[cur - o - k];
        v -= atO[i - 2 * o - k];
        v += emitI[i - o - k];
        c[row] = v;
      }
    } else {
      size_t n = nIns;
      // Index n - k (k = 1..n) of each array: insertion k's match score
      // o + k + 1 and o + k columns back, and its emissions, i - 2o - k
      // against i - o and i - o - k against i.
      const float *further = &hO[cur - o - 1 - n];
      const float *nearer = &hO[cur - o - n];
      const float *eBack = &atO[i - 2 * o - n];
      const float *eHere = &emitI[i - o - n];
      const float *pin = &p[r + 1];
      float *cin = &c[r + 1];
      size_t j = 0;
      while (j + L <= n) {
        // Cells j..j + L; the arrays hold them at n - L - j.. in reverse.
        size_t q = n - L - j;
        float v[STEADY_LANES], exit[STEADY_LANES];
        for (int l = 0; l < L; ++l) {
          int x = L - 1 - l;
          float t = pin[j + l];
          t -= further[q + x];
          t += nearer[q + x];
          t -= eBack[q + x];
          t += eHere[q + x];
          v[l] = t;
          exit[l] = pin[j + l] + tInMa + es;
        }
        for (int l = 0; l < L; ++l)
          cin[j + l] = v[l];
        TakeExits(exit, L, best, bestFrom, r + 1 + j);
        j += L;
      }
      for (; j < n; ++j) {
        size_t k = j + 1;
        float pr = pin[j];
        float v = pr;
        v -= hO[cur - o - k - 1];
        v += hO[cur - o - k];
        v -= atO[i - 2 * o - k];
        v += emitI[i - o - k];
        cin[j] = v;
        float exit = pr + tInMa + es;
        TakeExits(&exit, 1, best, bestFrom, r + 1 + j);
      }
    }

    if (nDel > 0 && nDel < (size_t)L) {
      const float *fromO = pe->From(seq[i - o]);
      const float *from2O = pe->From(seq[i - 2 * o]);
      for (size_t k = 1; k <= nDel; ++k) {
        size_t row = r + nIns + k;
        float pr = p[row];
        float exit = pr + tDeMa + es;
        if (exit > best) {
          best = exit;
          bestFrom = (int)row;
        }
        float v = pr;
        v -= parentO1;
        v += parentO;
        v -= atO[i - 2 * o + k];
        v += fromO[i - k];
        v -= from2O[i - k];
        v += emitI[i - 2 * o + k];
        c[row] = v;
      }
    } else if (nDel > 0) {
      size_t n = nDel;
      // Deletion k's emissions: i - 2o + k against i - o and i (index
      // k - 1), and i - k against i - o and i - 2o (index n - k).
      const float *eAt = &atO[i - 2 * o + 1];
      const float *eUp = &emitI[i - 2 * o + 1];
      const float *eByO = &pe->From(seq[i - o])[i - n];
      const float *eBy2O = &pe->From(seq[i - 2 * o])[i - n];
      size_t first = r + nIns + 1;
      const float *pdel = &p[first];
      float *cdel = &c[first];
      size_t j = 0;
      while (j + L <= n) {
        size_t q = n - L - j;
        float v[STEADY_LANES], exit[STEADY_LANES];
        for (int l = 0; l < L; ++l) {
          int x = L - 1 - l;
          float t = pdel[j + l];
          t -= parentO1;
          t += parentO;
          t -= eAt[j + l];
          t += eByO[q + x];
          t -= eBy2O[q + x];
          t += eUp[j + l];
          v[l] = t;
          exit[l] = pdel[j + l] + tDeMa + es;
        }
        for (int l = 0; l < L; ++l)
          cdel[j + l] = v[l];
        TakeExits(exit, L, best, bestFrom, first + j);
        j += L;
      }
      for (; j < n; ++j) {
        size_t q = n - 1 - j;
        float pr = pdel[j];
        float t = pr;
        t -= parentO1;
        t += parentO;
        t -= eAt[j];
        t += eByO[q];
        t -= eBy2O[q];
        t += eUp[j];
        cdel[j] = t;
        float exit = pr + tDeMa + es;
        TakeExits(&exit, 1, best, bestFrom, first + j);
      }
    }

    bestIn = best;
    bestFromIn = bestFrom;
  }
};

} // namespace

// Fills the current column for symbol nucIndex of sequence from the previous
// column, then moves the matrix forward one column. Columns past their
// warm-up (see SteadyState) are filled by the steady-state kernels, which do
// the same floating-point operations in the same order as the general path,
// so the result is identical.
void UModel::CalculateCurrentColumn(SequenceWindow *sequence, int nucIndex,
                                    bool *canBeRepetitive) {
  const symbol *seq = sequence->seq;
  size_t seqLen = sequence->maxLength + sequence->maxOverlap + 1;
  size_t maxOrder = matrix->maxPeriod > 0 ? matrix->maxPeriod - 1 : 0;
  if (fastColumns && maxOrder > 0 && nucIndex >= 0 &&
      (size_t)nucIndex < seqLen) {
    // Only the blocked kernel reads emissions by position
    if (UseBlockedKernel()) {
      size_t back = matrix->steady.maxSequenceBack > 0
                        ? (size_t)matrix->steady.maxSequenceBack
                        : 0;
      positionEmissions.Fill(seq, seqLen, (size_t)nucIndex, back, mscore,
                             bscore);
    }
    if (ColumnIsSteady(seqLen, nucIndex)) {
      SteadyColumn(seq, (size_t)nucIndex, maxOrder);
      matrix->MoveMatrixForward();
      return;
    }
    // A column still in its warm-up: the low orders may already be steady.
    // Fill them the fast way and the rest cell by cell, in row order, as the
    // general path would.
    size_t n = SteadyOrders(seqLen, nucIndex);
    if (n > 0) {
      SteadyColumn(seq, (size_t)nucIndex, n);
      CalculateGeneralRows(sequence, nucIndex, matrix->matchRows[n + 1]);
      matrix->RecordHistory(matrix->currentColumnIndex);
      matrix->MoveMatrixForward();
      return;
    }
  }
  CalculateGeneralRows(sequence, nucIndex, 0);
  matrix->RecordHistory(matrix->currentColumnIndex);
  matrix->MoveMatrixForward();
}

// How many orders, from order 1 up, have every cell past its warm-up in the
// current column (see SteadyState::byOrder). The bounds grow with the order,
// so these orders form a prefix.
size_t UModel::SteadyOrders(size_t seqLen, int nucIndex) const {
  long long prev = matrix->previousColumnIndex;
  long long cur = matrix->currentColumnIndex;
  long long length = matrix->length;
  long long i = nucIndex;
  if (i >= (long long)seqLen)
    return 0;
  const std::vector<OrderBounds> &b = matrix->steady.byOrder;
  size_t n = 0;
  for (size_t o = 1; o < b.size(); ++o) {
    if (prev > b[o].minPreviousColumn && cur >= b[o].maxLookback &&
        b[o].maxLookback < length && i >= b[o].maxSequenceBack)
      n = o;
    else
      break;
  }
  return n;
}

// Whether every cell of the current column is past its warm-up: no
// previous-column branch applies, every look-back column is at or before the
// current one without wrapping, and every sequence index is in range.
bool UModel::ColumnIsSteady(size_t seqLen, int nucIndex) const {
  const SteadyState &st = matrix->steady;
  return (long long)matrix->previousColumnIndex > st.minPreviousColumn &&
         (long long)matrix->currentColumnIndex >= st.maxLookback &&
         (unsigned long long)st.maxLookback < matrix->length &&
         (long long)nucIndex >= st.maxSequenceBack &&
         (long long)nucIndex < (long long)seqLen;
}

// The steady-state column for orders 1..nOrders: the blocked kernel, or the
// fused scalar one for models with few indel cells per order.
void UModel::SteadyColumn(const symbol *seq, size_t i, size_t nOrders) {
  if (UseBlockedKernel())
    CalculateSteadyColumn(seq, i, nOrders);
  else
    CalculateSteadyColumnFused(seq, i, nOrders);
}

// CalculateSteadyColumn for models with few indel cells per order: one pass
// over the orders, scalar, which is faster when the orders have fewer than two
// full blocks. The same operations in the same order, so the same result.
void UModel::CalculateSteadyColumnFused(const symbol *seq, size_t i,
                                        size_t nOrders) {
  const float tBgBg = tscore[CT_BACKGROUND][CT_BACKGROUND];
  const float tMaBg = tscore[CT_MATCH][CT_BACKGROUND];
  const float tBgMa = tscore[CT_BACKGROUND][CT_MATCH];
  const float tMaMa = tscore[CT_MATCH][CT_MATCH];
  const float tInMa = tscore[CT_INSERTION][CT_MATCH];
  const float tDeMa = tscore[CT_DELETION][CT_MATCH];

  UMatrix *m = matrix;
  const size_t cur = m->currentColumnIndex;
  const size_t maxIns = m->maxInsertions;
  const size_t maxDel = m->maxDeletions;
  int *ct = m->currentTracebackColumn;
  float *c = m->currentScoreColumn;
  const float *p = m->previousScoreColumn;
  const size_t stride = m->historyStride;
  float *history = m->history;
  // The background (order 0) or order-order match score of the column d >= 1
  // columns back
#define BACK(order, d) history[(order)*stride + cur - (d)]

  // diff[x][y] = mscore[x][y] - bscore[y]: the emission score of symbol y
  // against x, the same float as the general path computes
  float diff[5][5];
  for (int x = 0; x < 5; ++x)
    for (int y = 0; y < 5; ++y)
      diff[x][y] = mscore[x][y] - bscore[y];
  const symbol sI = seq[i];
  // EmissionScore(i, q) and PreviousEmissionScore(i, q, d)
#define EMIT(q) diff[seq[i - (q)]][sI]
#define EMIT_BACK(q, d) diff[seq[i - (d) - (q)]][seq[i - (d)]]

  c[0] = p[0] + tBgBg;
  ct[0] = 0;

  size_t r = 1;
  for (size_t o = 1; o <= nOrders; ++o) {
    int order = (int)o;

    // The match row: first its exit to the background row
    float score = p[r] + tMaBg;
    score += (order)*periodDecay;
    score += periodDecayOffset;
    if (score > c[0]) {
      c[0] = score;
      ct[0] = (int)r;
    }

    float zscore;
    if (immediateTransitionToRepeat) {
      zscore = p[0];
    } else {
      // Order 1 reads the current column's background row live
      zscore = o == 1 ? c[0] : BACK(0, o - 1);
      zscore += tBgBg * order;
    }
    float es = EMIT(o);
    float best = zscore + tBgMa + es;
    int bestFrom = 0;
    float stay = p[r] + tMaMa + es;
    if (stay > best) {
      best = stay;
      bestFrom = (int)r;
    }

    size_t nIns = std::min(maxIns, o - 1);
    size_t nDel = std::min(maxDel, o - 1);

    // This order's match score o and o + 1 columns back. Insertion k reads it
    // o + k and o + k + 1 columns back, so each iteration loads one new value
    // and reuses the previous one.
    float parentO = BACK(o, o);
    float parentO1 = BACK(o, o + 1);
    float nearer = parentO1;
    for (size_t k = 1; k <= nIns; ++k) {
      size_t row = r + k;
      size_t per = o + k;
      float pr = p[row];
      float exit = pr + tInMa + es;
      if (exit > best) {
        best = exit;
        bestFrom = (int)row;
      }
      float further = BACK(o, per + 1);
      float v = pr;
      v -= further;
      v += nearer;
      v -= EMIT_BACK(per, o);
      v += EMIT(per);
      c[row] = v;
      nearer = further;
    }

    for (size_t k = 1; k <= nDel; ++k) {
      size_t row = r + nIns + k;
      size_t delOrder = o - k;
      size_t eOrder = 2 * o - k;
      float pr = p[row];
      float exit = pr + tDeMa + es;
      if (exit > best) {
        best = exit;
        bestFrom = (int)row;
      }
      float v = pr;
      v -= parentO1;
      v += parentO;
      v -= EMIT_BACK(delOrder, o);
      v += EMIT_BACK(delOrder, k);
      v -= EMIT_BACK(eOrder, k);
      v += EMIT(eOrder);
      c[row] = v;
    }

    c[r] = best;
    history[o * stride + cur] = best;
    ct[o] = bestFrom;
    r += 1 + nIns + nDel;
  }
  history[cur] = c[0];
#undef BACK
#undef EMIT
#undef EMIT_BACK
}

// CalculateCurrentColumn for a steady column. The rows are walked in layout
// order (match row of order o, then its insertions, then its deletions)
// without looking up each cell's description, and the look-back columns are
// addressed directly. Every cell does the same floating-point operations in
// the same order as the general path, and the max into a parent row keeps the
// first strictly greater score, so the result is identical.
//
// Fills the background row and orders 1..nOrders, which must all be steady.
// With every order this is the whole column; with fewer, the caller fills the
// remaining rows with CalculateGeneralRows and records the history. Does not
// move the matrix forward.
void UModel::CalculateSteadyColumn(const symbol *seq, size_t i,
                                   size_t nOrders) {
  const PositionEmissions &pe = positionEmissions;
  const float tBgBg = tscore[CT_BACKGROUND][CT_BACKGROUND];
  const float tMaBg = tscore[CT_MATCH][CT_BACKGROUND];
  const float tBgMa = tscore[CT_BACKGROUND][CT_MATCH];
  const float tMaMa = tscore[CT_MATCH][CT_MATCH];
  const float tInMa = tscore[CT_INSERTION][CT_MATCH];
  const float tDeMa = tscore[CT_DELETION][CT_MATCH];

  UMatrix *m = matrix;
  const size_t cur = m->currentColumnIndex;
  const size_t maxIns = m->maxInsertions;
  const size_t maxDel = m->maxDeletions;
  int *ct = m->currentTracebackColumn;
  float *c = m->currentScoreColumn;
  const float *p = m->previousScoreColumn;
  const size_t stride = m->historyStride;
  float *history = m->history;
  const int *matchRows = m->matchRows;

  // EmissionScore(i, q) = diff[seq[i - q]][seq[i]]
  const float *emitI = pe.Against(seq[i]);

  // The background row's running max stays in locals and is written once at
  // the end, so the loop over orders does not go through memory for it.
  float bg = p[0] + tBgBg;
  int bgFrom = 0;

  if (orderBest.size() < nOrders + 1) {
    orderBest.resize(nOrders + 1);
    orderBestFrom.resize(nOrders + 1);
    orderEs.resize(nOrders + 1);
  }
  float *best = orderBest.data();
  int *bestFrom = orderBestFrom.data();
  float *es = orderEs.data();

  // Pass 1: every order's match row, from the background row and from
  // itself, and its exit to the background row, in order.
  for (size_t o = 1; o <= nOrders; ++o) {
    size_t r = matchRows[o];
    int order = (int)o;
    float score = p[r] + tMaBg;
    score += (order)*periodDecay;
    score += periodDecayOffset;
    if (score > bg) {
      bg = score;
      bgFrom = (int)r;
    }

    float zscore;
    if (immediateTransitionToRepeat) {
      zscore = p[0];
    } else {
      // Order 1 reads the current column's background row live
      zscore = o == 1 ? bg : history[cur - (o - 1)];
      zscore += tBgBg * order;
    }
    float e = emitI[i - o];
    float b = zscore + tBgMa + e;
    int from = 0;
    float stay = p[r] + tMaMa + e;
    if (stay > b) {
      b = stay;
      from = (int)r;
    }
    best[o] = b;
    bestFrom[o] = from;
    es[o] = e;
  }

  // Pass 2: the insertion and deletion rows, whose exits can still raise
  // their match row (after its own scores, as in the general path). Order 1
  // has none.
  if (maxIns + maxDel > 0) {
    for (size_t o = 2; o <= nOrders; ++o) {
      size_t nIns = std::min(maxIns, o - 1);
      size_t nDel = std::min(maxDel, o - 1);
      size_t r = matchRows[o];
      float e = es[o];
      const float *hO = &history[o * stride];
      if (nIns >= STEADY_LANES || nDel >= STEADY_LANES) {
        SteadyOrder cells = {seq, i,    o,     r,     cur, nIns,
                             nDel, e,   tInMa, tDeMa, hO,  &pe};
        cells.Fill(p, c, best[o], bestFrom[o]);
        continue;
      }
      // Fewer cells than a block: one at a time, inline
      const float *atO = pe.Against(seq[i - o]);
      float b = best[o];
      int fromRow = bestFrom[o];
      for (size_t k = 1; k <= nIns; ++k) {
        size_t row = r + k;
        float pr = p[row];
        float exit = pr + tInMa + e;
        if (exit > b) {
          b = exit;
          fromRow = (int)row;
        }
        float v = pr;
        v -= hO[cur - o - k - 1];
        v += hO[cur - o - k];
        v -= atO[i - 2 * o - k];
        v += emitI[i - o - k];
        c[row] = v;
      }
      if (nDel > 0) {
        float parentO = hO[cur - o];
        float parentO1 = hO[cur - o - 1];
        const float *byO = pe.From(seq[i - o]);
        const float *by2O = pe.From(seq[i - 2 * o]);
        for (size_t k = 1; k <= nDel; ++k) {
          size_t row = r + nIns + k;
          float pr = p[row];
          float exit = pr + tDeMa + e;
          if (exit > b) {
            b = exit;
            fromRow = (int)row;
          }
          float v = pr;
          v -= parentO1;
          v += parentO;
          v -= atO[i - 2 * o + k];
          v += byO[i - k];
          v -= by2O[i - k];
          v += emitI[i - 2 * o + k];
          c[row] = v;
        }
      }
      best[o] = b;
      bestFrom[o] = fromRow;
    }
  }

  // Pass 3: the match rows, their history and their tracebacks.
  for (size_t o = 1; o <= nOrders; ++o) {
    c[matchRows[o]] = best[o];
    history[o * stride + cur] = best[o];
    ct[o] = bestFrom[o];
  }
  c[0] = bg;
  ct[0] = bgFrom;
  history[cur] = bg;
}

// The general path of CalculateCurrentColumn for rows startRow.. of the
// current column, looking up each cell's description. Rows before startRow
// must already hold their scores for this column. Does not record the
// column's history or move the matrix forward.
void UModel::CalculateGeneralRows(SequenceWindow *sequence, int nucIndex,
                                  int startRow) {
  bool cbr = false; // can be repetitive - true if it is possible for next
                    // character to be repetitve
  float bestRepeatScore = NEG_INF;

  cell *desc = matrix->cellDescriptions;
  float *p = matrix->previousScoreColumn;
  float *c = matrix->currentScoreColumn;

  int *ct = matrix->currentTracebackColumn;

  int matind = matrix->previousColumnIndex;

  // CalculateEmissionScores(sequence->symbolFreqs[nucIndex]);

  for (int row = startRow; row < matrix->cellsPerColumn; ++row) {
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
        // The extensions of an insertion of more than one base. Skipped for
        // one base: adding t * 0 changes nothing for a finite t, and with
        // --ii 0 (t = -inf) it made the score NaN, so that even single-base
        // insertions could never be taken.
        if (desc[row].indelNumber > 1)
          score += tscore[CT_INSERTION][CT_INSERTION] *
                   (float)(desc[row].indelNumber - 1);
        // The bases after an indel are matches, so they pay the match to
        // match transition (it was log2(1 - tp_matchToZero), cheaper than a
        // match state, which made a path at period P - 1 with an insertion
        // in every copy outscore the true period P of a long unit)
        score += (double)tscore[CT_MATCH][CT_MATCH] * (float)order;

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
        // As for insertions (--dd 0)
        if (indelNum > 1)
          score += tscore[CT_DELETION][CT_DELETION] * (float)(indelNum - 1);
        // The bases after an indel are matches, so they pay the match to
        // match transition (it was log2(1 - tp_matchToZero), cheaper than a
        // match state, which made a path at period P - 1 with an insertion
        // in every copy outscore the true period P of a long unit)
        score += (double)tscore[CT_MATCH][CT_MATCH] * (float)order;

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
        // The extensions of an insertion of more than one base. Skipped for
        // one base: adding t * 0 changes nothing for a finite t, and with
        // --ii 0 (t = -inf) it made the score NaN, so that even single-base
        // insertions could never be taken.
        if (desc[row].indelNumber > 1)
          score += tscore[CT_INSERTION][CT_INSERTION] *
                   (float)(desc[row].indelNumber - 1);
        // The bases after an indel are matches, so they pay the match to
        // match transition (it was log2(1 - tp_matchToZero), cheaper than a
        // match state, which made a path at period P - 1 with an insertion
        // in every copy outscore the true period P of a long unit)
        score += (double)tscore[CT_MATCH][CT_MATCH] * (float)order;

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
        // As for insertions (--dd 0)
        if (indelNum > 1)
          score += tscore[CT_DELETION][CT_DELETION] * (float)(indelNum - 1);
        // The bases after an indel are matches, so they pay the match to
        // match transition (it was log2(1 - tp_matchToZero), cheaper than a
        // match state, which made a path at period P - 1 with an insertion
        // in every copy outscore the true period P of a long unit)
        score += (double)tscore[CT_MATCH][CT_MATCH] * (float)order;

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

  matrix->RecordHistory(matrix->currentColumnIndex);
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
