/*
 umatrix.cpp

 */

#include "umatrix.hpp"
#include <algorithm>
#include <stdio.h>

/*
 MoveMatrixForward moves the matrix forward while accounting for wrapping:
 1. Assign previousColumn index to current column
 2. Increment current column
 3. Check if we need to wrap
 4. Assign current/previous column pointers
 5. return wrap
 */
bool UMatrix::MoveMatrixForward() {

  bool wrap = false;

  // The current column is complete.
  if (allRows != NULL) {
    float *full = ScoreColumn(currentColumnIndex);
    float *dst = allRows + currentColumnIndex * cellsPerColumn;
    for (unsigned long long r = 0; r < cellsPerColumn; ++r)
      dst[r] = full[r];
  }

  // Move forward our column indexes (this never wraps)
  previousColumnIndex = currentColumnIndex;
  ++currentColumnIndex;

  // Assign our column pointers to the new columns
  currentScoreColumn = ScoreColumn(currentColumnIndex);
  previousScoreColumn = ScoreColumn(previousColumnIndex);

  currentTracebackColumn = tracebackColumns[currentColumnIndex];
  previousTracebackColumn = tracebackColumns[previousColumnIndex];

  return wrap;
}

void UMatrix::KeepAllRows() {
  if (allRows == NULL) {
    allRows = (float *)calloc(cellsPerColumn * (length + 1), sizeof(float));
    allRowsSince = currentColumnIndex;
  }
}

void UMatrix::RecordHistory(unsigned long long col) {
  float *full = ScoreColumn(col);
  for (unsigned long long o = 0; o < maxPeriod; ++o)
    history[o * historyStride + col] = full[matchRows[o]];
}

float UMatrix::ScoreAt(unsigned long long col, int row) const {
  cell d = cellDescriptions[row];
  if (d.type == CT_MATCH || d.type == CT_BACKGROUND)
    return HistoryAt(col, d.order);
  if (allRows != NULL && col >= allRowsSince)
    return allRows[col * cellsPerColumn + row];
  if (IsFullColumn(col))
    return ScoreColumn(col)[row];
  fprintf(stderr,
          "UMatrix: column %llu keeps only its background and match rows; "
          "call KeepAllRows() before computing the window to read its "
          "insertion and deletion rows\n",
          col);
  abort();
}

/*
 AdjustScoreToZero adjusts every column so that the scores are closer to zero
 (to ensure precision): every score is shifted by minus the previous
 column's background score. numberOfColumns is unused (every column is
 adjusted).
 */
void UMatrix::AdjustScoreToZero(int numberOfColumns) {
  // Calculate adjustment
  lastScoreAdjustment = previousScoreColumn[0] * -1;
  totalScoreAdjustment += lastScoreAdjustment;

  float adj = lastScoreAdjustment;
  for (unsigned long long o = 0; o < maxPeriod + 1; ++o) {
    float *h = history + o * historyStride;
    for (unsigned long long c = 0; c < length; ++c)
      h[c] += adj;
  }
  unsigned long long cols[2] = {previousColumnIndex, currentColumnIndex};
  for (int k = 0; k < 2; ++k) {
    if (cols[k] < length) {
      float *f = ScoreColumn(cols[k]);
      for (unsigned long long r = 0; r < cellsPerColumn; ++r)
        f[r] += adj;
    }
  }
  if (allRows != NULL) {
    for (unsigned long long i = 0; i < cellsPerColumn * length; ++i)
      allRows[i] += adj;
  }
}

// Resets numberOfColumns columns from startingAt to row 0 = 0 and every other
// row = NEG_INF
void UMatrix::AdjustScoreToZero(int startingAt, int numberOfColumns) {
  for (int i = 0; i < numberOfColumns; ++i) {
    unsigned long long p = i + startingAt;
    for (unsigned long long o = 0; o < maxPeriod + 1; ++o)
      history[o * historyStride + p] = o == 0 ? 0 : NEG_INF;
    if (IsFullColumn(p)) {
      float *col = ScoreColumn(p);
      for (unsigned long long j = 1; j < cellsPerColumn; ++j)
        col[j] = NEG_INF;
      col[0] = 0;
    }
    if (allRows != NULL) {
      float *col = allRows + p * cellsPerColumn;
      for (unsigned long long j = 1; j < cellsPerColumn; ++j)
        col[j] = NEG_INF;
      col[0] = 0;
    }
  }
}

// Adjust every row in <column> to be [row] = [row] + adjustment
void UMatrix::AdjustColumn(int column, float adjustment) {
  unsigned long long c = column;
  for (unsigned long long o = 0; o < maxPeriod + 1; ++o)
    history[o * historyStride + c] += adjustment;
  if (IsFullColumn(c)) {
    float *f = ScoreColumn(c);
    for (unsigned long long r = 0; r < cellsPerColumn; ++r)
      f[r] += adjustment;
  }
  if (allRows != NULL) {
    float *f = allRows + c * cellsPerColumn;
    for (unsigned long long r = 0; r < cellsPerColumn; ++r)
      f[r] += adjustment;
  }
}

// Return the score of row in the column d before the current one (d = 0 is
// the current column, read as it is being filled). Account for wrapping if
// currentColumnIndex - d < 0
float UMatrix::PreviousScore(int row, int d) {
  float s = NEG_INF;

  if (d < length) {
    unsigned long long c = PreviousColumn(d);
    if (c == currentColumnIndex)
      s = currentScoreColumn[row];
    else
      s = ScoreAt(c, row);
  }

  return s;
}

// Return the traceback in tracebackColumns[currentColumnIndex - d][row]
// Account for warpping if currentColumnIndex - d < 0
int UMatrix::PreviousTraceback(int row, int d) {
  int r = 0;
  if (d < length) {
    int c = PreviousColumn(d);

    r = tracebackColumns[c][row];
  }

  return r;
}

// Calculate currentColumnIndex - d while accounting for wrapping
int UMatrix::PreviousColumn(int d) {
  d = currentColumnIndex - d;
  if (d < 0) { // Check if we have to cycle to the other end of the matrix
    d += length;
  }
  return d;
}

/*
 ForwardTraceback returns a forward oriented traceback and
 makes the necessary adjustments for insertions and deletions:
 1. Create an array to hold the traceback if necessary
 2. Follow the traceback and assign the rows to the traceArray
 3. If we encounter an insertion, the actual insertion
        occurred (order + indelNum) previously, set the previous
        ((order + indelNum) - 1) values to be the parentIndex
        and then the pos - (order + indelNum) to be the correct indel row
 4. If we encounter a deletion, the actual deletion occured (order) previously.
    Continue as in (3), but using (order) instead of (order + indelnum)

 */
int *UMatrix::ForwardTraceback(int *traceArray, int windowLength, int row) {
  // The traceArray is forward oriented, and has correct insertion and
  // deletion positions
  if (windowLength > length) {
    printf("(ForwardTraceback) Length greater than matrix length, ret NULL\n");
    return NULL;
  }

  // If traceArray is NULl then we need to allocate the memory for our
  // traceArray
  if (traceArray == NULL) {
    traceArray = (int *)malloc(sizeof(int) * (windowLength + 1));
  }

  int pos = previousColumnIndex; // pos is where in the tracebackMatrix we are
                                 // looking

  for (int i = 0; i <= windowLength; ++i) {

    // Get some info on the current cell - for each tracebackColumn there
    // are only N rows, where N is the maximum peridoicity
    cell desc = cellDescriptions[row];
    int rowOrder = desc.order;

    // In case of C_MATCH/C_NONE we can procede normally
    if (desc.type == CT_MATCH || desc.type == CT_BACKGROUND) {
      traceArray[windowLength - i] = row;
      row = tracebackColumns[pos][rowOrder];
    }

    // In case of C_INSERTION we need to set the next (order + indelnum)
    // characters to be the cell of the insertions parent cell and then set
    // the last character to the insertion. We do things in this order
    // because insertions are detected (order + indelnum) after they occur,
    // and we are reading the traceback matrix in reverse order.
    else if (desc.type == CT_INSERTION) {

      int charactersBack = desc.order + desc.indelNumber;
      for (int j = 0; j < charactersBack; ++j) {

        traceArray[windowLength - i] = desc.parentIndex;
        ++i;

        if (i > windowLength) {
          return traceArray;
        }

        pos = pos - 1;
        if (pos < 0) {
          pos = length - 1;
        }
      }

      traceArray[windowLength - i] = row;
      row = desc.parentIndex;
    }

    // In case of C_DELETION we need to set the next (order) characters to
    // be the cell of the deletion's parent cell
    // and then set the last character to the actual deletion.
    else if (desc.type == CT_DELETION) {
      int charactersBack = desc.order;

      for (int j = 0; j < charactersBack; ++j) {

        traceArray[windowLength - i] = desc.parentIndex;
        ++i;

        if (i > windowLength) {

          return traceArray;
        }
        pos = pos - 1;
        if (pos < 0) {
          pos = length - 1;
        }
      }

      traceArray[windowLength - i] = row;
      row = desc.parentIndex;
    }

    // Whenever we decrease pos we need to check and make sure we haven't
    // gone past the matrix bounds if we have then we need to adjust for
    // that
    pos = pos - 1;
    if (pos < 0) {
      pos = length - 1;
    }
  }

  return traceArray;
}

// This has not been debugged yet
void UMatrix::CalculateTraceback(unsigned long long startColumn) {
  // Assume best row is 0
  unsigned long long row = 0;
  float best_value = HistoryAt(startColumn, 0);

  // The match rows in row order, which is order order.
  for (unsigned long long o = 1; o < maxPeriod; ++o) {
    float v = HistoryAt(startColumn, (int)o);
    if (v > best_value) {
      row = matchRows[o];
      best_value = v;
    }
  }
  // Do the normal calculations
  for (unsigned long long i = 0; i <= startColumn; ++i) {

    cell desc = cellDescriptions[row];
    int rowOrder = desc.order;

    if (desc.type == CT_MATCH || desc.type == CT_BACKGROUND) {
      traceback[startColumn - i] = row;
      row = tracebackColumns[startColumn - i][rowOrder];
    }

    else if (desc.type == CT_INSERTION) {
      int charactersBack = desc.order + desc.indelNumber;
      for (int j = 0; j < charactersBack; ++j) {
        traceback[startColumn - i] = desc.parentIndex;
        ++i;

        if (i >= startColumn)
          return;
      }

      traceback[startColumn - i] = row;
      row = desc.parentIndex;

    }

    else if (desc.type == CT_DELETION) {
      int charactersBack = desc.order;

      for (int j = 0; j < charactersBack; ++j) {
        traceback[startColumn - i] = desc.parentIndex;
        ++i;

        if (i >= startColumn)
          return;
      }

      traceback[startColumn - i] = row;
      row = desc.parentIndex;
    }
  }
}

void UMatrix::RestartMatrix() {
  currentColumnIndex = 1;
  previousColumnIndex = 0;

  previousScoreColumn = ScoreColumn(0);
  currentScoreColumn = ScoreColumn(1);

  previousTracebackColumn = tracebackColumns[0];
  currentTracebackColumn = tracebackColumns[1];

  for (unsigned long long i = 0; i < cellsPerColumn; ++i) {
    previousScoreColumn[i] = 0;
  }
  for (unsigned long long o = 0; o < maxPeriod + 1; ++o) {
    history[o * historyStride] = 0;
  }

  for (unsigned long long i = 0; i < maxPeriod + 1; ++i) {
    previousTracebackColumn[i] = 0;
  }

  if (allRows != NULL) {
    for (unsigned long long i = 0; i < cellsPerColumn; ++i)
      allRows[i] = 0;
    allRowsSince = 0;
  }
}

/** UMATRIX CLASS MANAGEMENT ***/

/*
 CreateMatrix allocates and initializes resources used by UMatrix:
 1. Count the number of cells in a column.
 2. Create the cellDescriptions array
 3. Allocate matrices
 4. Assign columns
 5. Initialize the 0th columns
 */

void UMatrix::CreateMatrix() {

  // Creation of cell descriptions
  int numberOfCells = 0;

  // Calculate the number of cells needed
  for (int i = 0; i < maxPeriod; ++i) {
    ++numberOfCells;

    if (i <= maxInsertions) {
      if (i > 1) {
        numberOfCells += i - 1;
      }
    }

    else {
      numberOfCells += maxInsertions;
    }

    if (i <= maxDeletions) {
      if (i > 1) {
        numberOfCells += i - 1;
      }
    }

    else {
      numberOfCells += maxDeletions;
    }
  }

  // Allocate the cell description
  cellsPerColumn = numberOfCells;
  cellDescriptions = (cell *)malloc(sizeof(cell) * (cellsPerColumn + 1));
  // printf("%llx to %llx\n", (unsigned long)cellDescriptions, (unsigned
  // long)cellDescriptions + sizeof(cell) * (cellsPerColumn + 1));
  // Calculate cell descriptions
  int cellIndex = 1;

  // 0th order cell description
  cellDescriptions[0].type = CT_BACKGROUND;
  cellDescriptions[0].parentIndex = 0;
  cellDescriptions[0].order = 0;
  cellDescriptions[0].indelNumber = 0;

  for (int i = 1; i < maxPeriod; ++i) {
    int matchIndex = cellIndex;

    cellDescriptions[cellIndex].type = CT_MATCH;
    cellDescriptions[cellIndex].parentIndex = 0;
    cellDescriptions[cellIndex].order = i;
    cellDescriptions[cellIndex].indelNumber = 0;

    ++cellIndex;

    // Calculate number of indels
    int numIns = maxInsertions;
    int numDel = maxDeletions;

    if (numIns >= i) {
      numIns = i - 1;
    }

    if (numDel >= i) {
      numDel = i - 1;
    }

    // Calculate insertion descriptions
    for (int j = 0; j < numIns; ++j) {
      cellDescriptions[cellIndex].type = CT_INSERTION;
      cellDescriptions[cellIndex].order = i;
      cellDescriptions[cellIndex].parentIndex = matchIndex;
      cellDescriptions[cellIndex].indelNumber = j + 1;
      cellIndex++;
    }

    // Calculate deletion descriptions
    for (int j = 0; j < numDel; ++j) {
      cellDescriptions[cellIndex].type = CT_DELETION;
      cellDescriptions[cellIndex].order = i;
      cellDescriptions[cellIndex].parentIndex = matchIndex;
      cellDescriptions[cellIndex].indelNumber = j + 1;
      cellIndex++;
    }
  }

  cellsPerColumn = cellIndex;

  /* for (int i = 0; i < cellsPerColumn; ++i) {
       printf("%i: (T: %i) (O: %i) (I: %i) (P: %i)\n",
              i,
              cellDescriptions[i].type,
              cellDescriptions[i].order,
              cellDescriptions[i].indelNumber,
              cellDescriptions[i].parentIndex);
   }*/

  // The row of each order's match cell
  matchRows = (int *)malloc(sizeof(int) * maxPeriod);
  matchRows[0] = 0;
  for (unsigned long long r = 1; r < cellsPerColumn; ++r) {
    if (cellDescriptions[r].type == CT_MATCH)
      matchRows[cellDescriptions[r].order] = (int)r;
  }

  // Steady-state bounds (see SteadyState)
  steady.minPreviousColumn = 0;
  steady.maxLookback = 0;
  steady.maxSequenceBack = 0;
  steady.byOrder.assign(maxPeriod > 0 ? maxPeriod : 1, OrderBounds{0, 0, 0});
  for (unsigned long long r = 0; r < cellsPerColumn; ++r) {
    cell d = cellDescriptions[r];
    long long o = d.order, k = d.indelNumber;
    long long minPrev = 0, lookback = 0, seqBack = 0;
    if (d.type == CT_MATCH) {
      minPrev = o - 1;
      lookback = o - 1;
      seqBack = o;
    } else if (d.type == CT_INSERTION) {
      minPrev = 2 * o + k + 1;
      lookback = o + k + 1;
      seqBack = 2 * o + k;
    } else if (d.type == CT_DELETION) {
      minPrev = 3 * o;
      lookback = o + 1;
      seqBack = 2 * o;
    }
    steady.minPreviousColumn = std::max(steady.minPreviousColumn, minPrev);
    steady.maxLookback = std::max(steady.maxLookback, lookback);
    steady.maxSequenceBack = std::max(steady.maxSequenceBack, seqBack);
    OrderBounds &b = steady.byOrder[d.order];
    b.minPreviousColumn = std::max(b.minPreviousColumn, minPrev);
    b.maxLookback = std::max(b.maxLookback, lookback);
    b.maxSequenceBack = std::max(b.maxSequenceBack, seqBack);
  }
  // Prefix maxima: entry o covers orders 1..o
  for (size_t o = 1; o < steady.byOrder.size(); ++o) {
    OrderBounds p = steady.byOrder[o - 1];
    OrderBounds &b = steady.byOrder[o];
    b.minPreviousColumn = std::max(b.minPreviousColumn, p.minPreviousColumn);
    b.maxLookback = std::max(b.maxLookback, p.maxLookback);
    b.maxSequenceBack = std::max(b.maxSequenceBack, p.maxSequenceBack);
  }

  // Allocate matrices
  historyStride = length + 1;
  fullColumns = (float *)calloc(2 * cellsPerColumn, sizeof(float));
  history = (float *)calloc((maxPeriod + 1) * historyStride, sizeof(float));
  tracebackMatrix = (int *)malloc(sizeof(int) * (maxPeriod + 1) * (length + 1));
  tracebackColumns = (int **)malloc(sizeof(int *) * (length + 1));

  // Assign columns of the traceback matrix
  for (unsigned long long i = 0; i <= length; ++i) {
    tracebackColumns[i] = &(tracebackMatrix[(i * (maxPeriod + 1))]);
  }

  previousScoreColumn = ScoreColumn(0);
  currentScoreColumn = ScoreColumn(1);

  previousTracebackColumn = tracebackColumns[0];
  currentTracebackColumn = tracebackColumns[1];

  for (int i = 0; i < maxPeriod + 1; ++i) {
    previousTracebackColumn[i] = 0;
  }

  traceback = (int *)malloc(sizeof(int) * (length + 2));

  for (int i = 0; i < length + 2; ++i) {
    traceback[i] = 0;
  }

  currentColumnIndex = 1;
  previousColumnIndex = 0;
}

UMatrix::UMatrix(int maximumPeriod, int maximumInsertions, int maximumDeletions,
                 int matrixLength) {
  maxPeriod = maximumPeriod;
  // A negative indel count builds no indel cells
  maxInsertions = maximumInsertions > 0 ? maximumInsertions : 0;
  maxDeletions = maximumDeletions > 0 ? maximumDeletions : 0;
  length = matrixLength;

  this->CreateMatrix();
}

UMatrix::~UMatrix() {
  free(fullColumns);
  fullColumns = NULL;
  free(history);
  history = NULL;
  free(matchRows);
  matchRows = NULL;
  free(allRows);
  allRows = NULL;

  if (tracebackMatrix != NULL) {
    free(tracebackMatrix);
    tracebackMatrix = NULL;
  }

  if (tracebackColumns != NULL) {
    free(tracebackColumns);
    tracebackColumns = NULL;
  }

  if (traceback != NULL) {
    free(traceback);
    traceback = NULL;
  }

  currentScoreColumn = NULL;
  previousScoreColumn = NULL;
  currentTracebackColumn = NULL;
  previousTracebackColumn = NULL;

  if (cellDescriptions != NULL) {
    free(cellDescriptions);
    cellDescriptions = NULL;
  }
}
