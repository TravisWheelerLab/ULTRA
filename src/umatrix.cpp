/*
 umatrix.cpp

 */

#include "umatrix.hpp"

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

  // Move forward our column indexes
  previousColumnIndex = currentColumnIndex;
  ++currentColumnIndex;

  // Check to see if we need to wrap around
  /*if (currentColumnIndex >= length) {
      currentColumnIndex = 0;
      wrap = true;
  }*/

  // Assign our column pointers to the new columns
  currentScoreColumn = scoreColumns[currentColumnIndex];
  previousScoreColumn = scoreColumns[previousColumnIndex];

  currentTracebackColumn = tracebackColumns[currentColumnIndex];
  previousTracebackColumn = tracebackColumns[previousColumnIndex];

  return wrap;
}

/*
 AdjustScoreToZero adjusts the previous numberOfColumns columns so that
 the scores are closer to zero (to ensure precision):
 1. Calculate how much we are adjusting the score by
 2. Adjust the previous numberOfColumns columns by that much
 */
void UMatrix::AdjustScoreToZero(int numberOfColumns) {
  // Calculate adjustment
  lastScoreAdjustment = previousScoreColumn[0] * -1;
  totalScoreAdjustment += lastScoreAdjustment;

  // Adjust columns
  for (int i = 0; i < length; ++i) {
    int c = PreviousColumn(i + 1);
    AdjustColumn(c, lastScoreAdjustment);
  }
}

void UMatrix::AdjustScoreToZero(int startingAt, int numberOfColumns) {
  for (int i = 0; i < numberOfColumns; ++i) {
    int p = i + startingAt;
    for (int j = 1; j < cellsPerColumn; ++j) {
      scoreColumns[p][j] = NEG_INF;
    }
    scoreColumns[p][0] = 0;
  }
}

// Adjust every row in <column> to be [row] = [row] + adjustment
void UMatrix::AdjustColumn(int column, float adjustment) {
  for (int i = 0; i < cellsPerColumn; ++i) {
    scoreColumns[column][i] += adjustment;
  }
}

// Return the score in scoreColumns[currentColumnIndex - d][row]
// Account for warpping if currentColumnIndex - d < 0
float UMatrix::PreviousScore(int row, int d) {
  float s = NEG_INF;

  if (d < length) {
    int c = PreviousColumn(d);

    s = scoreColumns[c][row];
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
  float best_value = scoreColumns[startColumn][0];

  for (unsigned long long i = 1; i < cellsPerColumn; ++i) {
    if (cellDescriptions[i].type == CT_MATCH)
      if (scoreColumns[startColumn][i] > best_value) {
        row = i;
        best_value = scoreColumns[startColumn][i];
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

  // AdjustScoreToZero(0, 100);

  /*for (int i = 0; i < maxPeriod; ++i) {
      tracebackColumns[0][i] = 0;
      tracebackColumns[1][i] = 0;
  }*/

  previousScoreColumn = scoreColumns[0];
  currentScoreColumn = scoreColumns[1];

  previousTracebackColumn = tracebackColumns[0];
  currentTracebackColumn = tracebackColumns[1];

  currentColumnIndex = 1;
  previousColumnIndex = 0;

  for (int i = 0; i < cellsPerColumn; ++i) {
    previousScoreColumn[i] = 0;
  }

  for (int i = 0; i < maxPeriod + 1; ++i) {
    previousTracebackColumn[i] = 0;
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

  // Allocate matrices
  scoreMatrix = (float *)malloc(sizeof(float) * cellsPerColumn * (length + 1));
  tracebackMatrix = (int *)malloc(sizeof(int) * (maxPeriod + 1) * (length + 1));
  scoreColumns = (float **)malloc(sizeof(float *) * (length + 1));
  tracebackColumns = (int **)malloc(sizeof(int *) * (length + 1));

  // Assign columns of matrices
  for (int i = 0; i < length; ++i) {
    scoreColumns[i] = &(scoreMatrix[(i * (cellsPerColumn))]);
    tracebackColumns[i] = &(tracebackMatrix[(i * (maxPeriod + 1))]);
  }

  previousScoreColumn = scoreColumns[0];
  currentScoreColumn = scoreColumns[1];

  previousTracebackColumn = tracebackColumns[0];
  currentTracebackColumn = tracebackColumns[1];

  // Zero out 0th columns
  for (int i = 0; i < cellsPerColumn; ++i) {
    previousScoreColumn[i] = 0;
  }

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
  maxInsertions = maximumInsertions;
  maxDeletions = maximumDeletions;
  length = matrixLength;

  this->CreateMatrix();
}

UMatrix::~UMatrix() {
  if (scoreMatrix != NULL) {
    free(scoreMatrix);
    scoreMatrix = NULL;
  }

  if (scoreColumns != NULL) {
    free(scoreColumns);
    scoreColumns = NULL;
  }

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
