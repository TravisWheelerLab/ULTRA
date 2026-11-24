/*
 umatrix.hpp

 Contains definitions for cell information and the general matrix used by ULTRA

 UMatrix contains two seperate matrix spaces:
    a float based score matrix used to store scores
    an int based transition matrix used to store transitions
        The transition matrix does not describe transitions from every cell
        instead it only describes transitions from CT_BACKGROUND cells and
 CT_MATCH this allows for some memory saving

 The length of the matrix does not correspond to the length of the
 sequence being analyzed. After ~<length> number of symbosl has been
 analyzed the matrix rotates so that the current column becomes zero again.

 */

#ifndef umatrix_hpp
#define umatrix_hpp

#include "Symbol.hpp"
#include <stdlib.h>
#include <vector>

/*** CELL AND CELL_TYPE ***/

#define NEG_INF (-10000000000000.0)

typedef unsigned char cell_type;

#define CT_BACKGROUND 0
#define CT_MATCH 1
#define CT_INSERTION 2
#define CT_DELETION 3

typedef struct {
  cell_type type;  // What type of cell is this
  int parentIndex; // Points to parent row
  int order;       // Period (or parent period in case of indels)
  int indelNumber; // Which indel does this state represent
} cell;

/*** STEADY-STATE BOUNDS ***/

// The bounds of SteadyState for the cells of a prefix of orders.
typedef struct {
  long long minPreviousColumn;
  long long maxLookback;
  long long maxSequenceBack;
} OrderBounds;

// Bounds past which no cell of a column takes a warm-up branch, every
// look-back column lies at or before the current one without wrapping, and
// every sequence index is in range. Derived from the cell layout. A match
// cell of order o needs previousColumnIndex >= o, looks back o - 1 columns
// and reads o symbols back. An insertion (o, k) needs previousColumnIndex >
// 2o + k + 1, looks back o + k + 1 and reads 2o + k back. A deletion (o, k)
// needs previousColumnIndex > 3o, looks back o + 1 and reads 2o back.
typedef struct {
  // The previous column index must be greater than this.
  long long minPreviousColumn;
  // The largest d any cell passes to PreviousScore.
  long long maxLookback;
  // The largest distance any cell reads back from the current symbol.
  long long maxSequenceBack;
  // The same bounds for the cells of orders 1..o only, entry o (entry 0 is
  // all zeros). Each is non-decreasing in o, so in a column that is not yet
  // steady the orders that are form a prefix.
  std::vector<OrderBounds> byOrder;
} SteadyState;

/*** UMATRIX CLASS ***/

// Scores are kept in two parts:
// - every row of the previous and the current column (fullColumns: column c
//   is slot c % 2), which is all the next column reads besides
// - the background and match rows of every column (history), stored by
//   order: order o (0 is the background row) of column c is at
//   history[o * historyStride + c], so that the look-backs of one order's
//   indel cells are contiguous.
// Nothing reads the insertion and deletion rows of older columns, unless
// KeepAllRows() is called (UModelDriver does), which keeps every row of every
// column as ULTRA 1.2 did. This is 10 to 25 times less memory than keeping
// every row.
class UMatrix {
public:
  // Class variables
  unsigned long long maxPeriod;     // Max period is included
  unsigned long long maxInsertions; // Max insertion is included
  unsigned long long maxDeletions;  // Max deletion is included

  unsigned long long length;         // Number of columns in matrix
  unsigned long long cellsPerColumn; // Number of rows in matrix

  float *fullColumns = NULL;     // 2 * cellsPerColumn
  float *history = NULL;         // (maxPeriod + 1) * historyStride
  unsigned long long historyStride = 0; // length + 1
  int *matchRows = NULL;         // The row of each order's match cell
                                 // (maxPeriod entries; entry 0 is row 0)
  float *allRows = NULL;         // Every row of every column, only after
                                 // KeepAllRows()
  unsigned long long allRowsSince = 0;

  SteadyState steady;

  int *tracebackMatrix = NULL;
  int **tracebackColumns = NULL;

  unsigned long long currentColumnIndex;
  unsigned long long previousColumnIndex;

  float *currentScoreColumn = NULL;
  float *previousScoreColumn = NULL;

  int *currentTracebackColumn = NULL;
  int *previousTracebackColumn = NULL;

  cell *cellDescriptions = NULL; // Array of size <cellsPerColumn>

  float lastScoreAdjustment = 0;
  float totalScoreAdjustment = 0;

  int *traceback = NULL;

  // Methods

  // Keep every row of every column from now on, so that ScoreAt can return
  // the insertion and deletion rows of any column. Costs cellsPerColumn *
  // (length + 1) floats.
  void KeepAllRows();
  bool KeepsAllRows() const { return allRows != NULL; }

  // Whether every row of column col is kept: the previous or the current
  // column.
  bool IsFullColumn(unsigned long long col) const {
    return col == currentColumnIndex || col == previousColumnIndex;
  }
  // Every row of the previous or the current column.
  float *ScoreColumn(unsigned long long col) const {
    return fullColumns + (col & 1) * cellsPerColumn;
  }
  // The background (order 0) or order-order match score of column col.
  float HistoryAt(unsigned long long col, int order) const {
    return history[(unsigned long long)order * historyStride + col];
  }
  // Copies the background and match rows of column col, which must be
  // complete and full, into history.
  void RecordHistory(unsigned long long col);
  // The score of row in column col. Background and match rows are kept for
  // every column; insertion and deletion rows only for the previous and the
  // current column (or every column after KeepAllRows).
  float ScoreAt(unsigned long long col, int row) const;

  bool MoveMatrixForward();

  void AdjustScoreToZero(int startingAt, int numberOfColumns);

  void AdjustScoreToZero(int numberOfColumns);
  void AdjustColumn(int column,
                    float adjustment); // row value = value + adjustment

  float PreviousScore(int row,
                      int d); // How far back we look
  int PreviousTraceback(int row, int d);
  int PreviousColumn(int d);

  int *ForwardTraceback(int *traceArray, int windowLength, int row);

  void RestartMatrix();
  void CalculateTraceback(unsigned long long startColumn);

  // Class management
  void CreateMatrix();

  UMatrix(int maximumPeriod, int maximumInsertions, int maximumDeletions,
          int matrixLength);
  ~UMatrix();
};

#endif /* umatrix_hpp */
