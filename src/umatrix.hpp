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

/*** UMATRIX CLASS ***/

class UMatrix {
public:
  // Class variables
  int maxPeriod;      // Max period is included
  int maxInsertions;  // Max insertion is included
  int maxDeletions;   // Max deletion is included

  int length;         // Number of columns in matrix
  int cellsPerColumn; // Number of rows in matrix

  float *scoreMatrix = NULL;
  float **scoreColumns = NULL;

  int *tracebackMatrix = NULL;
  int **tracebackColumns = NULL;

  int currentColumnIndex;
  int previousColumnIndex;

  float *currentScoreColumn = NULL;
  float *previousScoreColumn = NULL;

  int *currentTracebackColumn = NULL;
  int *previousTracebackColumn = NULL;

  cell *cellDescriptions = NULL; // Array of size <cellsPerColumn>

  float lastScoreAdjustment = 0;
  float totalScoreAdjustment = 0;

  int *traceback = NULL;

  // Methods

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
  void CalculateTraceback(int startColumn);

  // Class management
  void CreateMatrix();

  UMatrix(int maximumPeriod, int maximumInsertions, int maximumDeletions,
          int matrixLength);
  ~UMatrix();
};

#endif /* umatrix_hpp */
