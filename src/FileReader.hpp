//
//  FileReader.hpp
//  ultrax
//

#ifndef FileReader_hpp
#define FileReader_hpp

#include "FASTAReader.hpp"
#include <stdio.h>

typedef enum { UNKNOWN, FASTA, JSON, BED } file_type;

class FileReader {
public:
  int maxWindows;
  unsigned long maxSeqLength;
  unsigned long maxOverlapLength;

  bool multithread = true;
  file_type format = UNKNOWN;

  FASTAReader *fastaReader = NULL;

  SequenceWindow *GetReadyWindow();
  bool AddReadyWindow(SequenceWindow *window);

  SequenceWindow *GetWaitingWindow();
  bool AddWaitingWindow(SequenceWindow *window);

  bool FillWindows();

  bool IsReading();
  bool DoneReadingFile();

  unsigned long ReadyWindowsSize();

  void SetIsReading(bool value);

  FileReader(std::string FASTAFilePath, int maxWindows,
             unsigned long maxSeqLength, unsigned long maxOverlapLength,
             bool multithread);

  FileReader(std::string JSONFilePath, unsigned long maxSeqLength,
             unsigned long maxOverlapLength, bool multithread);
};

#endif /* FileReader_hpp */
