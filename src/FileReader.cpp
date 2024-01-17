//
//  FileReader.cpp
//  ultrax
//

#include "FileReader.hpp"

FileReader::FileReader(std::string FASTAFilePath, int mWindows,
                       unsigned long mSeqLength, unsigned long mOverlapLength,
                       bool mthread) {
  maxWindows = mWindows;
  maxSeqLength = mSeqLength;
  maxOverlapLength = mOverlapLength;
  multithread = mthread;
  format = FASTA;

  fastaReader = new FASTAReader(FASTAFilePath, maxWindows, maxSeqLength,
                                maxOverlapLength);
  fastaReader->multithread = multithread;
}

FileReader::FileReader(std::string JSONFilePath, unsigned long mSeqLength,
                       unsigned long mOverlapLength, bool mthread) {

  maxSeqLength = mSeqLength;
  maxOverlapLength = mOverlapLength;
  multithread = mthread;

  format = JSON;
  printf("JSON input no longer allowed.\n");
  exit(-1);
}

FileReader::~FileReader() {
  if (fastaReader) {
    delete fastaReader;
    fastaReader = nullptr;
  }

}

SequenceWindow *FileReader::GetReadyWindow() {
  if (format == FASTA)
    return fastaReader->GetReadyWindow();

  else {
    printf("SequenceWindow Error, invalid format type\n");
    exit(-1);
  }

  return NULL;
}

bool FileReader::AddReadyWindow(SequenceWindow *window) {
  if (format == FASTA)
    return fastaReader->AddReadyWindow(window);

  return false;
}

SequenceWindow *FileReader::GetWaitingWindow() {
  if (format == FASTA)
    return fastaReader->GetWaitingWindow();

  return NULL;
}

bool FileReader::AddWaitingWindow(SequenceWindow *window) {
  if (format == FASTA)
    return fastaReader->AddWaitingWindow(window);

  return false;
}

bool FileReader::IsReading() {
  if (format == JSON) {
    return false;
  }

  return fastaReader->isReading;
}

bool FileReader::DoneReadingFile() {
  if (format == JSON) {
    return true;
  }

  return fastaReader->doneReadingFile;
}

void FileReader::SetIsReading(bool value) {
  if (format == FASTA)
    fastaReader->isReading = value;
}

unsigned long FileReader::ReadyWindowsSize() {
  if (format == FASTA)
    return fastaReader->readyWindows.size();

  return 0;
}

bool FileReader::FillWindows() {
  if (format == FASTA)
    return fastaReader->FillWindows();

  else
    return true;
}
