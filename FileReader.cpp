//
//  FileReader.cpp
//  ultrax
//
//

#include "FileReader.hpp"


FileReader::FileReader(std::string FASTAFilePath,
           unsigned long mWindows,
           unsigned long mSeqLength,
           unsigned long mOverlapLength,
                    bool mthread)
{
    maxWindows = mWindows;
    maxSeqLength = mSeqLength;
    maxOverlapLength = mOverlapLength;
    multithread = mthread;
    format = FASTA;
    
    fastaReader = new FASTAReader(FASTAFilePath, maxWindows, maxSeqLength, maxOverlapLength);
    fastaReader->multithread = multithread;
    
}

FileReader::FileReader(std::string JSONFilePath,
           unsigned long mSeqLength,
           unsigned long mOverlapLength,
           bool mthread)
{
   
    maxSeqLength = mSeqLength;
    maxOverlapLength = mOverlapLength;
    multithread = mthread;
    
    format = JSON;
    jsonReader = new JSONReader(JSONFilePath, maxSeqLength, maxOverlapLength);
    jsonReader->multithread = multithread;
}


SequenceWindow *FileReader::GetReadyWindow() {
    if (format == FASTA)
        return fastaReader->GetReadyWindow();
    
    
    else if (format == JSON)
        return jsonReader->GetReadyWindow();
    
    return NULL;
}

bool FileReader::AddReadyWindow(SequenceWindow *window) {
    if (format == FASTA)
        return fastaReader->AddReadyWindow(window);
    
    
    else if (format == JSON)
        return jsonReader->AddReadyWindow(window);
    
    return false;
}

SequenceWindow *FileReader::GetWaitingWindow() {
    if (format == FASTA)
        return fastaReader->GetWaitingWindow();
    
    
    else if (format == JSON)
        return jsonReader->GetWaitingWindow();
    
    return NULL;
}

bool FileReader::AddWaitingWindow(SequenceWindow *window) {
    if (format == FASTA)
        return fastaReader->AddWaitingWindow(window);
    
    
    else if (format == JSON)
        return jsonReader->AddWaitingWindow(window);
    
    return false;
}

bool FileReader::IsReading() {
    if (format==JSON) {
        return false;
    }
    
    return fastaReader->isReading;
}

bool FileReader::DoneReadingFile() {
    if (format==JSON) {
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
    else
        return jsonReader->readyWindows.size();
}

bool FileReader::FillWindows() {
    if (format == FASTA)
        return fastaReader->FillWindows();
    
    else
        return true;
}
