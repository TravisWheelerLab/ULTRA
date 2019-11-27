//
//  FileReader.hpp
//  ultrax
//
//  Created by Daniel Olson on 10/16/19.
//  Copyright © 2019 Daniel Olson. All rights reserved.
//

#ifndef FileReader_hpp
#define FileReader_hpp

#include <stdio.h>
#include "FASTAReader.hpp"
#include "JSONReader.hpp"

typedef enum {
    UNKNOWN,
    FASTA,
    JSON
} file_type;

class FileReader {
public:
    
    unsigned long       maxWindows;
    unsigned long       maxSeqLength;
    unsigned long       maxOverlapLength;
    
    bool                multithread         = true;
    file_type           format              = UNKNOWN;
    
    FASTAReader*        fastaReader         = NULL;
    JSONReader*         jsonReader          = NULL;
    
    SequenceWindow *GetReadyWindow();
    bool AddReadyWindow(SequenceWindow* window);
    
    SequenceWindow *GetWaitingWindow();
    bool AddWaitingWindow(SequenceWindow* window);
    
    bool FillWindows();
    
    bool IsReading();
    bool DoneReadingFile();
    
    unsigned long ReadyWindowsSize();
    
    void SetIsReading(bool value);
    
    
    
    FileReader(std::string FASTAFilePath,
               unsigned long maxWindows,
               unsigned long maxSeqLength,
               unsigned long maxOverlapLength,
               bool multithread);
    
    FileReader(std::string JSONFilePath,
               unsigned long maxSeqLength,
               unsigned long maxOverlapLength,
               bool multithread);
    
};

#endif /* FileReader_hpp */
