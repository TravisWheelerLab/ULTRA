//
//  ultra.cpp
//  ultraP
//


#include "ultra.hpp"

void *UltraThreadLaunch(void *dat) {
   
    uthread *uth = (uthread *)dat;
 //   printf("Launching with %i\n", uth->id);
    uth->ultra->AnalyzeFileWithThread(dat);
    return NULL;
}

void Ultra::AnalyzeFile() {
    
    if (numberOfThreads <= 0) {
        printf("Number of threads must be set greater than 0.\n");
        exit(-1);
    }
    
   // printf("Setting number of threads to %i\n", numberOfThreads);
    reader->FillWindows();
  //  printf("Windows filled\n");
    
    
    
    if (pthread_mutex_init(&outerLock, NULL) != 0) {
        printf("Failed to create outer mutex lock. Exiting.\n");
        exit(-1);
    }
    
    if (pthread_mutex_init(&innerLock, NULL) != 0) {
        printf("Failed to create outer mutex lock. Exiting.\n");
        exit(-1);
    }
    
    if (pthread_mutex_init(&repeatLock, NULL) != 0) {
        printf("Failed to create repeat lock. Exiting.\n");
        exit(-1);
    }
    
    fprintf(out, "Sequence Name,");
    if (outputReadID) {
        fprintf(out, "Window ID,");
    }
    fprintf(out, "Start,End,Period,Score,Mismatches,Insertions,Deletions,Consensus");
    if (outputRepeatSequence) {
        fprintf(out, ",Sequence\n");
    }
    
    
    for (int i = 1; i < numberOfThreads; ++i) {
        //printf("Creating thread with id %i\n", i);
        pthread_create(&threads[i]->p_thread, NULL, UltraThreadLaunch, threads[i]);
    }
    
    UltraThreadLaunch(threads[0]);
    
}

SequenceWindow *Ultra::GetSequenceWindow(SequenceWindow *seq) {
    
    SequenceWindow *retval = NULL;
    bool shouldRead = false;
    
    if (seq != NULL)
        reader->AddWaitingWindow(seq);
    
    if (multithreading)
        pthread_mutex_lock(&outerLock);
    
    // Check to see if we need to read more sequence first
    
    if (!reader->doneReadingFile) {
        if (!reader->isReading) {
            if (reader->readyWindows.size() < minReaderSize) {
                shouldRead = true;
                reader->isReading = true;
            }
        }
    }
    
    if (multithreading)
        pthread_mutex_unlock(&outerLock);
    
    if (shouldRead) {
        if (multithreading)
            pthread_mutex_lock(&innerLock);
        
        reader->FillWindows();
        reader->isReading = false;
        
        if (multithreading)
            pthread_mutex_unlock(&innerLock);
    }
    
    retval = reader->GetReadyWindow();
    
    /*
    SequenceWindow *retval = NULL;
    pthread_mutex_lock(&outerLock);
    bool shouldRead = false;
    
    if (seq != NULL)
        reader->AddWaitingWindow(seq);
    
    // We aren't done yet
    if (reader->doneReadingFile == false || reader->readyWindows.size() > 0) {
        if (reader->doneReadingFile == false &&
            reader->readyWindows.size() < minReaderSize)
        {
            //         printf("%i: Filling windows in critical section\n", tid);
            if (reader->isReading == false) {
                reader->isReading = true;
                shouldRead = true;
            }
            
        }
        
        else {
            // Busy wait un til there is a sequence window ready
            while (reader->doneReadingFile != false && reader->readyWindows.size() == 0) {
                
            }
            
            if (reader->readyWindows.size() > 0) {
                retval = reader->GetReadyWindow();
            }
        }
    }
    
    if (!shouldRead) {
        retval = reader->GetReadyWindow();
    }
    
    
    pthread_mutex_unlock(&outerLock);
    
    if (shouldRead) {
        pthread_mutex_lock(&innerLock);
        
        reader->FillWindows();
        
        pthread_mutex_unlock(&innerLock);
    }
    
    
    */
    return retval;
}

void Ultra::AnalyzeFileWithThread(void *dat) {
    
    uthread *uth = (uthread *)dat;
    int tid = uth->id;
    
    
    SequenceWindow *currentWindow = GetSequenceWindow(NULL);
    while (currentWindow != NULL || !reader->doneReadingFile) {
        //printf("%i: %llx %i\n", i++, (unsigned long long)currentWindow, reader->doneReadingFile);
        
        if (currentWindow == NULL) {
            //printf("Current windowi s NULL...\n");
        }
        
        if (currentWindow != NULL)
            AnalyzeSequenceWindow(currentWindow, uth);
        
        currentWindow = GetSequenceWindow(currentWindow);
    }
    
    if (multithreading) {
        if (tid == 0) {
            for (int i = 1; i < numberOfThreads; ++i) {
                pthread_join(threads[i]->p_thread, NULL);
            }
        }
        
        else {
            pthread_exit(NULL);
        }
        
    }
}

void Ultra::AnalyzeSequenceWindow(SequenceWindow *sequence, uthread *uth) {
    
    int sleng = (int)sequence->length + (int)sequence->overlap;
    //matrix->AdjustScoreToZero(0, 2);
    UModel *model = uth->model;
    model->matrix->RestartMatrix();
    
    /* UMatrix *newmat = new UMatrix(matrix->maxPeriod,
     matrix->maxInsertions,
     matrix->maxDeletions,
     matrix->length);
     delete matrix;
     matrix = newmat;
     model->matrix = matrix;*/
    
    /* printf("Analyzing sequence: ");
     for (int i = 0; i < 20; ++i) {
     printf("%c", CharForSymbol(sequence->seq[i]));
     }
     printf("\n");*/
    
    for (int i = 0; i < sleng; ++i) {
        // printf("%i\n", i);
        model->CalculateCurrentColumn(sequence, i);
    }
    
    model->matrix->CalculateTraceback(model->matrix->previousColumnIndex);
    int i = 0;
    RepeatRegion *r = GetNextRepeat(sequence, model->matrix, &i);
    
    
    while (r != NULL) {
        
        if (outputRepeatSequence) {
            r->StoreSequence(sequence);
        }
        
        uth->repeats.push_back(r);

    
        r = GetNextRepeat(sequence, model->matrix, &i);
        
    }
    
    if (primaryThread == uth->id) {
        
        outRepeats.insert(outRepeats.end(), uth->repeats.begin(), uth->repeats.end());
        uth->repeats.clear();
        
        if (outRepeats.size() > repeatBuffer) {
            OutputRepeats();
        }
        
        primaryThread = (primaryThread + 1) % numberOfThreads;
       // printf("%i\n", primaryThread);
    }
    
    
}

void Ultra::OutputRepeats(bool flush) {
    
    
    // unsigned long long symbolsMasked = 0;
    // unsigned long long lastSeq = 0;
    
    
    int min = 10 * numberOfThreads;
    
    if (flush) {
        for (int i = 0; i < numberOfThreads; ++i) {
            if (threads[i]->repeats.size() > 0) {
                outRepeats.insert(outRepeats.end(), threads[i]->repeats.begin(), threads[i]->repeats.end());
            }
        }
        
        min = 0;
    }
    
    SortRepeatRegions();
    
    for (int i = (int)outRepeats.size() - 1; i >= min; --i) {
        
        RepeatRegion *next = NULL;
        if (i > 0)
            next = outRepeats[i - 1];
        
        RepeatRegion *r = outRepeats[i];
        outRepeats.pop_back();
        
        if (r->regionScore < scoreThreshold) {
            delete r;
            continue;
        }
        
        fprintf(out, "%s,", r->sequenceName.c_str());
        if (outputReadID) {
            fprintf(out, "%llu,", r->readID);
        }
        fprintf(out, "%llu,", r->sequenceStart);
        fprintf(out, "%llu,", r->sequenceStart + r->repeatLength);
        
        fprintf(out, "%i,", r->repeatPeriod);
        fprintf(out, "%f,", r->regionScore);
        
        fprintf(out, "%i,%i,%i,", r->mismatches, r->insertions, r->deletions);
        fprintf(out, "%s,", r->GetConsensus().c_str());
        
        if (outputRepeatSequence) {
            fprintf(out, "%s", r->sequence.c_str());
        }
        
        fprintf(out, "\n");
        
        delete r;
    }
    
    fflush(out);
}

void Ultra::AddRepeat(RepeatRegion *region) {
    /*
    pthread_mutex_lock(&repeatLock);
    bool output = false;
    repeats.push_back(region);
    
    
    if (repeatBuffer < repeats.size()) {
        {
            for (int i = 0; i < repeats.size(); ++i) {
                outRepeats.push_back(repeats[i]);
            }
            
            outRepeats.clear();
            
            output = true;
            
        }
    }
    
    pthread_mutex_unlock(&repeatLock);
    
    if (output)
        OutputRepeats();*/
    
}

void Ultra::SortRepeatRegions() {
    std::sort(outRepeats.begin(), outRepeats.end(), CompareRepeatOrder());
}

Ultra::Ultra(UModel*        m,
             FASTAReader*   r)
{
    out = stdout;
    
    reader = r;
    outRepeats.reserve(repeatBuffer + 100);
}

Ultra::Ultra(Settings* s, int n) {
    settings = s;
    
    if (settings->v_outFilePath == "")
        out = stdout;
    
    else
        out = fopen(settings->v_outFilePath.c_str(), "w");
    
    numberOfThreads = settings->v_numberOfThreads;
    
    if (settings->v_readWholeFile)
        settings->v_windowSize = -1;
    
    scoreThreshold = settings->v_scoreThreshold;
    outputReadID = settings->v_showWindowID;
    outputRepeatSequence = settings->v_outputRepeatSequence;
    
    
    //printf("Creating reader.\n");
    reader = new FASTAReader(settings->v_filePath,
                             settings->v_numberOfWindows,
                             settings->v_windowSize,
                             settings->v_overlapSize);
    
    int leng = settings->v_windowSize + settings->v_overlapSize + 2;
    
    //printf("Creating threads.\n");
    for (int i = 0; i < numberOfThreads; ++i) {
        //printf("Creating model %i\n", i);
        UModel *mod = new UModel(settings->v_maxPeriod,
                                 settings->v_maxInsertion,
                                 settings->v_maxDeletion,
                                 leng);
        
        mod->SetMatchProbabilities(settings->v_matchProbability);
        
        mod->tp_zeroToMatch = settings->v_zeroToMatch;
        mod->tp_matchToZero = settings->v_matchToZero;
        
        mod->tp_matchToInsertion  = settings->v_matchToInsertion;
        mod->tp_matchToDeletion   = settings->v_matchToDeletion;
        
        mod->tp_consecutiveInsertion  = settings->v_consecutiveInsertion;
        mod->tp_consecutiveDeletion   = settings->v_consecutiveDeletion;
        
        mod->CalculateScores();
        
        models.push_back(mod);
    }
    
    
    for (int i = 0; i < numberOfThreads; ++i) {
        uthread *newThread = new uthread;
        newThread->id = i;
        newThread->ultra = this;
        newThread->model = models[i];
        newThread->repeats.reserve(64);
        threads.push_back(newThread);
    }
    
    
    if (numberOfThreads == 1) {
        reader->multithread = false;
        multithreading = false;
    }
    
    else {
        reader->multithread = true;
        multithreading = true;
    }
    
}


bool CompareRepeatOrder::operator() (RepeatRegion *lhs, RepeatRegion *rhs) {
   // printf("%i vs %i and %i vs %i\n", lhs->readID, rhs->readID, lhs->windowStart, rhs->windowStart);
    if (lhs->readID != rhs->readID) {
        return lhs->readID > rhs->readID;
    }
    
    return lhs->windowStart > rhs->windowStart;
}

