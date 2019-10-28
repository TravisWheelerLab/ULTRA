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
    
    OutputJSONStart();
    
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
    
    if (!reader->DoneReadingFile()) {
        if (!reader->IsReading()) {
            if (reader->ReadyWindowsSize() < minReaderSize) {
                shouldRead = true;
                reader->SetIsReading(true);
                
            }
        }
    }
    
    if (multithreading)
        pthread_mutex_unlock(&outerLock);
    
    if (shouldRead) {
        if (multithreading)
            pthread_mutex_lock(&innerLock);
        
        reader->FillWindows();
        reader->SetIsReading(false);
        
        if (multithreading)
            pthread_mutex_unlock(&innerLock);
    }
    
    retval = reader->GetReadyWindow();
    
    return retval;
}

void Ultra::AnalyzeFileWithThread(void *dat) {
    
    uthread *uth = (uthread *)dat;
    int tid = uth->id;
    
    
    SequenceWindow *currentWindow = GetSequenceWindow(NULL);
    while (currentWindow != NULL || !reader->DoneReadingFile()) {
        //printf("%i: %llx %i\n", i++, (unsigned long)currentWindow, reader->doneReadingFile);
        
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
        
        
        if (sequence->jrepeat) {
            JSONRepeat *jrep = new JSONRepeat();
            
            jrep->passID = passID;
            jrep->start = r->sequenceStart;
            jrep->length = r->repeatLength;
            jrep->consensus = r->GetConsensus();
            
            if (settings->v_outputRepeatSequence) {
                r->StoreSequence(sequence);
                jrep->sequence = r->sequence;
            }
            
            if (settings->v_showTraceback) {
                r->StoreTraceback(model->matrix);
                jrep->traceback = r->traceback;
            }
            
            jrep->insertions = r->insertions;
            jrep->deletions = r->deletions;
            jrep->substitutions = r->mismatches;
            jrep->score = r->regionScore;
            
            sequence->jrepeat->subrepeats.push_back(jrep);
            
        }
        
        if (outputRepeatSequence) {
            r->StoreSequence(sequence);
        }
        
        if (settings->v_showTraceback) {
            r->StoreTraceback(model->matrix);
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

bool Ultra::FixRepeatOverlap() {
    // We will return true if we combined repeats
    // return false;
    
    RepeatRegion *c = outRepeats.back();
    outRepeats.pop_back();
    
    if (outRepeats.empty()) {
        outRepeats.push_back(c);
        return false;
    }
    
    RepeatRegion *n = outRepeats.back();
    
    unsigned long cSeqEnd = c->sequenceStart + c->repeatLength;
    // unsigned long nSeqEnd = n->sequenceStart + n->repeatLength;
    
    double cScorePerSymbol = c->regionScore / (double)(c->repeatLength - c->repeatPeriod);
    double nScorePerSymbol = n->regionScore / (double)(n->repeatLength - n->repeatPeriod);
    
    if (c->sequenceID == n->sequenceID) {
        
        //Check to see if n is entirely contained in c
        
        if (c->sequenceStart + c->repeatLength > n->sequenceStart + n->repeatLength) {
            delete n;
            outRepeats.pop_back();
            c->overlapCorrection = OC_PERIOD_LEFT;
            outRepeats.push_back(c);
            return true;
        }
        
        // Deal with small overlaps
        if (c->readID == n->readID) {
            // Just make sure that n doesn't dwarf c
            if (n->sequenceStart <= c->sequenceStart + 3) {
                delete c;
                n->overlapCorrection = OC_PERIOD_RIGHT;
                return true;
            }
            
            // Check to see if there is partial overlap between c and n
            
            else if (n->sequenceStart < cSeqEnd) {
                
                //if (n->se)
                
                n->overlapCorrection = OC_PERIOD_SPLIT;
                
                unsigned long overlap = (cSeqEnd - n->sequenceStart);
                n->sequenceStart += overlap;
                n->repeatLength -= overlap;
                
                if (outputRepeatSequence) {
                    n->sequence = n->sequence.substr(overlap);
                }
                
                if (settings->v_showTraceback) {
                    n->traceback = n->traceback.substr(overlap);
                }
            }
            
            
        }
        
        // Check for overlap problems across the sequence window overlap
        else if (n->readID == c->readID + 1 && cSeqEnd > n->sequenceStart) {
            
            unsigned long overlap = (cSeqEnd - n->sequenceStart);
            
            
            // If C starts in the window overlap we can always let n take care of it
            if (c->windowStart >= reader->maxSeqLength - n->winOverlapSize) {
                delete c;
                n->overlapCorrection = OC_OVERLAP_RIGHT;
                return true;
            }
            
            // If C starts before window overlap and n ends before window overlap, we can let c take care of it
            else if (n->windowStart + n->repeatLength < n->winOverlapSize) {
                outRepeats.pop_back();
                delete n;
                c->overlapCorrection = OC_OVERLAP_LEFT;
                outRepeats.push_back(c);
                return false;
            }
            
            // The repeat goes across the repeat window, we'll split it up
            else {
                
                // Check to see if we can just join the repeats
                if (c->repeatPeriod == n->repeatPeriod) {
                    n->sequenceStart = c->sequenceStart;
                    
                    n->combinedRepeat = true;
                    n->regionScore = (double)(c->repeatLength - c->repeatPeriod) * cScorePerSymbol;
                    n->regionScore += (double)(n->repeatLength - n->repeatPeriod - overlap) * nScorePerSymbol;
                    
                    if (outputRepeatSequence) {
                        n->sequence = c->sequence.substr(0, c->repeatLength - overlap).append(n->sequence);
                    }
                    
                    if (settings->v_showTraceback) {
                        n->traceback = c->traceback.substr(0, c->repeatLength - overlap).append(n->traceback);
                    }
                    
                    n->overlapCorrection = OC_OVERLAP_RIGHT;
                    delete c;
                    return true;
                    
                }
                
                // Split things down the middle
                else {
                    c->repeatLength -= overlap / 2;
                    c->regionScore -= (double)(overlap / 2) * cScorePerSymbol;
                    
                    if (outputRepeatSequence) {
                        c->sequence = c->sequence.substr(0, c->repeatLength - (overlap / 2));
                    }
                    
                    if (settings->v_showTraceback) {
                        c->traceback = c->traceback.substr(0, c->repeatLength - (overlap / 2));
                    }
                    
                    overlap = overlap + 1;
                    
                    n->repeatLength -= (overlap / 2);
                    n->sequenceStart += (overlap / 2);
                    n->regionScore -= (double)(overlap / 2) * nScorePerSymbol;
                    
                    if (outputRepeatSequence) {
                        n->sequence = n->sequence.substr(overlap / 2);
                    }
                    
                    if (settings->v_showTraceback) {
                        n->traceback = n->traceback.substr(overlap / 2);
                    }
                    
                    c->overlapCorrection = OC_OVERLAP_SPLIT;
                    
                    
                }
                
            }
            
        }
    }
    
    outRepeats.push_back(c);
    return false;
}

void Ultra::OutputRepeats(bool flush) {
    
    if (AnalyzingJSON) {
        if (!flush)
            return;
        
        else {
            OutputJSONRepeats();
            return;
        }
    }
    // unsigned long symbolsMasked = 0;
    // unsigned long lastSeq = 0;
    
    
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
    
    
    while (outRepeats.size() > min) {
        
        
        while (FixRepeatOverlap() && outRepeats.size() > min);
        
        RepeatRegion *r = outRepeats.back();
        outRepeats.pop_back();
        
        
        if (r->regionScore < scoreThreshold ||
            r->repeatLength < (r->repeatPeriod * settings->v_repeatThreshold) ||
            r->repeatLength < settings->v_lengthThreshold) {
            
            delete r;
            r = NULL;
            continue;
        }
        
        OutputRepeat(r);
        //fprintf(out, "%f to %f, ", r->startScore, r->endScore);
        /*fprintf(out, "%s", r->sequenceName.c_str());
         
         if (outputReadID) {
         fprintf(out, ",%lu", r->readID);
         }
         
         if (settings->v_debugOverlapCorrection) {
         fprintf(out, ",%i", r->overlapCorrection);
         }
         
         fprintf(out, ",%lu", r->sequenceStart);
         fprintf(out, ",%lu", r->sequenceStart + r->repeatLength);
         
         fprintf(out, ",%i", r->repeatPeriod);
         fprintf(out, ",%f", r->regionScore);
         
         fprintf(out, ",%i,%i,%i", r->mismatches, r->insertions, r->deletions);
         fprintf(out, ",%s", r->GetConsensus().c_str());
         
         if (outputRepeatSequence) {
         fprintf(out, ",%s", r->sequence.c_str());
         }
         
         if (settings->v_showTraceback) {
         fprintf(out, ",%s", r->traceback.c_str());
         }
         
         fprintf(out, "\n");*/
        
        delete r;
        r = NULL;
    }
    
    if (flush) {
        fprintf(out, "]\n}\n");
    }
    
    fflush(out);
}

void Ultra::OutputJSONRepeats() {
    for (int i = 0; i < reader->jsonReader->repeats.size(); ++i) {
        reader->jsonReader->repeats[i]->OutputRepeat(out, i == 0);
    }
    
    fprintf(out, "]\n}\n");
}



void Ultra::OutputJSONStart() {
    fprintf(out, "{\"Passes\":[{\n");
    fprintf(out, "\"Pass ID\": %i,\n", passID);
    fprintf(out, "\"Version\": \"%s\",\n", ULTRA_VERSION_STRING);
    fprintf(out, "\"Parameters\": {\n");
    fprintf(out, "%s\n}\n}],\n", settings->JSONString().c_str());
    
    fprintf(out, "\"Repeats\": [");
}

void Ultra::OutputJSONKey(std::string key) {
    fprintf(out, "\"%s\": ", key.c_str());
}

void Ultra::OutputRepeat(RepeatRegion *r) {
    if (!firstRepeat) {
        fprintf(out, ",\n\n");
    }
    
    else {
        firstRepeat = false;
    }
    
    fprintf(out, "{");
    
    OutputJSONKey("PassID");
    fprintf(out, "%i", passID);
    
    fprintf(out, ",\n");
    OutputJSONKey("SequenceName");
    fprintf(out, "\"%s\"", r->sequenceName.c_str());
    
    fprintf(out, ",\n");
    OutputJSONKey("Start");
    fprintf(out, "%lu", r->sequenceStart);
    
    fprintf(out, ",\n");
    OutputJSONKey("Length");
    fprintf(out, "%lu", r->repeatLength);
    
    fprintf(out, ",\n");
    OutputJSONKey("Period");
    fprintf(out, "%i", r->repeatPeriod);
    
    fprintf(out, ",\n");
    OutputJSONKey("Score");
    fprintf(out, "%f", r->regionScore);
    
    fprintf(out, ",\n");
    OutputJSONKey("Substitutions");
    fprintf(out, "%i", r->mismatches);
    
    fprintf(out, ",\n");
    OutputJSONKey("Insertions");
    fprintf(out, "%i", r->insertions);
    
    fprintf(out, ",\n");
    OutputJSONKey("Deletions");
    fprintf(out, "%i", r->deletions);
    
    
    
    fprintf(out, ",\n");
    OutputJSONKey("Consensus");
    fprintf(out, "\"%s\"", r->GetConsensus().c_str());

    if (outputReadID) {
        fprintf(out, ",\n");
        OutputJSONKey("ReadID");
        fprintf(out, "\"%lu\"", r->readID);
        
    }
    
    if (outputRepeatSequence) {
        
        fprintf(out, ",\n");
        OutputJSONKey("Sequence");
        fprintf(out, "\"%s\"", r->sequence.c_str());
        
    }
    
    if (settings->v_showTraceback) {
        fprintf(out, ",\n");
        OutputJSONKey("Traceback");
        fprintf(out, "\"%s\"", r->traceback.c_str());
        
    }
    
    if (settings->v_debugOverlapCorrection) {
        fprintf(out, ",\n");
        OutputJSONKey("OC");
        fprintf(out, "\"%i\"", r->overlapCorrection);
        
    }
    
    fprintf(out, "}");
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
    
    passID = settings->v_passID;
    AnalyzingJSON = settings->v_JSONInput;
    
    if (settings->v_JSONInput) {
        reader = new FileReader(settings->v_filePath,
                                 settings->v_windowSize,
                                 settings->v_overlapSize,
                                settings->v_numberOfThreads > 1);
        
        reader->jsonReader->ReadFile();
        
        if (passID == -1) {
            unsigned long largestPass = 1;
            for (int i = 0; i < reader->jsonReader->passes.size(); ++i) {
                int t = reader->jsonReader->passes[i]->passID;
                
                if (t >= largestPass)
                    largestPass = t + 1;
            }
            passID = largestPass;
        }
        
    }
    
    
    else {
        reader = new FileReader(settings->v_filePath,
                                settings->v_numberOfWindows,
                                 settings->v_windowSize,
                                 settings->v_overlapSize,
                                settings->v_numberOfThreads > 1);
    }
    
    if (passID < 0)
    passID = 0;
   
    int leng = settings->v_windowSize + settings->v_overlapSize + 2;
    
    //printf("Creating threads.\n");
    for (int i = 0; i < numberOfThreads; ++i) {
        //printf("Creating model %i\n", i);
        UModel *mod = new UModel(settings->v_maxPeriod,
                                 settings->v_maxInsertion,
                                 settings->v_maxDeletion,
                                 leng);
        
        mod->periodDecay = settings->v_repeatPeriodDecay;
        
        mod->SetMatchProbabilities(settings->v_matchProbability);
        
        mod->tp_zeroToMatch = settings->v_zeroToMatch;
        mod->tp_matchToZero = settings->v_matchToZero;
        
        mod->tp_matchToInsertion  = settings->v_matchToInsertion;
        mod->tp_matchToDeletion   = settings->v_matchToDeletion;
        
        mod->tp_consecutiveInsertion  = settings->v_consecutiveInsertion;
        mod->tp_consecutiveDeletion   = settings->v_consecutiveDeletion;
        
        mod->CalculateScores();
        /*
         for (int i = 0; i < 5; ++i) {
         for (int j = 0; j < 5; ++j) {
         printf("%f ", mod->mscore[i][j]);
         }
         printf("\n");
         }*/
        
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

