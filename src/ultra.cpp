//
//  ultra.cpp
//  ultraP
//

#include <cmath>

#include "ultra.hpp"

void *UltraThreadLaunch(void *dat) {
  uthread *uth = (uthread *)dat;

  uth->ultra->AnalyzeFileWithThread(dat);
  return nullptr;
}

void Ultra::AnalyzeFile() {

  if (numberOfThreads <= 0) {
    printf("Number of threads must be set greater than 0.\n");
    exit(-1);
  }

  reader->FillWindows();

  if (pthread_mutex_init(&outerLock, nullptr) != 0) {
    printf("Failed to create outer mutex lock. Exiting.\n");
    exit(-1);
  }

  if (pthread_mutex_init(&innerLock, nullptr) != 0) {
    printf("Failed to create inner mutex lock. Exiting.\n");
    exit(-1);
  }

  if (pthread_mutex_init(&repeatLock, nullptr) != 0) {
    printf("Failed to create repeat lock. Exiting.\n");
    exit(-1);
  }


  OutputULTRASettings();
  InitializeWriter();

 // double Log2PvalForScore(double score, double period);


  for (int i = 1; i < numberOfThreads; ++i) {
    pthread_create(&threads[i]->p_thread, NULL, UltraThreadLaunch, threads[i]);
  }

  UltraThreadLaunch(threads[0]);
}

void Ultra::InitializeWriter() {
  writer->InitializeWriter(this);
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

double Ultra::Log2PvalForScore(double score, double period) {
  double loc = (settings->v_exponLocM * period) + settings->v_exponLocB;
  double scale = (settings->v_exponScaleM * period) + settings->v_exponScaleB;

  // Cap location
  if (loc < 0.2)
    loc = 0.2;

  return (-1.0 * (score - loc) / scale) / log2(2.71828);
}

void Ultra::AnalyzeSequenceWindow(SequenceWindow *sequence, uthread *uth) {

  int sleng = (int)sequence->length + (int)sequence->overlap;

  UModel *model = uth->model;
  model->matrix->RestartMatrix();

  for (int i = 0; i < sleng; ++i) {
    model->CalculateCurrentColumn(sequence, i);
  }

  model->matrix->CalculateTraceback(model->matrix->previousColumnIndex);
  int i = 0;
  RepeatRegion *r = GetNextRepeat(sequence, model->matrix, &i);

  while (r != NULL) {

    // Calculate P val
    r->logPVal = Log2PvalForScore(r->regionScore, r->repeatPeriod);

    if (storeSequence) {
      r->StoreSequence(sequence);
    }

    //if (storeTraceback) {
    //  r->storeTraceback(uth->model->matrix);
    //}

    if (storeScores) {
      r->StoreScores(model->matrix);
    }

    if (storeProfileNumbers) {
      r->GetLogoNumbers();
    }

    uth->repeats.push_back(r);

    r = GetNextRepeat(sequence, model->matrix, &i);
  }

  if (primaryThread == uth->id) {

    outRepeats.insert(outRepeats.end(), uth->repeats.begin(),
                      uth->repeats.end());
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

  double cScorePerSymbol =
      c->regionScore / (double)(c->repeatLength - c->repeatPeriod);
  double nScorePerSymbol =
      n->regionScore / (double)(n->repeatLength - n->repeatPeriod);

  if (c->sequenceID == n->sequenceID) {

    if (c->sequenceStart + c->repeatLength >
        n->sequenceStart + n->repeatLength) {
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

        // if (n->se)

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

      // If C starts in the window overlap we can always let n take care
      // of it
      if (c->windowStart >= reader->maxSeqLength - n->winOverlapSize) {
        delete c;
        n->overlapCorrection = OC_OVERLAP_RIGHT;
        return true;
      }

      // If C starts before window overlap and n ends before window
      // overlap, we can let c take care of it
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
          n->regionScore =
              (double)(c->repeatLength - c->repeatPeriod) * cScorePerSymbol;
          n->regionScore +=
              (double)(n->repeatLength - n->repeatPeriod - overlap) *
              nScorePerSymbol;

          if (outputRepeatSequence) {
            n->sequence = c->sequence.substr(0, c->repeatLength - overlap)
                              .append(n->sequence);
          }

          if (settings->v_showTraceback) {
            n->traceback = c->traceback.substr(0, c->repeatLength - overlap)
                               .append(n->traceback);
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
            c->sequence =
                c->sequence.substr(0, c->repeatLength - (overlap / 2));
          }

          if (settings->v_showTraceback) {
            c->traceback =
                c->traceback.substr(0, c->repeatLength - (overlap / 2));
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

  /*if (AnalyzingJSON) {
    if (!flush)
      return;

    else {
      OutputJSONRepeats();
      return;
    }
  }*/
  // unsigned long symbolsMasked = 0;
  // unsigned long lastSeq = 0;


  int min = numberOfThreads * numberOfThreads;

  if (flush) {
    for (int i = 0; i < numberOfThreads; ++i) {
      if (threads[i]->repeats.size() > 0) {
        outRepeats.insert(outRepeats.end(), threads[i]->repeats.begin(),
                          threads[i]->repeats.end());
      }
    }

    min = 0;
  }

  SortRepeatRegions();

  while (outRepeats.size() > min) {

   // if (settings->v_correctOverlap)
   //   while (FixRepeatOverlap() && outRepeats.size() > min)
   //     ;

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
    delete r;
    r = NULL;
  }

  if (flush) {
    writer->EndWriter();
    //fprintf(out, "]\n}\n");
  }

  fflush(out);
}

void Ultra::OutputJSONRepeats() {
  for (int i = 0; i < reader->jsonReader->repeats.size(); ++i) {

    reader->jsonReader->repeats[i]->OutputRepeat(out, scoreThreshold, i == 0);
  }

  fprintf(out, "]\n}\n");
}

void Ultra::OutputULTRASettings() {
  fprintf(settings_out, "{\"Version\": \"%s\", \n", settings->StringVersion().c_str());
  fprintf(settings_out, "\"Parameters\": {\n");
  fprintf(settings_out, "%s}}\n", settings->JSONString().c_str());
  if (settings_out != stdout)
    fclose(settings_out);
}

void Ultra::OutputJSONStart() {
  /*
  fprintf(out, "{\"Passes\":[\n");

  if (reader->format == JSON) {
    if (reader->jsonReader->passes.size() > 0) {

      for (int i = 0; i < reader->jsonReader->passes.size(); ++i) {
        reader->jsonReader->passes[i]->OutputPass(out);
        fprintf(out, ",\n");
      }
    }
  }

  fprintf(out, "{");
  fprintf(out, "\"Pass ID\": %i,\n", passID);
  fprintf(out, "\"Version\": \"%s\",\n", settings->StringVersion().c_str());
  fprintf(out, "\"Parameters\": {\n");
  fprintf(out, "%s\n}\n}],\n", settings->JSONString().c_str());
*/
  fprintf(out, "{\"Repeats\": [");
}

void Ultra::OutputJSONKey(std::string key) {
  fprintf(out, "\"%s\": ", key.c_str());
}

void Ultra::OutputRepeat(RepeatRegion *r, bool isSubRep) {
  writer->WriteRepeat(r);
  /*
  if (!firstRepeat && !isSubRep) {
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

  if (settings->v_calculateLogPVal) {
    fprintf(out, ",\n");
    OutputJSONKey("Log2 Pval");
    fprintf(out, "%f", r->logPVal);
  }

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

  if (settings->v_showLogoNumbers) {
    fprintf(out, ",\n");
    OutputJSONKey("Logo Numbers");

    fprintf(out, "\"");

    fprintf(out, "%i", r->logoNumbers[0]);
    for (int i = 1; i < r->repeatLength; ++i) {
      fprintf(out, ",%i", r->logoNumbers[i]);
    }

    fprintf(out, "\"");
  }

  if (settings->v_showScores) {
    fprintf(out, ",\n");
    OutputJSONKey("PositionScoreDelta");
    //  fprintf(out, "\"%s\"", r->traceback.c_str());
    std::string scores = "";

    for (unsigned long i = 0; i < r->repeatLength; ++i) {
      if (i > 0)
        scores += ":";
      scores += std::to_string(r->scores[i]);
    }

    fprintf(out, "\"%s\"", scores.c_str());
  }

  if (settings->v_debugOverlapCorrection) {
    fprintf(out, ",\n");
    OutputJSONKey("OC");
    fprintf(out, "\"%i\"", r->overlapCorrection);
  }

  if (!isSubRep && settings->v_splitRepeats &&
      r->repeatPeriod < settings->v_maxSplitPeriod) {
    std::vector<RepeatRegion *> *subReps =
        r->SplitRepeats(settings->v_splitDepth, settings->v_splitCutoff);

    if (subReps != NULL) {
      fprintf(out, ",\n");
      OutputJSONKey("Subrepeats");
      fprintf(out, "[");
      for (int i = 0; i < subReps->size(); ++i) {
        if (i > 0)
          fprintf(out, ",\n");
        subReps->at(i)->logPVal = Log2PvalForScore(
            subReps->at(i)->regionScore, subReps->at(i)->repeatPeriod);
        OutputRepeat(subReps->at(i), true);
        delete subReps->at(i);
      }
      fprintf(out, "]");

      subReps->clear();
      delete subReps;
    }
  }

  fprintf(out, "}");*/

}

void Ultra::SortRepeatRegions() {
  std::sort(outRepeats.begin(), outRepeats.end(), CompareRepeatOrder());
}

Ultra::Ultra(Settings *s, int n) {
  settings = s;

  if (settings->v_outFilePath == "") {
    out = stdout;
    settings_out = stdout;
  }

  else {
    out = fopen(settings->v_outFilePath.c_str(), "w");
    std::string settings_file_path = settings->v_outFilePath + ".settings";
    settings_out = fopen(settings_file_path.c_str(), "w");
  }

  if (settings->v_outputFormat == JSON) {
    writer = new JSONFileWriter();
  }

  else if (settings->v_outputFormat == BED) {
    writer = new BEDFileWriter();
  }


  numberOfThreads = settings->v_numberOfThreads;

  if (settings->v_readWholeFile)
    settings->v_windowSize = -1;

  scoreThreshold = settings->v_scoreThreshold;
  outputReadID = settings->v_showWindowID;
  outputRepeatSequence = settings->v_outputRepeatSequence;

  passID = settings->v_passID;
  AnalyzingJSON = settings->v_JSONInput;

  if (settings->v_JSONInput) {
    reader = new FileReader(settings->v_filePath, settings->v_windowSize,
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
      passID = (int)largestPass;
    }
  }

  else {
    reader = new FileReader(settings->v_filePath, settings->v_numberOfWindows,
                            settings->v_windowSize, settings->v_overlapSize,
                            settings->v_numberOfThreads > 1);
  }

  if (passID < 0)
    passID = 0;

  int leng = settings->v_windowSize + settings->v_overlapSize + 2;

  if (settings->v_splitRepeats) {
    storeTraceAndSequence = true;
    printf("STORE TRACE AND SEQUENCE\n");
  }

  if (settings->v_showLogoNumbers) {
    storeTraceAndSequence = true;
  }

  // printf("Creating threads.\n");
  for (int i = 0; i < numberOfThreads; ++i) {
    // We now are making the v_maxPeriod setting more intuitive, by adding 1 to
    // it. This makes a v_maxPeriod of 10 able to detect repeats of length 10.
    UModel *mod =
        new UModel(settings->v_maxPeriod + 1, settings->v_maxInsertion,
                   settings->v_maxDeletion, leng);

    mod->periodDecay = settings->v_repeatPeriodDecay;

    mod->SetMatchProbabilities(settings->v_matchProbability);
    mod->SetATCGProbabilities(settings->v_Apctg, settings->v_Tpctg,
                              settings->v_Cpctg, settings->v_Gpctg);

    // printf(mod->backgroundProbabilties)

    mod->tp_zeroToMatch = settings->v_zeroToMatch;
    mod->tp_matchToZero = settings->v_matchToZero;

    mod->tp_matchToInsertion = settings->v_matchToInsertion;
    mod->tp_matchToDeletion = settings->v_matchToDeletion;

    mod->tp_consecutiveInsertion = settings->v_consecutiveInsertion;
    mod->tp_consecutiveDeletion = settings->v_consecutiveDeletion;

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

bool CompareRepeatOrder::operator()(RepeatRegion *lhs, RepeatRegion *rhs) {
  if (lhs->readID != rhs->readID) {
    return lhs->readID > rhs->readID;
  }

  return lhs->windowStart > rhs->windowStart;
}
