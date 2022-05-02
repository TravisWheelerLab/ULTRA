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
    pthread_create(&threads[i]->p_thread, nullptr, UltraThreadLaunch,
                   threads[i]);
  }

  UltraThreadLaunch(threads[0]);
}

void Ultra::InitializeWriter() { writer->InitializeWriter(this); }

SequenceWindow *Ultra::GetSequenceWindow(SequenceWindow *seq) {

  SequenceWindow *retval = nullptr;
  bool shouldRead = false;

  if (seq != nullptr)
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

// NOTE:
// THIS IS NOT THREAD SAFE
// WE ASSUME THATS OK (FOR NOW)
int Ultra::SmallestReadID() {
  int smallest = 100000000;

  for (int i = 0; i < threads.size(); ++i) {
    if (threads[i]->activeReadID < smallest) {
      smallest = threads[i]->activeReadID;
    }
  }

  return smallest;
}

void Ultra::AnalyzeFileWithThread(void *dat) {

  uthread *uth = (uthread *)dat;
  int tid = uth->id;

  SequenceWindow *currentWindow = GetSequenceWindow(nullptr);
  while (currentWindow != nullptr || !reader->DoneReadingFile()) {
    if (currentWindow != nullptr)
      AnalyzeSequenceWindow(currentWindow, uth);

    currentWindow = GetSequenceWindow(currentWindow);
  }

  if (multithreading) {
    if (tid == 0) {
      for (int i = 1; i < numberOfThreads; ++i) {
        pthread_join(threads[i]->p_thread, nullptr);
      }
    }

    else {
      pthread_exit(nullptr);
    }
  }
}

double Ultra::Log2PvalForScore(double score, double period) const {
  double loc = (settings->v_exponLocM * period) + settings->v_exponLocB;
  double scale = (settings->v_exponScaleM * period) + settings->v_exponScaleB;

  // Cap location
  if (loc < 0.2)
    loc = 0.2;

  return (-1.0 * (score - loc) / scale) / log2(2.71828);
}

void Ultra::AnalyzeSequenceWindow(SequenceWindow *sequence, uthread *uth) {

  int sleng = (int)sequence->length + (int)sequence->overlap;

  if (uth->repeats.size() == 0) {
    uth->activeReadID = sequence->readID;
  }

  UModel *model = uth->model;
  model->matrix->RestartMatrix();

  for (int i = 0; i < sleng; ++i) {
    model->CalculateCurrentColumn(sequence, i);
  }

  // WE'RE RIGHT HERE
  // GOING TO TRY TO PUSH ALL THE CODE IN !!!

  model->matrix->CalculateTraceback(model->matrix->previousColumnIndex);
  int i = 0;
  RepeatRegion *r = GetNextRepeat(sequence, model->matrix, &i);

  while (r != nullptr) {

    // Calculate P val
    r->logPVal = Log2PvalForScore(r->regionScore, r->repeatPeriod);

    /*
    if (storeSequence) {
      r->StoreSequence(sequence);
    }
    if (storeScores) {
      r->StoreScores(model->matrix);
    }
    if (storeProfileNumbers) {
      r->GetLogoNumbers();
    }*/

    r->StoreSequence(sequence);
    r->StoreScores(model->matrix);
    r->GetLogoNumbers();

    uth->repeats.push_back(r);

    r = GetNextRepeat(sequence, model->matrix, &i);
  }

  if (primaryThread == uth->id) {

    uth->activeReadID = sequence->readID;
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

  int maxReadID = SmallestReadID() - 20;

  if (flush) {
    for (int i = 0; i < numberOfThreads; ++i) {
      if (threads[i]->repeats.size() > 0) {
        outRepeats.insert(outRepeats.end(), threads[i]->repeats.begin(),
                          threads[i]->repeats.end());
      }
    }

    maxReadID = 100000000;
  }

  SortRepeatRegions();

  while (outRepeats.size() > 0) {

    RepeatRegion *r = outRepeats.back();

    if (r->readID >= maxReadID) {
      break;
    }

    outRepeats.pop_back();
    OutputRepeat(r);
    std::vector<int> *splits = SplitRepeat(r, 0.2, 5, 5, 3);
    if (splits != nullptr) {
      printf("Splits!=null\n");
      for (int i = 0; i < splits->size(); ++i) {
          printf("split %i\n", splits->at(i));
        }
    }
    delete r;
    r = nullptr;
  }

  if (flush) {
    writer->EndWriter();
    // fprintf(out, "]\n}\n");
  }

  fflush(out);
}

void Ultra::OutputULTRASettings() {
  fprintf(settings_out, "{\"Version\": \"%s\", \n",
          settings->StringVersion().c_str());
  fprintf(settings_out, "\"Parameters\": {\n");
  fprintf(settings_out, "%s}}\n", settings->JSONString().c_str());
  if (settings_out != stdout)
    fclose(settings_out);
}

void Ultra::OutputRepeat(RepeatRegion *r, bool isSubRep) {
  writer->WriteRepeat(r);
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
