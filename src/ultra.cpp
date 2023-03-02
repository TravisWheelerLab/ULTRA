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
  // printf("Calculating traceback\n");
  model->matrix->CalculateTraceback(model->matrix->previousColumnIndex);
  int i = 0;
  RepeatRegion *r = GetNextRepeat(sequence, model->matrix, &i);

  while (r != nullptr) {

    // Calculate P val

    r->logPVal = Log2PvalForScore(r->regionScore, r->repeatPeriod);

    if (storeTraceAndSequence) {
      // printf("Storing trace...\n");
      r->StoreSequence(sequence);
      r->StoreTraceback(model->matrix);
      r->StoreScores(model->matrix);

      r->GetLogoNumbers();
    }
    int splitWindow = std::max(r->repeatPeriod * settings->v_splitDepth,
                               settings->v_minSplitWindow);
    if (r->repeatPeriod <= settings->v_maxSplitPeriod) {
      if (r->repeatLength > 2 * splitWindow) {

        // This almost feels right.
        // Probably want to change it later though.
        float join_threshold = 1.0 - (1.0 / settings->v_splitThreshold);

        r->splits = uth->splitter->SplitsForRegion(r, splitWindow,
                                                   settings->v_splitThreshold);
        r->consensi = uth->splitter->ConsensiForSplit(r, r->splits, 0.65);
        ValidateSplits(r->consensi, r->splits, 0.85);
      } else {
        r->splits = new std::vector<int>;
        r->consensi = new std::vector<std::string>;
        r->consensi->push_back(r->string_consensus);
      }
    }

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

void Ultra::CorrectOverlap(int maxReadID) {}

void Ultra::OutputRepeats(bool flush) {

  int maxReadID = SmallestReadID() - (2 * numberOfThreads);

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

  while (!outRepeats.empty()) {

    RepeatRegion *r = outRepeats.back();

    if (settings->v_correctOverlap) {
      CorrectOverlap(maxReadID);
    }

    if (r->readID >= maxReadID) {
      break;
    }

    outRepeats.pop_back();

    // Check if we need to correct overlap
    // Correct overlap while possible
    if (settings->v_correctOverlap) {
      while (!outRepeats.empty()) {
        if (outRepeats.back()->readID >= maxReadID) {
          outRepeats.push_back(r);
          return;
        }

        if (repeats_overlap(r, outRepeats.back())) {
          r = joint_repeat_region(r, outRepeats.back());
          if (r->splits != nullptr)
            ValidateSplits(r->consensi, r->splits, 0.85);
          outRepeats.pop_back();
        }

        else
          break;
      }
    }

    r->SortConsensi();
    OutputRepeat(r);

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
  if (!settings->v_suppress_out)
    writer->WriteRepeat(r);
  if (settings->v_maskWithCase ||
      settings->v_maskWithN) {

    this->StoreMaskForRegion(r);
  }
}

void Ultra::StoreMaskForRegion(RepeatRegion *r) {
  if (masks_for_seq.find(r->sequenceName) == masks_for_seq.end())
    masks_for_seq[r->sequenceName] = new std::vector<mregion>();
  masks_for_seq[r->sequenceName]->push_back(mregion{r->sequenceStart, r->repeatLength + r->sequenceStart});
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

  reader = new FileReader(settings->v_filePath, settings->v_numberOfWindows,
                          settings->v_windowSize, settings->v_overlapSize,
                          settings->v_numberOfThreads > 1);

  int leng = settings->v_windowSize + settings->v_overlapSize + 2;

  if (settings->v_maxSplitPeriod > 0) {
    storeTraceAndSequence = true;
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
    newThread->splitter = new SplitWindow();
    newThread->splitter->AllocateSplitWindow(4, settings->v_maxPeriod + 1);
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
