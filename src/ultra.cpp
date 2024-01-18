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

  // float Log2PvalForScore(float score, float period);

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

float Ultra::Log2PvalForScore(float score, float period) const {
  float loc =
      (settings->pval_exponent_loc_m * period) + settings->pval_exponent_loc_b;
  float scale = (settings->pval_exponent_scale_m * period) +
                settings->pval_exponent_scale_b;

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

  if (shuffleSequence) {
    ShuffleSequenceWindow(sequence);
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
    int splitWindow = std::max(r->repeatPeriod * settings->split_depth,
                               settings->min_split_window);
    if (r->repeatPeriod <= settings->max_split) {
      if (r->repeatLength > 2 * splitWindow) {

        // This almost feels right.
        // Probably want to change it later though.
        float join_threshold = 1.0 - (1.0 / settings->split_threshold);

        r->splits = uth->splitter->SplitsForRegion(r, splitWindow,
                                                   settings->split_threshold);
        r->consensi = uth->splitter->ConsensiForSplit(r, r->splits, 0.65);
        ValidateSplits(r->consensi, r->splits, 0.85);
      } else {
        r->splits = new std::vector<int>;
        r->consensi = new std::vector<std::string>;
        r->consensi->push_back(r->string_consensus);
      }
    }

    // Make sure we have at least 1 unit of repeat here...
    if (r->repeatLength >= r->repeatPeriod)
      uth->repeats.push_back(r);
    else {
      delete r;
      r = nullptr;
    }

    r = GetNextRepeat(sequence, model->matrix, &i);
  }

  if (primaryThread == uth->id) {

    uth->activeReadID = sequence->readID;
    outRepeats.insert(outRepeats.end(), uth->repeats.begin(),
                      uth->repeats.end());
    uth->repeats.clear();

    if (outRepeats.size() > repeatBuffer && !this->settings->disable_streaming_out) {
      OutputRepeats();
    }

    primaryThread = (primaryThread + 1) % numberOfThreads;
    // printf("%i\n", primaryThread);
  }
}

void Ultra::OutputRepeats(bool flush) {

  if (!flush && settings->disable_streaming_out)
    return;
  int maxReadID = SmallestReadID() - (3 * numberOfThreads);

  if (flush) {
    for (int i = 0; i < numberOfThreads; ++i) {
      if (threads[i]->repeats.size() > 0) {
        outRepeats.insert(outRepeats.end(), threads[i]->repeats.begin(),
                          threads[i]->repeats.end());
        threads[i]->repeats.clear();
      }
    }
    maxReadID = 100000000;
  }

  SortRepeatRegions();

  while (!outRepeats.empty()) {

    RepeatRegion *r = outRepeats.back();

    if (r->readID >= maxReadID) {
      break;
    }

    outRepeats.pop_back();

    // Check if we need to correct overlap
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

    r->SortConsensi();
    if (r->repeatLength / r->repeatPeriod >= settings->min_units)
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
  if (!settings->hide_settings) {
    fprintf(settings_out, "{\"Version\": \"%s\", \n", ULTRA_VERSION_STRING);
    fprintf(settings_out, "\"Parameters\": {\n");
    fprintf(settings_out, "%s}}\n", settings->json_string().c_str());
  }

  if (settings_out != stdout)
    fclose(settings_out);
}

void Ultra::OutputRepeat(RepeatRegion *r, bool isSubRep) {
  if (!settings->suppress_out)
    writer->WriteRepeat(r);
  if (settings->produce_mask) {
    this->StoreMaskForRegion(r);
  }
}

void Ultra::StoreMaskForRegion(RepeatRegion *r) {
  if (masks_for_seq.find(r->sequenceID) == masks_for_seq.end())
    masks_for_seq[r->sequenceID] = new std::vector<mregion>();
  masks_for_seq[r->sequenceID]->push_back(
      mregion{r->sequenceStart, r->repeatLength + r->sequenceStart});
}

void Ultra::SortRepeatRegions() {
  std::sort(outRepeats.begin(), outRepeats.end(), CompareRepeatOrder());
}

unsigned long long Ultra::Coverage() {
  unsigned long long coverage = 0;

  // Itterate through each sequence and find the coverage for that sequence
  for (auto [seq_id, regions] : this->masks_for_seq) {
    auto cleaned_regions = CleanedMasks(regions);

    for (auto region : *cleaned_regions) {
      coverage += region.end - region.start;
    }

    cleaned_regions->clear();
    delete cleaned_regions;

  }
  return coverage;
}

Ultra::Ultra(Settings *s) {
  settings = s;

  out = stdout;
  settings_out = stdout;

  if (!settings->out_file.empty()) {
    out = fopen(settings->out_file.c_str(), "w");
    std::string settings_file_path = settings->out_file + ".settings";
    if (!settings->hide_settings) {
      settings_out = fopen(settings_file_path.c_str(), "w");
    }
  }

  if (settings->json) {
    writer = new JSONFileWriter();
  }

  else {
    writer = new BEDFileWriter();
  }

  numberOfThreads = settings->threads;
  scoreThreshold = settings->min_score;
  outputReadID = settings->show_wid;
  outputRepeatSequence = !settings->hide_seq;

  passID = 0;

  reader = new FileReader(settings->in_file, settings->windows,
                          settings->window_size, settings->overlap,
                          settings->threads > 1);

  int leng = settings->window_size + (settings->overlap + 2);
  storeTraceAndSequence = true;
  if (settings->max_split > 0) {
    storeTraceAndSequence = true;
  }

  if (settings->show_logo_nums) {
    storeTraceAndSequence = true;
  }

  // printf("Creating threads.\n");
  for (int i = 0; i < numberOfThreads; ++i) {
    // We now are making the v_maxPeriod setting more intuitive, by adding 1 to
    // it. This makes a v_maxPeriod of 10 able to detect repeats of length 10.
    UModel *mod = new UModel(settings->max_period + 1, settings->max_insert,
                             settings->max_delete, leng);

    mod->periodDecay = settings->period_decay;

    mod->SetMatchProbabilities(settings->match_probability);
    mod->SetATCGProbabilities(settings->a_freq, settings->t_freq,
                              settings->c_freq, settings->g_freq);

    // printf(mod->backgroundProbabilties)

    mod->tp_zeroToMatch = settings->transition_nr;
    mod->tp_matchToZero = settings->transition_rn;

    mod->tp_matchToInsertion = settings->transition_ri;
    mod->tp_matchToDeletion = settings->transition_rd;

    mod->tp_consecutiveInsertion = settings->transition_ii;
    mod->tp_consecutiveDeletion = settings->transition_dd;

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
    if (!settings->no_split)
      newThread->splitter->AllocateSplitWindow(4, settings->max_period + 1);
    else
      newThread->splitter->maxPeriod = 0;
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

Ultra::~Ultra() {
  settings = nullptr;
  delete reader;
  reader = nullptr;

  delete writer;
  writer = nullptr;

  if (out != stdout && out != nullptr) {
    fclose(out);
    out = nullptr;
  }

  for (auto& pair : masks_for_seq) {
    delete pair.second;  // pair.second is a std::vector<mregion> *
  }
  masks_for_seq.clear();

  for (auto val : outRepeats) {
    delete val;  // pair.second is a std::vector<mregion> *
  }
  outRepeats.clear();

  for (auto val : models) {
    delete val;  // pair.second is a std::vector<mregion> *
  }
  models.clear();

  for (auto val : threads) {
    for (auto val2 : val->repeats) {
      delete val2;
    }
    val->splitter->DeallocSplitWindow();
    delete val->splitter;
    delete val;  // pair.second is a std::vector<mregion> *
  }
  threads.clear();
}

bool CompareRepeatOrder::operator()(RepeatRegion *lhs, RepeatRegion *rhs) {
  if (lhs->sequenceID != rhs->sequenceID) {
    return lhs->sequenceID > rhs->sequenceID;
  }

  return lhs->sequenceStart > rhs->sequenceStart;
}

