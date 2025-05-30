//
//  ultra.cpp
//  ultraP
//

#include <cmath>

#include "ultra.hpp"
#include <cmath>
void *UltraThreadLaunch(void *dat) {
  uthread *uth = (uthread *)dat;

  uth->ultra->AnalyzeFileWithThread(dat);
  return nullptr;
}

void Ultra::AnalyzeFile() {

  if (numberOfThreads <= 0) {
    fprintf(stderr, "Number of threads must be greater than 0.\n");
    exit(-1);
  }
  reader->fastaReader->shuffle = shuffleSequence;
  reader->FillWindows();

  if (pthread_mutex_init(&outerLock, nullptr) != 0) {
    fprintf(stderr, "Failed to create outer mutex lock. Exiting.\n");
    exit(-1);
  }

  if (pthread_mutex_init(&innerLock, nullptr) != 0) {
    fprintf(stderr, "Failed to create inner mutex lock. Exiting.\n");
    exit(-1);
  }

  if (pthread_mutex_init(&repeatLock, nullptr) != 0) {
    fprintf(stderr, "Failed to create repeat lock. Exiting.\n");
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

void Ultra::InitializeWriter() {
  for (int i = 0; i < writers.size(); ++i) {
    writers[i]->InitializeWriter(this, outs[i]);
  }
}

SequenceWindow *Ultra::GetSequenceWindow(SequenceWindow *seq, uthread *uth) {

  SequenceWindow *retval = nullptr;
  bool shouldRead = false;

  if (seq != nullptr)
    reader->AddWaitingWindow(seq);

  if (multithreading)
    pthread_mutex_lock(&outerLock);

  // Check to see if we need to read more sequence first

  if (!reader->DoneReadingFile()) {
    if (!reader->IsReading()) {
      if (reader->ReadyWindowsSize() == 0) {
        shouldRead = true;
        reader->SetIsReading(true);
      }
    }
  }

  if (shouldRead) {
    reader->FillWindows();
    reader->SetIsReading(false);
  }

  retval = reader->GetReadyWindow();

  if (multithreading)
    pthread_mutex_unlock(&outerLock);

  return retval;
}

// NOTE:
// THIS IS NOT THREAD SAFE
// WE ASSUME THAT'S OK (FOR NOW)
int Ultra::SmallestReadID() {
  int smallest = 100000000;

  for (int i = 0; i < threads.size(); ++i) {
    // printf("%i: %i\n", i, threads[i]->smallestReadID);
    if (threads[i]->smallestReadID < smallest) {
      smallest = threads[i]->smallestReadID;
    }
  }
  /*
  printf("%i\n", smallest);
  printf("====================\n");*/
  return smallest;
}

void Ultra::AnalyzeFileWithThread(void *dat) {

  uthread *uth = (uthread *)dat;
  int tid = uth->id;

  SequenceWindow *currentWindow = GetSequenceWindow(nullptr, uth);
  while (currentWindow != nullptr || !reader->DoneReadingFile()) {
    if (currentWindow != nullptr)
      AnalyzeSequenceWindow(currentWindow, uth);

    currentWindow = GetSequenceWindow(currentWindow, uth);
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

double Ultra::Log2PvalForScore(float score, float period) const {
  double loc = settings->p_value_loc;
  double scale = settings->p_value_scale;
  double freq = settings->p_value_freq;
  return log2(exp(-1.0 * (score - loc) / scale) * freq);
}

double Ultra::PvalForScore(float score) const {
  double loc = settings->p_value_loc;
  double scale = settings->p_value_scale;
  double freq = settings->p_value_freq;
  return exp(-1.0 * (score - loc) / scale) * freq;
}

std::vector<RepeatRegion *> *
Ultra::FindRepeatsInString(const std::string &seq) {
  // Make sure that the seq window can fit in the DP matrix
  uthread *uth = this->threads[0];
  if (seq.length() > uth->model->matrix->length) {
    fprintf(stderr,
            "ULTRA model has maximum size %llu but string has length %zu\n",
            uth->model->matrix->length, seq.length());

    return nullptr;
  }

  // Make sure that uth isn't holding any repeats right now
  if (uth->repeats.size() > 0) {
    fprintf(stderr, "ULTRA repeat array is not empty.\n");
    return nullptr;
  }

  // Create sequence window from the string
  SequenceWindow *seq_window = new SequenceWindow(seq.length(), 0);
  seq_window->length = seq.length();
  seq_window->start = 0;
  seq_window->end = seq.length();

  // Fill sequence window
  for (int i = 0; i < seq.length(); ++i) {
    seq_window->seq[i] = SymbolForChar(seq[i]);
  }

  // Store and change primary thread
  int pthread = this->primaryThread;
  this->primaryThread = -10;

  // Run AnalyzeSequenceWindow
  this->AnalyzeSequenceWindow(seq_window, uth);

  // Gather repeats
  std::vector<RepeatRegion *> *repeats = new std::vector<RepeatRegion *>();
  for (int i = 0; i < uth->repeats.size(); ++i) {
    repeats->push_back(uth->repeats[i]);
  }

  // Clean up after ourselves
  uth->repeats.clear();
  this->primaryThread = pthread;
  delete seq_window;

  return repeats;
}

void Ultra::AnalyzeSequenceWindow(SequenceWindow *sequence, uthread *uth) {

  int sleng = (int)sequence->length + (int)sequence->overlap;

  uth->activeReadID = sequence->readID;
  if (uth->activeReadID < uth->smallestReadID)
    uth->smallestReadID = uth->activeReadID;

  UModel *model = uth->model;
  model->matrix->RestartMatrix();

  for (int i = 0; i < sleng; ++i) {
    model->CalculateCurrentColumn(sequence, i);
  }
  model->CalculateCurrentColumnWithoutEmission();

  // WE'RE RIGHT HERE
  // GOING TO TRY TO PUSH ALL THE CODE IN !!!
  // printf("Calculating traceback\n");
  model->matrix->CalculateTraceback(model->matrix->previousColumnIndex);
  int i = 0;
  RepeatRegion *r = GetNextRepeat(sequence, model->matrix, &i);

  while (r != nullptr) {

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

    uth->repeats.push_back(r);
    r = GetNextRepeat(sequence, model->matrix, &i);
  }

  if (primaryThread == uth->id) {
    outRepeats.insert(outRepeats.end(), uth->repeats.begin(),
                      uth->repeats.end());
    uth->repeats.clear();
    uth->smallestReadID = uth->activeReadID;

    if (outRepeats.size() > repeatBuffer &&
        !this->settings->disable_streaming_out) {
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
  if (maxReadID < 0)
    maxReadID = 0;

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
        auto old_r = r;
        r = joint_repeat_region(r, outRepeats.back());
        if (r->splits != nullptr)
          ValidateSplits(r->consensi, r->splits, 0.85);

        delete old_r;
        delete outRepeats.back();
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
    for (auto writer : writers)
      writer->EndWriter();
    // fprintf(out, "]\n}\n");
  }

  for (auto out : outs)
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

  if (r->readID > last_read_id) {
    last_rep_end = 0;
  }

  auto rep_start = r->sequenceStart;
  if (rep_start < last_rep_end) {
    rep_start = last_rep_end;
  }

  if (rep_start < r->sequenceStart + r->repeatLength) {
    total_coverage += (r->sequenceStart + r->repeatLength) - rep_start;
  }

  if (r->sequenceStart + r->repeatLength > last_rep_end) {
    last_rep_end = r->sequenceStart + r->repeatLength;
  }

  if (!settings->suppress_out)
    for (auto writer : writers)
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

unsigned long long Ultra::Coverage() { return total_coverage; }

Ultra::Ultra(Settings *s) {
  settings = s;

  settings_out = stdout;

  if (!settings->out_file.empty()) {

    int c = 0;
    if (settings->ultra_out)
      c++;
    if (settings->json_out)
      c++;
    if (settings->bed_out)
      c++;
    if (c > 1) {
      if (settings->ultra_out) {
        std::string ultra_path = settings->out_file + ".tsv";
        FILE *out = fopen(ultra_path.c_str(), "w");
        if (out == NULL) {
          fprintf(stderr, "Unable to open output file %s\n",
                  ultra_path.c_str());
          exit(-1);
        }
        outs.push_back(out);
        writers.push_back(new TabFileWriter());
      }

      if (settings->json_out) {
        std::string json_path = settings->out_file + ".json";
        FILE *out = fopen(json_path.c_str(), "w");
        if (out == NULL) {
          fprintf(stderr, "Unable to open output file %s\n", json_path.c_str());
          exit(-1);
        }
        outs.push_back(out);
        writers.push_back(new JSONFileWriter());
      }

      if (settings->bed_out) {
        std::string bed_path = settings->out_file + ".bed";
        FILE *out = fopen(bed_path.c_str(), "w");
        if (out == NULL) {
          fprintf(stderr, "Unable to open output file %s\n", bed_path.c_str());
          exit(-1);
        }
        outs.push_back(out);
        writers.push_back(new BEDFileWriter());
      }
    }

    else {
      FILE *out = fopen(settings->out_file.c_str(), "w");
      if (out == NULL) {
        fprintf(stderr, "Unable to open output file %s\n",
                settings->out_file.c_str());
        exit(-1);
      }
      outs.push_back(out);
      if (settings->ultra_out)
        writers.push_back(new TabFileWriter());
      else if (settings->json_out)
        writers.push_back(new JSONFileWriter());
      else if (settings->bed_out)
        writers.push_back(new BEDFileWriter());
    }

    std::string settings_file_path = settings->out_file + ".settings";
    if (!settings->hide_settings) {
      settings_out = fopen(settings_file_path.c_str(), "w");
      if (settings_out == NULL) {
        fprintf(stderr, "Unable to open settings output file %s\n",
                settings_file_path.c_str());
        exit(-1);
      }
    }

  } else {
    if (settings->ultra_out)
      writers.push_back(new TabFileWriter());
    else if (settings->json_out)
      writers.push_back(new JSONFileWriter());
    else if (settings->bed_out)
      writers.push_back(new BEDFileWriter());
    outs.push_back(stdout);
  }

  numberOfThreads = settings->threads;
  scoreThreshold = settings->min_score;
  outputReadID = settings->show_wid;
  outputRepeatSequence = settings->show_seq;

  passID = 0;
  if (!settings->run_without_reader) {
    reader = new FileReader(settings->in_file, settings->windows,
                            settings->window_size, settings->overlap,
                            settings->threads > 1);
    reader->fastaReader->shuffle = shuffleSequence;
  }

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
                             settings->max_delete, leng + 1);

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
    if (!settings->run_without_reader) {
      reader->multithread = false;
    }
    multithreading = false;
  }

  else {
    if (!settings->run_without_reader) {
      reader->multithread = true;
    }
    multithreading = true;
  }
}

Ultra::~Ultra() {
  settings = nullptr;

  delete reader;
  reader = nullptr;

  for (auto writer : writers) {
    delete writer;
  }
  writers.clear();

  for (auto out : outs) {
    if (out != stdout && out != nullptr) {
      fclose(out);
    }
  }
  outs.clear();

  for (auto &pair : masks_for_seq) {
    delete pair.second; // pair.second is a std::vector<mregion> *
  }
  masks_for_seq.clear();

  for (auto val : outRepeats) {
    delete val; // pair.second is a std::vector<mregion> *
  }
  outRepeats.clear();

  for (auto val : models) {
    delete val; // pair.second is a std::vector<mregion> *
  }
  models.clear();

  for (auto val : threads) {
    for (auto val2 : val->repeats) {
      delete val2;
    }
    val->splitter->DeallocSplitWindow();
    delete val->splitter;
    delete val; // pair.second is a std::vector<mregion> *
  }
  threads.clear();
}

bool CompareRepeatOrder::operator()(RepeatRegion *lhs, RepeatRegion *rhs) {
  if (lhs->sequenceID != rhs->sequenceID) {
    return lhs->sequenceID > rhs->sequenceID;
  }

  return lhs->sequenceStart > rhs->sequenceStart;
}
