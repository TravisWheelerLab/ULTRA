//
//  ultra.cpp
//  ultraP
//

#include <climits>
#include <cmath>
#include <condition_variable>
#include <deque>
#include <map>
#include <mutex>
#include <thread>

#include "ultra.hpp"
void Ultra::AnalyzeFile() {

  if (numberOfThreads <= 0) {
    fprintf(stderr, "Number of threads must be greater than 0.\n");
    exit(-1);
  }
  reader->fastaReader->shuffle = shuffleSequence;
  if (numberOfThreads == 1) {
    reader->FillWindows();
  } else {
    // The workers start on the first windows while the rest are read; see
    // AnalyzeFileWithWorkers.
    reader->FillWindows(2 * (size_t)numberOfThreads);
  }

  OutputULTRASettings();
  InitializeWriter();

  if (numberOfThreads == 1)
    AnalyzeFileWithThread();
  else
    AnalyzeFileWithWorkers();
}

void Ultra::InitializeWriter() {
  for (int i = 0; i < writers.size(); ++i) {
    writers[i]->InitializeWriter(this, outs[i]);
  }
}

// Returns seq to the pool of free windows, refills the pool when no ready
// window is left, and takes the ready window with the smallest read ID.
SequenceWindow *Ultra::GetSequenceWindow(SequenceWindow *seq) {
  if (seq != nullptr)
    reader->AddWaitingWindow(seq);

  // Check to see if we need to read more sequence first
  if (!reader->DoneReadingFile() && !reader->IsReading() &&
      reader->ReadyWindowsSize() == 0) {
    reader->SetIsReading(true);
    reader->FillWindows();
    reader->SetIsReading(false);
  }

  return reader->GetReadyWindow();
}

// The smallest read ID of a window that may still have repeats to come.
unsigned long long Ultra::SmallestReadID() { return smallestReadID; }

// One thread: windows are analysed in read order.
void Ultra::AnalyzeFileWithThread() {
  uthread *uth = threads[0];
  SequenceWindow *currentWindow = GetSequenceWindow(nullptr);
  while (currentWindow != nullptr || !reader->DoneReadingFile()) {
    if (currentWindow != nullptr)
      AnalyzeSequenceWindow(currentWindow, uth);

    currentWindow = GetSequenceWindow(currentWindow);
  }
}

// -t N: gives the output of AnalyzeFileWithThread. Windows are read in file
// order and handed to the workers, and their results are applied in read
// order, each window returning to the pool after its results, as the
// one-thread loop does. Free windows are refilled a few at a time as they
// come back, keeping two windows per worker queued, so the workers never wait
// for a whole pool to be applied. Window contents and read ids depend only on
// file order, so this does not change the output.
//
// (ULTRA 1.2 let every thread output the repeats of a shared buffer in turn,
// so the output, and which repeats were joined across window boundaries,
// changed from run to run.)
void Ultra::AnalyzeFileWithWorkers() {
  struct Done {
    SequenceWindow *window;
    std::vector<RepeatRegion *> repeats;
  };

  std::mutex jobLock;
  std::condition_variable jobReady;
  std::deque<std::pair<size_t, SequenceWindow *>> jobs;
  bool noMoreJobs = false;

  std::mutex doneLock;
  std::condition_variable doneReady;
  std::map<size_t, Done> done;

  const size_t workers = (size_t)numberOfThreads;
  const size_t queueTarget = 2 * workers;

  std::vector<std::thread> pool;
  for (size_t w = 0; w < workers; ++w) {
    uthread *uth = threads[w];
    pool.emplace_back([&, uth]() {
      while (true) {
        std::pair<size_t, SequenceWindow *> job;
        {
          std::unique_lock<std::mutex> lock(jobLock);
          jobReady.wait(lock, [&]() { return noMoreJobs || !jobs.empty(); });
          if (jobs.empty())
            return;
          job = jobs.front();
          jobs.pop_front();
        }
        Done d;
        d.window = job.second;
        FindRepeatsInWindow(job.second, uth, d.repeats);
        {
          std::lock_guard<std::mutex> lock(doneLock);
          done[job.first] = std::move(d);
        }
        doneReady.notify_one();
      }
    });
  }

  // Read order of the next window to hand out and of the next one whose
  // results are applied.
  size_t nextSent = 0;
  size_t nextApplied = 0;
  while (true) {
    // Keep two windows per worker queued. Windows are read as soon as the
    // pool has free ones; reading order and read ids do not depend on when
    // this happens.
    while (nextSent - nextApplied < queueTarget) {
      if (reader->ReadyWindowsSize() == 0) {
        if (reader->DoneReadingFile())
          break;
        reader->SetIsReading(true);
        reader->FillWindows(workers);
        reader->SetIsReading(false);
      }
      bool sent = false;
      SequenceWindow *w;
      while ((w = reader->GetReadyWindow()) != nullptr) {
        {
          std::lock_guard<std::mutex> lock(jobLock);
          jobs.push_back(std::make_pair(nextSent, w));
        }
        jobReady.notify_one();
        ++nextSent;
        sent = true;
      }
      if (!sent) // The pool has no free window, or the file is done.
        break;
    }

    if (nextApplied == nextSent) {
      if (reader->DoneReadingFile() && reader->ReadyWindowsSize() == 0)
        break;
      continue;
    }

    // Apply results strictly in read order, as the one-thread loop does,
    // returning each window to the pool afterwards.
    Done d;
    {
      std::unique_lock<std::mutex> lock(doneLock);
      doneReady.wait(lock, [&]() { return done.count(nextApplied) > 0; });
      auto it = done.find(nextApplied);
      d = std::move(it->second);
      done.erase(it);
    }
    BeginWindow(d.window);
    FinishWindow(d.repeats);
    reader->AddWaitingWindow(d.window);
    ++nextApplied;
  }

  {
    std::lock_guard<std::mutex> lock(jobLock);
    noMoreJobs = true;
  }
  jobReady.notify_all();
  for (auto &t : pool)
    t.join();
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

// Library entry point: analyses seq as one window (it must fit the model, so
// set window_size to at least the longest string) and returns its repeats,
// which are not filtered, joined or written.
std::vector<RepeatRegion *> *
Ultra::FindRepeatsInString(const std::string &seq) {
  // A window of n symbols fills matrix columns 1..n and then column n + 1
  // without emission, so n must be below the matrix length.
  uthread *uth = this->threads[0];
  if (seq.length() >= uth->model->matrix->length) {
    fprintf(stderr,
            "ULTRA model has maximum size %llu but string has length %zu\n",
            uth->model->matrix->length - 1, seq.length());

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

  std::vector<RepeatRegion *> *repeats = new std::vector<RepeatRegion *>();
  FindRepeatsInWindow(seq_window, uth, *repeats);
  delete seq_window;

  return repeats;
}

void Ultra::AnalyzeSequenceWindow(SequenceWindow *sequence, uthread *uth) {
  BeginWindow(sequence);
  std::vector<RepeatRegion *> repeats;
  FindRepeatsInWindow(sequence, uth, repeats);
  FinishWindow(repeats);
}

// The read-ID bookkeeping for the window about to be applied.
void Ultra::BeginWindow(SequenceWindow *sequence) {
  activeReadID = sequence->readID;
  if (activeReadID < smallestReadID)
    smallestReadID = activeReadID;
}

// Hands a window's repeats to outRepeats and writes out what is complete once
// the buffer is over repeatBuffer.
void Ultra::FinishWindow(std::vector<RepeatRegion *> &repeats) {
  outRepeats.insert(outRepeats.end(), repeats.begin(), repeats.end());
  repeats.clear();
  smallestReadID = activeReadID;

  if (outRepeats.size() > repeatBuffer &&
      !this->settings->disable_streaming_out) {
    OutputRepeats();
  }
}

// Finds the repeats of one window with the model and splitter of uth. Safe to
// run on several threads at once, each with its own uthread.
void Ultra::FindRepeatsInWindow(SequenceWindow *sequence, uthread *uth,
                                std::vector<RepeatRegion *> &repeats) {

  int sleng = (int)sequence->length + (int)sequence->overlap;

  UModel *model = uth->model;
  model->matrix->RestartMatrix();

  for (int i = 0; i < sleng; ++i) {
    model->CalculateCurrentColumn(sequence, i);
  }
  model->CalculateCurrentColumnWithoutEmission();

  model->matrix->CalculateTraceback(model->matrix->previousColumnIndex);
  int i = 0;
  RepeatRegion *r = GetNextRepeat(sequence, model->matrix, &i);

  while (r != nullptr) {

    if (storeTraceAndSequence) {
      r->StoreSequence(sequence);
      r->StoreTraceback(model->matrix);
      r->StoreScores(model->matrix);

      r->GetLogoNumbers();
    }
    int splitWindow = std::max(r->repeatPeriod * settings->split_depth,
                               settings->min_split_window);
    if (r->repeatPeriod <= settings->max_split) {
      if (r->repeatLength > 2 * splitWindow) {
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

    repeats.push_back(r);
    r = GetNextRepeat(sequence, model->matrix, &i);
  }
}

// Sorts outRepeats, then joins and writes every repeat from windows more than
// 3 reads back (all of them when flush, which also ends the writers).
void Ultra::OutputRepeats(bool flush) {

  if (!flush && settings->disable_streaming_out)
    return;
  // Repeats from windows at least this far on may still be joined with the
  // next windows' repeats, so they wait.
  unsigned long long maxReadID = SmallestReadID();
  maxReadID = maxReadID > 3 ? maxReadID - 3 : 0;

  if (flush) {
    // Everything is written. (ULTRA 1.2 used 100,000,000 here, which dropped
    // every repeat from a later window.)
    maxReadID = ULLONG_MAX;
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

      // The repeat to join r with: the next one if it overlaps r with the
      // same period. Otherwise the first later one that still starts within
      // r and has r's period: a repeat of another period found inside r,
      // such as a microsatellite in a satellite's unit, sorts between r and
      // its continuation from the next window and kept the two from being
      // joined, so they were written as two overlapping records.
      long take = -1;
      if (repeats_overlap(r, outRepeats.back())) {
        take = (long)outRepeats.size() - 1;
      } else {
        unsigned long rEnd = r->sequenceStart + r->repeatLength;
        for (long k = (long)outRepeats.size() - 2; k >= 0; --k) {
          RepeatRegion *c = outRepeats[k];
          if (c->sequenceID != r->sequenceID || c->sequenceStart > rEnd)
            break;
          if (c->readID >= maxReadID) {
            // Not final yet: wait, as for the next repeat.
            outRepeats.push_back(r);
            return;
          }
          if (repeats_overlap(r, c)) {
            take = k;
            break;
          }
        }
      }
      if (take < 0)
        break;

      auto old_r = r;
      r = joint_repeat_region(r, outRepeats[take]);
      if (r->splits != nullptr)
        ValidateSplits(r->consensi, r->splits, 0.85);

      delete old_r;
      delete outRepeats[take];
      outRepeats.erase(outRepeats.begin() + take);
    }

    r->SortConsensi();
    if (r->regionScore >= scoreThreshold &&
        r->repeatLength >= settings->min_length &&
        r->repeatLength / r->repeatPeriod >= settings->min_units)
      OutputRepeat(r);

    delete r;
    r = nullptr;
  }

  if (flush) {
    for (auto writer : writers)
      writer->EndWriter();
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

  // Coverage is counted per sequence
  if (r->sequenceID != last_seq_id) {
    last_seq_id = r->sequenceID;
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
                            false);
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
    newThread->splitter = new SplitWindow();
    if (!settings->no_split)
      newThread->splitter->AllocateSplitWindow(4, settings->max_period + 1);
    else
      newThread->splitter->maxPeriod = 0;
    threads.push_back(newThread);
  }

  // Only the main thread reads windows, so the reader needs no locks.
  if (!settings->run_without_reader) {
    reader->multithread = false;
    reader->fastaReader->multithread = false;
  }
  multithreading = numberOfThreads > 1;
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
