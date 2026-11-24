//
//  repeat.cpp
//  ultraP
//

#include "repeat.hpp"
#include "RepeatSplitter.hpp"
#include <algorithm>
#include <vector>
// The indel test compares two periods of bases, and at least this many.
#define INDEL_TEST_MIN_BASES 12
// Added to the agreement score of the indel test for a base that equals the
// first copy's base in its column. A matching base counts the same in the
// logo.
#define INDEL_TEST_FIRST_COPY_WEIGHT 4

// Scores for aligning a repeat to its consensus, and how many times the
// alignment is repeated with the consensus it produced. In a simulation of
// 5,200 noisy tandem arrays, these settings gave a consensus that was a
// rotation of the planted unit in 2,798. Eight other score settings gave
// between 2,791 and 2,798, and a third pass gave no more.
#define REALIGN_MATCH 2
#define REALIGN_MISMATCH -2
#define REALIGN_GAP -5
#define REALIGN_PASSES 2

// Builds the logo (symbol counts per unit position) and accumulates
// mismatches, insertions and deletions.
//
// The three counters come from one walk along the matrix traceback that takes
// every insertion and deletion the model calls as real, so they agree with the
// stored traceback. That walk puts later copies in the wrong columns when the
// model calls a deletion in a copy that is full length, or misses an indel, so
// the logo is built in two steps:
//
// 1. The walk is repeated, and at each insertion or deletion the bases that
//    follow are compared with the logo so far, once at the column the indel
//    implies and once at the column without it. The indel is ignored when the
//    column without it agrees better.
// 2. The whole repeat is aligned to the consensus of that logo, treating the
//    consensus as a circle, and each base is counted in the column it aligns
//    to. This is done REALIGN_PASSES times.
//
// In a simulation of 5,200 noisy tandem arrays, this gave a consensus that was
// a rotation of the planted unit in 2,798, against 2,534 for the single walk.
void RepeatRegion::CreateLogo(SequenceWindow *window, UMatrix *matrix) {
  std::vector<int> walkLogo;
  int walkMismatches = 0, walkInsertions = 0, walkDeletions = 0;
  WalkLogo(window, matrix, false, walkLogo, walkMismatches, walkInsertions,
           walkDeletions, &events);
  mismatches += walkMismatches;
  insertions += walkInsertions;
  deletions += walkDeletions;

  // With no indel on the path there is nothing to test, and the second walk
  // would repeat the first.
  if (walkInsertions != 0 || walkDeletions != 0) {
    int m = 0, ins = 0, del = 0;
    WalkLogo(window, matrix, true, walkLogo, m, ins, del, nullptr);
  }
  for (int pass = 0; pass < REALIGN_PASSES; ++pass) {
    walkLogo = RealignLogo(walkLogo, window);
  }

  logoMemory = (int *)malloc(sizeof(int) * repeatPeriod * (NUM_SYMBOLS + 1));
  logo = (int **)malloc(sizeof(int *) * repeatPeriod);
  for (int i = 0; i < repeatPeriod; ++i) {
    logo[i] = &logoMemory[i * NUM_SYMBOLS];
    for (int j = 0; j < NUM_SYMBOLS; ++j) {
      logo[i][j] = walkLogo[i * NUM_SYMBOLS + j];
    }
  }
}

// One walk along the matrix traceback from windowStart, comparing each base
// with the base the model aligned it to in the previous copy. A matching base
// counts 4 in its column of the logo and any other base counts 1. testIndels
// turns on the indel test described at CreateLogo. When outEvents is not
// null it receives the substitutions, insertions and deletions at each
// position of the repeat (3 counts per position).
void RepeatRegion::WalkLogo(SequenceWindow *window, UMatrix *matrix,
                            bool testIndels, std::vector<int> &outLogo,
                            int &outMismatches, int &outInsertions,
                            int &outDeletions,
                            std::vector<unsigned int> *outEvents) {
  const long period = repeatPeriod;
  const unsigned long periodU = (unsigned long)(long)repeatPeriod;
  const symbol *seq = window->seq;
  const unsigned long seqCap = window->maxLength + window->maxOverlap + 1;
  const unsigned long repeatEnd = windowStart + repeatLength;

  std::vector<int> lg((size_t)repeatPeriod * NUM_SYMBOLS, 0);
  if (outEvents != nullptr)
    outEvents->assign(3 * (size_t)repeatLength, 0);

  // Adds n to count kind (0 substitution, 1 insertion, 2 deletion) at window
  // position at.
  auto note = [&](unsigned long at, int kind, int n) {
    if (outEvents == nullptr || at < windowStart)
      return;
    unsigned long k = at - windowStart;
    if (k < repeatLength)
      (*outEvents)[3 * k + kind] += (unsigned int)n;
  };

  unsigned long endpoint = repeatLength - periodU;
  unsigned long fplace = periodU + windowStart;
  unsigned long bplace = windowStart;

  unsigned long winLen = window->length + window->overlap;
  if (endpoint > winLen)
    endpoint = winLen;

  int p = 0;

  // Traceback columns of the indels the test ignored, in increasing order,
  // and how far the previous-copy side has read through them.
  std::vector<unsigned long> ignored;
  size_t ignoredSeen = 0;

  // How well the bases from start on agree with the logo so far and with the
  // first copy, if the base at start sits in column col0.
  auto agreement = [&](unsigned long start, long col0) -> int {
    unsigned long len = std::max(2 * periodU, (unsigned long)INDEL_TEST_MIN_BASES);
    int score = 0;
    for (unsigned long j = 0; j < len; ++j) {
      unsigned long pos = start + j;
      if (pos >= repeatEnd || pos >= seqCap)
        break;
      long col = (col0 + (long)j) % period;
      if (col < 0)
        col += period;
      int sym = seq[pos];
      score += lg[col * NUM_SYMBOLS + sym];
      if (seq[windowStart + col] == sym)
        score += INDEL_TEST_FIRST_COPY_WEIGHT;
    }
    return score;
  };

  int mism = 0, ins = 0, del = 0;

  for (unsigned long i = 0; i < endpoint; ++i) {
    int frow = matrix->traceback[fplace];
    int brow = matrix->traceback[bplace];

    cell fdesc = matrix->cellDescriptions[frow];
    cell bdesc = matrix->cellDescriptions[brow];

    if (i < periodU)
      bdesc = matrix->cellDescriptions[0];

    if (testIndels) {
      // A deletion of k moves the next base k columns on; an insertion of n
      // leaves the base after it where it was.
      bool keep = true;
      if (fdesc.type == CT_DELETION) {
        int k = fdesc.indelNumber;
        keep = agreement(fplace, p + k) >= agreement(fplace, p);
      } else if (fdesc.type == CT_INSERTION) {
        int n = fdesc.indelNumber;
        unsigned long after = fplace + (unsigned long)(long)n;
        keep = agreement(after, p) >= agreement(after, p + n);
      }
      if (!keep) {
        ignored.push_back(fplace);
        fdesc = matrix->cellDescriptions[0];
      }

      // The same indel, met again one copy later on the other side.
      while (ignoredSeen < ignored.size() && ignored[ignoredSeen] < bplace)
        ++ignoredSeen;
      if (ignoredSeen < ignored.size() && ignored[ignoredSeen] == bplace &&
          (bdesc.type == CT_DELETION || bdesc.type == CT_INSERTION))
        bdesc = matrix->cellDescriptions[0];
    }

    if (fdesc.type == CT_INSERTION) {
      ins += fdesc.indelNumber;
      note(fplace, 1, fdesc.indelNumber);
      i += fdesc.indelNumber;
      fplace += fdesc.indelNumber;
    } else if (fdesc.type == CT_DELETION) {
      del += fdesc.indelNumber;
      note(fplace, 2, fdesc.indelNumber);
      bplace += fdesc.indelNumber;
      p += fdesc.indelNumber;
    }

    if (bdesc.type == CT_INSERTION) {
      bplace += bdesc.indelNumber;
    } else if (bdesc.type == CT_DELETION) {
      // The bases skipped here sit at the unit positions the previous copy
      // deleted. They are counted in those columns (weight 1, as there is no
      // base to compare them with) and the phase moves past them.
      for (int j = 0; j < bdesc.indelNumber; ++j) {
        unsigned long fp = fplace + j;
        if (fp < repeatEnd)
          lg[((p + j) % repeatPeriod) * NUM_SYMBOLS + seq[fp]] += 1;
      }
      fplace += bdesc.indelNumber;
      i += bdesc.indelNumber;
      p += bdesc.indelNumber;
    }

    if (fplace >= repeatEnd)
      break;

    symbol s = seq[fplace];
    symbol sb = seq[bplace];

    p = p % repeatPeriod;

    if (s != sb) {
      mism += 1;
      note(fplace, 0, 1);
      lg[p * NUM_SYMBOLS + s] += 1;
    } else {
      lg[p * NUM_SYMBOLS + s] += 4;
    }

    ++p;
    ++fplace;
    ++bplace;
  }

  outLogo.swap(lg);
  outMismatches = mism;
  outInsertions = ins;
  outDeletions = del;
}

// Aligns the repeat to the consensus of logo and returns a logo with one
// count for each aligned base, in the column it aligns to.
//
// The consensus is treated as a circle, with column 0 following the last
// column, so the alignment runs through as many copies as the repeat has and
// may start in any column. A base is aligned to a column (REALIGN_MATCH or
// REALIGN_MISMATCH), or is an insertion, or a column is skipped (both
// REALIGN_GAP). An unknown base (N_UNK) scores 0 in every column. Bases
// called insertions and unknown bases are not counted, and a column that
// ends with no counts keeps its counts from logo. The traceback takes one
// byte per base and column.
std::vector<int> RepeatRegion::RealignLogo(const std::vector<int> &logo,
                                           SequenceWindow *window) {
  const unsigned char DIAGONAL = 0;
  const unsigned char INSERTED_BASE = 1;
  const unsigned char SKIPPED_COLUMN = 2;

  std::vector<int> newLogo(logo.size(), 0);
  size_t period = repeatPeriod > 0 ? (size_t)repeatPeriod : 0;
  size_t start = windowStart;
  size_t seqCap = window->maxLength + window->maxOverlap + 1;
  size_t n = start < seqCap ? std::min((size_t)repeatLength, seqCap - start) : 0;
  if (period == 0 || n == 0)
    return newLogo;
  const symbol *seq = window->seq;

  // Arg-max symbol per column, lowest symbol on ties, as in
  // CreateConsensusFromLogo.
  std::vector<int> consensus(period, 0);
  for (size_t j = 0; j < period; ++j) {
    const int *row = &logo[j * NUM_SYMBOLS];
    int best = 0;
    for (int t = 1; t < NUM_SYMBOLS; ++t) {
      if (row[t] > row[best])
        best = t;
    }
    consensus[j] = best;
  }

  std::vector<unsigned char> traceback(n * period, DIAGONAL);
  std::vector<int> prev(period, 0);
  std::vector<int> cur(period, 0);
  for (size_t i = 0; i < n; ++i) {
    int sym = seq[start + i];
    unsigned char *moves = &traceback[i * period];
    for (size_t j = 0; j < period; ++j) {
      size_t before = (j + period - 1) % period;
      int diagonal = prev[before] + (sym == N_UNK               ? 0
                                     : sym == consensus[j]      ? REALIGN_MATCH
                                                                : REALIGN_MISMATCH);
      int inserted = prev[j] + REALIGN_GAP;
      if (diagonal >= inserted) {
        cur[j] = diagonal;
        moves[j] = DIAGONAL;
      } else {
        cur[j] = inserted;
        moves[j] = INSERTED_BASE;
      }
    }
    // A run of skipped columns may pass from the last column to column 0, so
    // the row is swept twice.
    for (int sweep = 0; sweep < 2; ++sweep) {
      for (size_t j = 0; j < period; ++j) {
        size_t before = (j + period - 1) % period;
        int skipped = cur[before] + REALIGN_GAP;
        if (skipped > cur[j]) {
          cur[j] = skipped;
          moves[j] = SKIPPED_COLUMN;
        }
      }
    }
    prev.swap(cur);
  }

  // Best column for the last base (lowest column on ties), then back.
  size_t j = 0;
  for (size_t t = 1; t < period; ++t) {
    if (prev[t] > prev[j])
      j = t;
  }
  size_t i = n;
  size_t skippedInARow = 0;
  while (i > 0) {
    unsigned char move = traceback[(i - 1) * period + j];
    if (move == DIAGONAL) {
      int sym = seq[start + i - 1];
      if (sym != N_UNK)
        newLogo[j * NUM_SYMBOLS + sym] += 1;
      j = (j + period - 1) % period;
      --i;
      skippedInARow = 0;
    } else if (move == INSERTED_BASE) {
      --i;
      skippedInARow = 0;
    } else {
      j = (j + period - 1) % period;
      ++skippedInARow;
      // Skipping every column gains nothing, so the alignment never does it;
      // this only guards the loop.
      if (skippedInARow > period)
        break;
    }
  }

  for (size_t c = 0; c < period; ++c) {
    bool empty = true;
    for (int t = 0; t < NUM_SYMBOLS; ++t) {
      if (newLogo[c * NUM_SYMBOLS + t] != 0) {
        empty = false;
        break;
      }
    }
    if (empty) {
      for (int t = 0; t < NUM_SYMBOLS; ++t)
        newLogo[c * NUM_SYMBOLS + t] = logo[c * NUM_SYMBOLS + t];
    }
  }
  return newLogo;
}

void RepeatRegion::CreateLogoWithoutMatrix() {

  if (sequence.length() == 0 || traceback.length() == 0) {
    return;
  }

  logoMemory = (int *)malloc(sizeof(int) * repeatPeriod * (NUM_SYMBOLS + 1));
  logo = (int **)malloc(sizeof(int *) * repeatPeriod);

  for (int i = 0; i < repeatPeriod; ++i) {
    int p = i * NUM_SYMBOLS;
    logo[i] = &logoMemory[p];
    for (int j = 0; j < NUM_SYMBOLS; ++j) {
      logo[i][j] = 0;
    }
  }

  int p = 0;
  for (int i = 0; i < repeatLength; ++i) {
    if (traceback[i] == '.') {
      logo[p % repeatPeriod][SymbolForChar(sequence[i])]++;
      ++p;
    }

    else if (traceback[i] == 'I') {
      // Do nothing
    }

    else if (traceback[i] == 'D') {
      // Count how many deletions
      int dels = 1;
      for (int j = 1; traceback[i + j] == 'D'; ++j) {
        ++dels;
      }

      for (int j = 0; j < dels; ++j) {
        logo[(p + dels) % repeatPeriod][SymbolForChar(sequence[i])]++;
        ++p;
      }
    }
  }
}

void RepeatRegion::CreateConsensusFromLogo() {
  consensus = (symbol *)malloc(sizeof(symbol) * repeatPeriod);
  for (int i = 0; i < repeatPeriod; ++i) {
    symbol s = 0;

    for (int j = 1; j < NUM_SYMBOLS; ++j) {
      if (logo[i][j] > logo[i][s])
        s = j;
    }
    consensus[i] = s;
  }

  this->string_consensus = GetConsensus();
}

void RepeatRegion::StoreSequence(SequenceWindow *window) {
  sequence = "";
  sequence.reserve(repeatLength + 1);

  for (int i = 0; i < repeatLength; ++i) {
    sequence.push_back(CharForSymbol(window->seq[i + windowStart]));
  }
  sequence[repeatLength] = '\0';
}

RepeatRegion::RepeatRegion() {
  sequenceName = "";
  sequence = "";
  sequenceID = 0;
  windowStart = 0;
  sequenceStart = 0;
  repeatLength = 0;
  repeatPeriod = 0;
  startScore = 0;
  endScore = 0;
  regionScore = 0;
  readID = 0;
  falseStart = 0;

  mismatches = 0;
  insertions = 0;
  deletions = 0;

  winOverlapSize = 0;
  overlapCorrection = OC_NONE;

  scores = NULL;
  logoNumbers = NULL;
  logoMemory = NULL;
  logo = NULL;
  consensus = NULL;

  consensi = nullptr;
  splits = nullptr;

  lookForward = NULL;
  lookBack = NULL;

  forwardCounts = NULL;
  backCounts = NULL;

  combinedRepeat = false;
}

RepeatRegion::~RepeatRegion() {
  if (logo != nullptr) {
    free(logo);
    logo = nullptr;
  }

  if (consensus != nullptr) {
    free(consensus);
    consensus = nullptr;
  }

  if (logoMemory != nullptr) {
    free(logoMemory);
    logoMemory = nullptr;
  }

  if (logoNumbers != nullptr) {
    free(logoNumbers);
    logoNumbers = nullptr;
  }

  if (lookBack != nullptr) {
    free(lookBack);
    lookBack = nullptr;
  }

  if (lookForward != nullptr) {
    free(lookForward);
    lookForward = nullptr;
  }

  if (forwardCounts != nullptr) {
    free(forwardCounts);
    forwardCounts = nullptr;
  }

  if (backCounts != nullptr) {
    free(backCounts);
    backCounts = nullptr;
  }

  if (scores != nullptr) {
    free(scores);
    scores = nullptr;
  }

  if (splits != nullptr) {
    splits->erase(splits->begin(), splits->end());
    delete splits;
    splits = nullptr;
  }

  if (consensi != nullptr) {
    consensi->erase(consensi->begin(), consensi->end());
    delete consensi;
    consensi = nullptr;
  }
}

std::string RepeatRegion::GetConsensus() {
  std::string con = "";
  con.reserve(repeatPeriod);
  for (int i = 0; i < repeatPeriod; ++i) {
    con.push_back(CharForSymbol(consensus[i]));
  }
  // con[repeatPeriod] = '\0';
  return con;
}

void RepeatRegion::StoreTraceback(UMatrix *matrix) {

  for (int i = 0; i < repeatPeriod; ++i) {
    traceback.push_back('.');
  }

  unsigned long end = windowStart + repeatLength;
  for (unsigned long i = windowStart + repeatPeriod; i < end; ++i) {
    cell c = matrix->cellDescriptions[matrix->traceback[i]];

    if (c.type == CT_INSERTION) {
      i += c.indelNumber - 1;
      for (int j = 0; j < c.indelNumber; ++j) {
        traceback.push_back('I');
      }
    }

    else if (c.type == CT_DELETION) {
      i += c.indelNumber - 1;
      // traceback[traceback.length() - 1] = 'D';
      traceback.push_back('D');
      for (int j = 1; j < c.indelNumber; ++j) {
        traceback.push_back('D');
      }
      // traceback.push_back('.');
    }

    else {
      traceback.push_back('.');
    }
  }
}

void RepeatRegion::StoreScores(UMatrix *matrix) {

  scores = (float *)malloc(sizeof(float) * repeatLength);

  for (unsigned long i = 1; i <= repeatLength; ++i) {

    unsigned long p = i + windowStart;
    int cell = matrix->traceback[p];
    if (matrix->cellDescriptions[cell].type != CT_MATCH)
      cell = matrix->cellDescriptions[cell].parentIndex;
    float s = matrix->ScoreAt(p, cell);

    int prevCell = matrix->traceback[p - 1];
    if (matrix->cellDescriptions[prevCell].type != CT_MATCH)
      prevCell = matrix->cellDescriptions[prevCell].parentIndex;
    float ps = matrix->ScoreAt(p - 1, prevCell);

    scores[i - 1] = s - ps;
  }
}

void RepeatRegion::LookBackDistance() {
  // We can't look back before the index repeatPeriod
  // so set those values to 0
  for (int i = 0; i < repeatPeriod; ++i) {
    lookBack[i] = 0;
  }

  int indelNum = 0;

  int inPeriod = 0;
  int delPeriod = 0;

  int indelStop = 0;
  int delSwitch = 0;

  for (int i = repeatPeriod; i < repeatLength; ++i) {
    // Assume we aren't in an indel state and correct if we are
    lookBack[i] = repeatPeriod;

    // We are not in an indel state...
    // Check to see if we are starting one
    if (indelNum == 0) {

      if (traceback[i] == 'D') {
        // Count all the Ds
        indelStop = i + repeatPeriod;

        delPeriod = repeatPeriod - 1;
        inPeriod = repeatPeriod + repeatPeriod - 1;

        indelNum = 1;

        for (int j = 1; traceback[i + j] == 'D'; ++j) {
          --delPeriod;
          --inPeriod;
          ++indelNum;
        }

        lookBack[i] = delPeriod;

      }

      else if (traceback[i] == 'I') {
        indelStop = i + repeatPeriod;

        indelNum = 1;
        inPeriod = repeatPeriod + 1;

        lookBack[i] = 0;
      }
    }

    // We are in an indel state!
    else {
      // Check to see if there are more insertions
      if (traceback[i] == 'I') {
        ++indelStop;

        ++inPeriod;
        ++indelNum;

        lookBack[i] = 0;
      }

      // We are in a deletion state
      else if (delPeriod > 0) {
        lookBack[i] = delPeriod;

        if (i >= indelStop) {
          lookBack[i] = repeatPeriod;
          indelNum = 0;
          inPeriod = 0;
          delPeriod = 0;
        }

        else if (i >= indelStop - indelNum) {
          lookBack[i] = inPeriod;
        }

      }

      else if (traceback[i] == '.') {
        // Check to see if we're done with the indel
        if (i > indelStop) {
          indelNum = 0;
          inPeriod = 0;
          delPeriod = 0;
        }

        // We are in an insertion
        else if (inPeriod > 0) {
          lookBack[i] = inPeriod;
        }
      }
    }
  }
}

void RepeatRegion::LookForwardDistance() {
  // We can't look forward after the index repeatLength - repeatPeriod
  // so set those values to 0
  for (int i = repeatLength - repeatPeriod; i < repeatLength; ++i) {
    lookForward[i] = 0;
  }

  int indelNum = 0;

  int inPeriod = 0;
  int delPeriod = 0;

  int indelStop = 0;
  int delSwitch = 0;

  for (int i = repeatLength - repeatPeriod - 1; i >= 0; --i) {

    // Assume that there hasn't been an indel
    lookForward[i] = repeatPeriod;

    // We are not in an indel state
    if (indelNum == 0) {
      if (traceback[i] == 'D') {
        // Count all the Ds

        delPeriod = repeatPeriod;
        inPeriod = repeatPeriod + repeatPeriod;

        indelNum = 0;
        for (; traceback[i] == 'D'; --i) {
          --delPeriod;
          --inPeriod;
          ++indelNum;
          lookForward[i] = repeatPeriod;
        }

        indelStop = i - repeatPeriod + 1;
        lookForward[i] = delPeriod;
      }

      else if (traceback[i] == 'I') {
        indelStop = i - repeatPeriod;
        indelNum = 1;
        inPeriod = repeatPeriod + 1;
        lookForward[i] = 0;
      }
    }

    // We are in a deletion state
    else if (delPeriod > 0) {
      lookForward[i] = delPeriod;

      if (i < indelStop) {
        lookForward[i] = repeatPeriod;
        indelNum = 0;
        inPeriod = 0;
        delPeriod = 0;
      }

      else if (i < indelStop + indelNum) {
        lookForward[i] = inPeriod;
      }

    }

    // We are in an insertion state
    else {

      if (traceback[i] == '.') {
        // We're done
        if (i < indelStop) {
          indelNum = 0;
          inPeriod = 0;
          delPeriod = 0;
        }

        // We are in a J state
        else if (inPeriod > 0) {
          lookForward[i] = inPeriod;
        }
      }

      else if (traceback[i] == 'I') {
        --indelStop;
        ++inPeriod;
        ++indelNum;
        lookForward[i] = 0;
      }
    }
  }
}

void RepeatRegion::CountFromDistance(int *distance, char *counts, const int dir,
                                     const int depth) {

  for (int i = 0; i < repeatLength; ++i) {

    char count = 0;

    // We don't perform computations on insertions
    if (distance[i] == 0) {
      counts[i] = -1;
      continue;
    }

    int p = i;
    for (int j = 0; j < depth; ++j) {

      p += (distance[p] * dir);

      if (p < 0 || p >= repeatLength)
        break;

      if (sequence[i] == sequence[p])
        ++count;

      if (distance[p] == 0)
        break;
    }

    counts[i] = count;
  }
}

void RepeatRegion::CountForwardBackwardConsensus(const int depth) {
  lookBack = (int *)malloc(sizeof(int) * repeatLength);
  lookForward = (int *)malloc(sizeof(int) * repeatLength);

  LookBackDistance();
  LookForwardDistance();
}

std::vector<RepeatRegion *> *RepeatRegion::SplitRepeats(const int depth,
                                                        const int cutoff) {
  CountForwardBackwardConsensus(depth);

  const int invCutoff = depth - cutoff;

  int repStart = 0;

  std::vector<RepeatRegion *> *subReps = NULL;

  for (int i = 0; i < repeatLength; ++i) {
    // Don't split on insertions or at
    // the beginning/end of the sequence
    if (backCounts[i] < 0 || forwardCounts[i] < 0 || lookForward[i] == 0)
      continue;

    // Check to see if the ith position has the right conditions
    if (backCounts[i] > invCutoff && forwardCounts[i] < cutoff) {
      int f = i + lookForward[i];

      // Don't split on insertions or at the beginning/end of the sequence
      if (backCounts[f] < 0 || forwardCounts[f] < 0 || lookBack[i] == 0)
        continue;

      // Split it!
      if (backCounts[f] < cutoff && forwardCounts[f] > invCutoff) {
        if (subReps == NULL)
          subReps = new std::vector<RepeatRegion *>();

        subReps->push_back(SubRepeat(repStart, f - repStart));
        repStart = f;
        i = f; // skip ahead to next unit
      }
    }
  }

  if (repStart > 0) {
    subReps->push_back(SubRepeat(repStart, (int)(repeatLength - repStart)));
  }

  return subReps;
}

RepeatRegion *RepeatRegion::SubRepeat(int start, int length) {
  RepeatRegion *region = new RepeatRegion();

  region->sequenceName = sequenceName;

  if (traceback.length() > 0)
    region->traceback = traceback.substr(start, length);
  if (sequence.length() > 0)
    region->sequence = sequence.substr(start, length);

  region->sequenceID = sequenceID;
  region->readID = readID;

  region->winOverlapSize = winOverlapSize;
  region->winTotalLength = winTotalLength;

  region->combinedRepeat = combinedRepeat;

  // Don't know if I should add start to falseStart...
  region->falseStart = falseStart;

  region->windowStart = windowStart + start;
  region->sequenceStart = sequenceStart + start;

  region->repeatLength = length;

  region->repeatPeriod = repeatPeriod;

  // Add real score capabilities at some point...
  float scorePerSymbol = regionScore / ((float)repeatLength);

  region->startScore = startScore + (scorePerSymbol * ((float)start));
  region->endScore = region->startScore + (scorePerSymbol * ((float)length));
  region->regionScore = region->endScore - region->startScore;

  if (region->sequence.length() > 0 && region->traceback.length() > 0) {
    region->CreateLogoWithoutMatrix();
    region->CreateConsensusFromLogo();
  }

  return region;
}

// This assumes that matrix has already calculated the traceback
RepeatRegion *GetNextRepeat(SequenceWindow *window, UMatrix *matrix, int *pos) {
  if (*pos == 0)
    *pos = 1;

  bool foundRepeat = false;
  int i;

  unsigned long seqLength = window->length + window->overlap;

  RepeatRegion *region = nullptr;

  for (i = *pos; i < seqLength; ++i) {

    if (matrix->traceback[i] > 0) {
      foundRepeat = true;
      break;
    }
  }

  if (foundRepeat) {
    int row = matrix->traceback[i];
    cell desc = matrix->cellDescriptions[row];

    region = new RepeatRegion();
    region->winOverlapSize = window->overlap;
    region->winTotalLength = window->length + window->overlap;
    region->sequenceName = window->sequenceName;
    region->sequenceID = window->seqID;
    region->readID = window->readID;
    region->repeatPeriod = desc.order;

    region->windowStart = i - region->repeatPeriod - 1;
    region->sequenceStart =
        window->start + (i - window->overlap - region->repeatPeriod) - 1;

    region->startScore = matrix->ScoreAt(i, row);

    int length = region->repeatPeriod;

    // We use a <= here because we calculate one additional column
    for (; i <= seqLength; ++i) {
      if (matrix->traceback[i] != 0) {
        ++length;
      }

      else {
        ++i; // always make i = end + 1
        break;
      }
    }

    region->repeatLength = length;

    row = matrix->traceback[i - 1];

    region->endScore = matrix->ScoreAt(i - 1, row);
    region->regionScore = region->endScore - region->startScore;
  }

  *pos = i;

  if (region != nullptr) {
    region->readID = window->readID;
    region->CreateLogo(window, matrix);
    region->CreateConsensusFromLogo();
  }

  return region;
}

void RepeatRegion::GetLogoNumbers() {

  unsigned long start = this->windowStart;
  unsigned long length = this->repeatLength;
  logoNumbers = (int *)malloc(sizeof(int) * (length + 1));

  int num = 0;
  for (unsigned long i = 0; i < length; ++i) {

    if (traceback[i] == '.') {
      logoNumbers[i] = num % repeatPeriod;
      ++num;
    }

    if (traceback[i] == 'I') {
      logoNumbers[i] = -1;
    }

    int delCount = 0;
    for (int j = 0; traceback[i + j] == 'D'; ++j) {
      ++num;
      ++delCount;
    }

    for (int j = 0; j < delCount; ++j) {
      logoNumbers[i] = num % repeatPeriod;

      ++num;
      ++i;
    }
    if (delCount > 0)
      --i;
  }
}

bool repeats_overlap(RepeatRegion *r1, RepeatRegion *r2,
                     bool require_same_period) {
  if (r1->sequenceID == r2->sequenceID) {
    if (require_same_period) {
      if (r1->repeatPeriod == r2->repeatPeriod) {
        return (r1->sequenceStart + r1->repeatLength >= r2->sequenceStart);
      } else {
        return false;
      }
    }
    return (r1->sequenceStart + r1->repeatLength >= r2->sequenceStart);
  }
  return false;
}

// This does not check if r1 is dominated by r2
// This does not calculate pval for region
RepeatRegion *joint_repeat_region(RepeatRegion *r1, RepeatRegion *r2) {

  if (!repeats_overlap(r1, r2, false))
    return nullptr;

  RepeatRegion *joint_rep = new RepeatRegion();

  joint_rep->string_consensus = r1->string_consensus;
  joint_rep->sequenceStart = r1->sequenceStart;

  // The joint region is the first s1_seq_len positions of r1 followed by r2
  // without its first r2_skip positions. If r2 ends inside r1 the joint
  // region is all of r1 and none of r2; taking r1 up to the start of r2 and
  // then all of r2 would cut r1 short.
  bool r2_inside = r2->sequenceStart + r2->repeatLength <=
                   r1->sequenceStart + r1->repeatLength;
  int s1_seq_len = r2->sequenceStart - r1->sequenceStart;
  int r2_skip = 0;
  if (r2_inside) {
    s1_seq_len = r1->repeatLength;
    r2_skip = r2->repeatLength;
  }
  int overlap = (r1->sequenceStart + r1->repeatLength) - r2->sequenceStart;

  joint_rep->repeatLength = s1_seq_len + (r2->repeatLength - r2_skip);

  joint_rep->sequenceName = r1->sequenceName;
  joint_rep->sequenceID = r1->sequenceID;
  joint_rep->readID = r1->readID;
  joint_rep->winOverlapSize = r1->winOverlapSize;
  joint_rep->winTotalLength = r1->winTotalLength;

  joint_rep->combinedRepeat = true;

  joint_rep->repeatPeriod = r1->repeatPeriod;

  joint_rep->overlapCorrection = r1->overlapCorrection + 1;

  // We have to join ->sequence
  if (!r1->sequence.empty() && !r2->sequence.empty()) {
    joint_rep->sequence =
        r1->sequence.substr(0, s1_seq_len) +
        r2->sequence.substr(std::min((size_t)r2_skip, r2->sequence.size()));
  }

  // Check if we have to join ->traceback
  if (!r1->traceback.empty() && !r2->traceback.empty()) {
    joint_rep->traceback =
        r1->traceback.substr(0, s1_seq_len) +
        r2->traceback.substr(std::min((size_t)r2_skip, r2->traceback.size()));
  }

  // Check to see if we have to join scores
  if (r1->scores != nullptr && r2->scores != nullptr) {
    joint_rep->scores =
        (float *)malloc(sizeof(float) * joint_rep->repeatLength);

    // copy from repeat 1
    int i = 0;
    for (; i < s1_seq_len; ++i) {
      joint_rep->scores[i] = r1->scores[i];
    }

    // Copy from repeat 2
    for (int j = r2_skip; j < r2->repeatLength; ++j) {
      joint_rep->scores[i] = r2->scores[j];
      ++i;
    }
  }

  // Check to see if we have to join logo numbers
  if (r1->logoNumbers != nullptr && r2->logoNumbers != nullptr) {
    joint_rep->logoNumbers =
        (int *)malloc(sizeof(int) * joint_rep->repeatLength);

    // copy from repeat 1
    int i = 0;
    for (; i < s1_seq_len; ++i) {
      joint_rep->logoNumbers[i] = r1->logoNumbers[i];
    }

    // Copy from repeat 2
    for (int j = r2_skip; j < r2->repeatLength; ++j) {
      joint_rep->logoNumbers[i] = r2->logoNumbers[j];
      ++i;
    }
  }

  // We have to join score
  // Perform weighted average of scores
  float pct_seq1 = (float)s1_seq_len / (float)joint_rep->repeatLength;
  float pct_seq2 = 1.0 - pct_seq1;

  joint_rep->regionScore =
      (r1->regionScore * pct_seq1) + (r2->regionScore * pct_seq2);

  // When both regions were counted and scored, the joined region is r1
  // extended by the part of r2 past r1's end: its counts are r1's plus r2's
  // events there, and its score is r1's plus r2's per-position score changes
  // there. (The weighted average above, which stays for regions without
  // them, gave a repeat that crosses a window boundary about half its score
  // and no substitutions, insertions or deletions.)
  if (!r1->events.empty() && !r2->events.empty() && r2->scores != nullptr) {
    unsigned long r1End = r1->sequenceStart + r1->repeatLength;
    size_t past =
        r1End > r2->sequenceStart ? (size_t)(r1End - r2->sequenceStart) : 0;
    past = std::min(past, r2->events.size() / 3);
    joint_rep->events = r1->events;
    joint_rep->events.insert(joint_rep->events.end(),
                             r2->events.begin() + 3 * past, r2->events.end());
    unsigned int totals[3] = {0, 0, 0};
    for (size_t k = 0; k < joint_rep->events.size(); ++k)
      totals[k % 3] += joint_rep->events[k];
    joint_rep->mismatches = (int)totals[0];
    joint_rep->insertions = (int)totals[1];
    joint_rep->deletions = (int)totals[2];
    float score = r1->regionScore;
    for (size_t k = std::min(past, (size_t)r2->repeatLength);
         k < r2->repeatLength; ++k)
      score += r2->scores[k];
    joint_rep->regionScore = score;
  }

  // Join splits
  // Right now we do this lazily by combining split vectors
  // A better way of doing this is to join the logo nums
  // and then pass it through

  // r2 adds nothing, so the splits of r1 stand
  if (r2_inside) {
    if (r1->splits != nullptr) {
      joint_rep->splits = new std::vector<int>(*r1->splits);
      if (r1->consensi != nullptr)
        joint_rep->consensi = new std::vector<std::string>(*r1->consensi);
    }
    return joint_rep;
  }

  // don't bother joining repeat splits if diff repeat periods
  if (r1->repeatPeriod != r2->repeatPeriod)
    return joint_rep;
  // don't bother joining repeat splits if one of them doesn't have splits
  if (r1->splits == nullptr || r2->splits == nullptr)
    return joint_rep;

  joint_rep->splits = new std::vector<int>(*r1->splits);
  joint_rep->consensi = new std::vector<std::string>(*r1->consensi);

  int split_pos = r2->repeatLength;

  int split_start = -1;

  // find first split
  for (int i = 0; i < r2->splits->size(); ++i) {
    // Check if the split (should be) contained in r1
    if (r2->splits->at(i) < overlap)
      continue;

    else {
      split_start = i;
      break;
    }
  }

  if (split_start == -1)
    return joint_rep;

  for (int i = split_start; i < r2->splits->size(); ++i) {
    // Splits are positions within their own repeat; move r2's to where r2
    // sits in the joint region (-1 marks a split that was removed)
    int split_i = r2->splits->at(i);
    joint_rep->splits->push_back(split_i < 0 ? split_i : split_i + s1_seq_len);
    joint_rep->consensi->push_back(r2->consensi->at(i + 1));
  }

  return joint_rep;
}

std::string RepeatRegion::PermutationForString(const std::string &str,
                                               int offset) {
  std::string new_string;
  new_string.reserve(str.length());
  for (int i = 0; i < str.length(); ++i) {
    int p = (offset + i) % str.size();
    new_string.push_back(str[p]);
  }

  return new_string;
}

// Outputs 0 if equal, 1 if perm2<perm1, -1 if perm1<perm2
int RepeatRegion::CompareStrPerm(const std::string &str, int perm1, int perm2) {
  int size = str.length();
  for (int i = 0; i < size; ++i) {
    int p1 = (i + perm1) % size;
    int p2 = (i + perm2) % size;

    if (str[p1] < str[p2])
      return -1;
    if (str[p2] < str[p1])
      return 1;
  }

  return 0;
}

void RepeatRegion::SortConsensus(int index) {
  int best_perm = 0;
  for (int i = 1; i < consensi->at(index).length(); ++i) {
    int c = CompareStrPerm(consensi->at(index), best_perm, i);
    // If two permutations are identical, we can terminate
    if (c == 0)
      break;

    if (c == 1)
      best_perm = i;
  }

  if (best_perm == 0)
    return;
  this->consensi->at(index) =
      PermutationForString(this->consensi->at(index), best_perm);
}

void RepeatRegion::SortConsensi() {

  if (string_consensus.length() > 1) {
    int best_perm = 0;
    for (int i = 1; i < string_consensus.length(); ++i) {
      int c = CompareStrPerm(string_consensus, best_perm, i);
      // If two permutations are identical, we can terminate
      if (c == 0)
        break;

      if (c == 1)
        best_perm = i;
    }
    string_consensus = PermutationForString(string_consensus, best_perm);
  }

  if (this->consensi == nullptr) {
    return;
  }
  if (this->consensi->empty()) {
    return;
  }

  for (int i = 0; i < consensi->size(); ++i) {
    SortConsensus(i);
  }
}

RepeatRegion *dead_func(RepeatRegion *r1) {
  return joint_repeat_region(r1, r1);
}
