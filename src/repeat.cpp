//
//  repeat.cpp
//  ultraP
//

#include "repeat.hpp"
#include "RepeatSplitter.h"
#include <vector>
void RepeatRegion::CreateLogo(SequenceWindow *window, UMatrix *matrix) {
  logoMemory = (int *)malloc(sizeof(int) * repeatPeriod * (NUM_SYMBOLS + 1));
  logo = (int **)malloc(sizeof(int *) * repeatPeriod);

  for (int i = 0; i < repeatPeriod; ++i) {
    int p = i * NUM_SYMBOLS;
    logo[i] = &logoMemory[p];
    for (int j = 0; j < NUM_SYMBOLS; ++j) {
      logo[i][j] = 0;
    }
  }

  unsigned long backset = repeatPeriod;
  unsigned long endpoint = repeatLength - repeatPeriod;

  unsigned long fplace = repeatPeriod + windowStart;
  unsigned long bplace = windowStart;

  if (endpoint > window->length + window->overlap)
    endpoint = window->length + window->overlap;

  int insertionReset = -1;
  int deletionReset = -1;

  int p = 0;

  // float lastUnitScore = 0;
  int modp = 0;

  for (unsigned long i = 0; i < endpoint; ++i) {
    int frow = matrix->traceback[fplace];
    int brow = matrix->traceback[bplace];

    cell fdesc = matrix->cellDescriptions[frow];
    cell bdesc = matrix->cellDescriptions[brow];

    if (i < repeatPeriod)
      bdesc = matrix->cellDescriptions[0];

    if (fdesc.type == CT_INSERTION) {
      insertions += fdesc.indelNumber;
      i += fdesc.indelNumber;
      fplace += fdesc.indelNumber;
      /*     for (int j = 0; j < fdesc.indelNumber; ++j) {
               lookBack[i + j] = 0;
           }*/

    }

    else if (fdesc.type == CT_DELETION) {
      deletions += fdesc.indelNumber;
      bplace += fdesc.indelNumber;
      p += fdesc.indelNumber;
      modp += fdesc.indelNumber + 1;
    }

    else {
      ++modp;
    }

    modp = modp % repeatPeriod;

    if (bdesc.type == CT_INSERTION) {
      bplace += bdesc.indelNumber;
    }

    else if (bdesc.type == CT_DELETION) {
      fplace += bdesc.indelNumber;
      i += bdesc.indelNumber;
    }

    symbol s = window->seq[fplace];
    symbol sb = window->seq[bplace];

    p = p % repeatPeriod;

    if (s != sb) {
      mismatches += 1;
      logo[p][s] += 1;
    }

    else {
      logo[p][s] += 4;
    }

    ++p;
    ++fplace;
    ++bplace;
    --insertionReset;
    --deletionReset;

    if (deletionReset == 0) {
      backset = backset + repeatPeriod;
    }

    if (insertionReset == 0) {
      backset = repeatPeriod;
    }
  }
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

  logPVal = 0;

  winOverlapSize = 0;
  overlapCorrection = OC_NONE;

  scores = NULL;
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
    float s = matrix->scoreColumns[p][cell];

    int prevCell = matrix->traceback[p - 1];
    if (matrix->cellDescriptions[prevCell].type != CT_MATCH)
      prevCell = matrix->cellDescriptions[cell].parentIndex;
    float ps = matrix->scoreColumns[p - 1][prevCell];

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

  RepeatRegion *region = NULL;

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

    region->startScore = matrix->scoreColumns[i][row];

    int length = region->repeatPeriod;

    for (; i < seqLength; ++i) {
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

    region->endScore = matrix->scoreColumns[i - 1][row];
    region->regionScore = region->endScore - region->startScore;
  }

  *pos = i;

  if (region != NULL) {
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
  joint_rep->repeatLength =
      (r2->sequenceStart + r2->repeatLength) - r1->sequenceStart;

  int s1_seq_len = r2->sequenceStart - r1->sequenceStart;
  int overlap = (r1->sequenceStart + r1->repeatLength) - r2->sequenceStart;

  joint_rep->sequenceName = r1->sequenceName;
  joint_rep->sequenceID = r1->sequenceID;
  joint_rep->readID = r2->readID;
  joint_rep->winOverlapSize = r1->winOverlapSize;
  joint_rep->winTotalLength = r1->winTotalLength;

  joint_rep->combinedRepeat = true;

  joint_rep->repeatPeriod = r1->repeatPeriod;

  joint_rep->overlapCorrection = OC_TRUE;

  // We have to join ->sequence
  if (!r1->sequence.empty() && !r2->sequence.empty()) {
    joint_rep->sequence = r1->sequence.substr(0, s1_seq_len) + r2->sequence;
  }
  // We have to join score
  // Perform weighted average of scores
  float pct_seq1 = (float)s1_seq_len / (float)joint_rep->repeatLength;
  float pct_seq2 = 1.0 - pct_seq1;

  joint_rep->regionScore =
      (r1->regionScore * pct_seq1) + (r2->regionScore * pct_seq2);

  // Join splits
  // Right now we do this lazily by combining split vectors
  // A better way of doing this is to join the logo nums
  // and then pass it through

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
    joint_rep->splits->push_back(r2->splits->at(i));
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
