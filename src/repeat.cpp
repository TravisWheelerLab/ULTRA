//
//  repeat.cpp
//  ultraP
//

#include "repeat.hpp"
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

  // printf("%llu %llu %llu\n", backset, fplace, endpoint);
  if (endpoint > window->length + window->overlap)
    endpoint = window->length + window->overlap;

  // printf("\n%i vs %i\n", windowStart, endpoint);

  int insertionReset = -1;
  int deletionReset = -1;

  int p = 0;

  // double lastUnitScore = 0;
  int modp = 0;

  for (unsigned long i = 0; i < endpoint; ++i) {
    int frow = matrix->traceback[fplace];
    int brow = matrix->traceback[bplace];
    //   printf("%f\n", matrix->scoreColumns[fplace][frow]);
    // printf("%i(%i,%i,%i).", row, insertions, deletions, mismatches);

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

      // printf("Insertion.\n");
    }

    else if (fdesc.type == CT_DELETION) {
      deletions += fdesc.indelNumber;
      bplace += fdesc.indelNumber;
      p += fdesc.indelNumber;
      modp += fdesc.indelNumber + 1;
      // printf("Deletion.\n");
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
    //  printf("%i (%i): %c %c\n", fplace, bplace, CharForSymbol(s),
    //  CharForSymbol(sb));
    p = p % repeatPeriod;

    if (s != sb) {
      mismatches += 1;
      logo[p][s] += 1;
    }

    else {
      // printf("\n**%i %i**\n", p, s);
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

  if (sequence.length() == 0 || traceback.length() == 0)
    return;
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
      //   printf("%i ", logo[i][j]);
    }
    //  printf("\n");
    consensus[i] = s;
  }
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

  for (int i = 0; i < repeatPeriod; ++i) {
    con.push_back(CharForSymbol(consensus[i]));
  }
  con[repeatPeriod] = '\0';

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

  scores = (double *)malloc(sizeof(double) * repeatLength);

  for (unsigned long i = 1; i <= repeatLength; ++i) {

    unsigned long p = i + windowStart;
    int cell = matrix->traceback[p];
    if (matrix->cellDescriptions[cell].type != CT_MATCH)
      cell = matrix->cellDescriptions[cell].parentIndex;
    double s = matrix->scoreColumns[p][cell];

    int prevCell = matrix->traceback[p - 1];
    if (matrix->cellDescriptions[prevCell].type != CT_MATCH)
      prevCell = matrix->cellDescriptions[cell].parentIndex;
    double ps = matrix->scoreColumns[p - 1][prevCell];

    scores[i - 1] = s - ps;

    //    printf("[%i %i]: (%f, %f) %i %f %i %f \n", i, p,
    //    matrix->scoreColumns[p - 1][27], matrix->scoreColumns[p][27],
    //    cell, s, prevCell, ps);
  }
  // exit(0);
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

  /*printf("\n------SEQUENCE------\n");

  for (int i = 0; i < repeatLength; ++i) {
      printf("%c", sequence[i]);
  }

  printf("\n------END SEQUENCE------\n");

  printf("\n------TRACEBACK------\n");

  for (int i = 0; i < repeatLength; ++i) {
      printf("%c", traceback[i]);
  }

  printf("\n------END TRACEBACK------\n");


  printf("\n------LOOK BACK------\n");

  for (int i = 0; i < repeatLength; ++i) {
      printf("%i", lookBack[i]);
  }

  printf("\n------END LOOK BACK------\n");

  printf("\n------LOOK FORWARD------\n");

  for (int i = 0; i < repeatLength; ++i) {
      printf("%i", lookForward[i]);
  }

  printf("\n------END LOOK FORWARD------\n");*/

  backCounts = (char *)malloc(sizeof(char) * repeatLength);
  forwardCounts = (char *)malloc(sizeof(char) * repeatLength);

  CountFromDistance(lookBack, backCounts, -1, depth);
  CountFromDistance(lookForward, forwardCounts, 1, depth);

  /* printf("\n------BACK COUNTS------\n");

   for (int i = 0; i < repeatLength; ++i) {
       int c = backCounts[i];
       if (c < 0)
           printf("*");
       else
           printf("%i", backCounts[i]);
   }

   printf("\n------END BACK COUNTS------\n");

   printf("\n------FORWARD COUNTS------\n");

   for (int i = 0; i < repeatLength; ++i) {
       int c = forwardCounts[i];
       if (c < 0)
           printf("*");
       else
           printf("%i", forwardCounts[i]);
   }

   printf("\n------END FORWARD COUNTS------\n");*/

  // CountFromDistance(lookForward, forwardCounts, 1, depth);
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
      // printf("i:%i f:%i bf:%i ff:%i")

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
  double scorePerSymbol = regionScore / ((double)repeatLength);

  region->startScore = startScore + (scorePerSymbol * ((double)start));
  region->endScore = region->startScore + (scorePerSymbol * ((double)length));
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

    // printf("%i: %llx, %llx, %llx\n", i, (unsigned long)matrix, (unsigned
    // long)matrix->traceback, (unsigned long)&matrix->traceback[i]);
    if (matrix->traceback[i] > 0) {
      foundRepeat = true;
      break;
    }

    else if (matrix->traceback[i] < 0) {
      //     printf("%i %i\n", i, matrix->traceback[i]);
    }
  }

  if (foundRepeat) {
    int row = matrix->traceback[i];
    // printf("%i %i\n", matrix->traceback[i], matrix->traceback[i-1]);
    cell desc = matrix->cellDescriptions[row];

    region = new RepeatRegion();
    region->winOverlapSize = window->overlap;
    region->winTotalLength = window->length + window->overlap;
    region->sequenceName = window->sequenceName;
    region->sequenceID = window->seqID;
    region->readID = window->readID;
    region->repeatPeriod = desc.order;

    // printf("%i %i %i %i\n", window->start, i, window->overlap,
    // region->repeatPeriod);
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

  // for (int i = 0; i < )
  // printf("\n\n\n\n");
  // printf("Get logo numbers for %i %i\n", start, length);
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

      /*printf("D:%i %i %c %c %i %i\n", i, desc.type, sequence[i], traceback[i],
       * matrix->traceback[p], logoNumbers[i]);*/
      ++num;
      ++i;
    }
    if (delCount > 0)
      --i;
    /*
     if (delCount == 0)
     printf("%i %i %c %c %i %i\n", i, desc.type, sequence[i], traceback[i],
     matrix->traceback[p], logoNumbers[i]);*/
  }
}

bool repeats_overlap(RepeatRegion *r1, RepeatRegion *r2,
                     bool require_same_period) {
  if (r1->sequenceID == r2->sequenceID) {
    if (require_same_period) {
      if (r1->repeatPeriod == r2->repeatPeriod) {
        return (r1->sequenceStart + r1->repeatLength >= r2->sequenceStart);
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
  joint_rep->sequenceStart = r1->sequenceStart;
  joint_rep->repeatLength = (r2->sequenceStart + r2->repeatLength) - r1->sequenceStart;

  int s1_seq_len = r2->sequenceStart - r1->sequenceStart;

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

  joint_rep->regionScore = (r1->regionScore * pct_seq1) + (r2->regionScore * pct_seq2);

  // If possible we should join subrepeats

  return joint_rep;
}

RepeatRegion * dead_func(RepeatRegion *r1) {
  return joint_repeat_region(r1, r1);
}
