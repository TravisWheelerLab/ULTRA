//
// Created by Daniel Olson on 10/20/21.
//

#include "repeatSplitter.hpp"
#include "repeat.hpp"

#define NUMSYM (4)

// Note, it is necessary to give each PairWindow a partner,
// and it is also necessary to call calculateTerms() on each window
PairWindow::PairWindow(int period, float pseudo_count) {
  this->period = period;

  this->start = 0;
  this->end = 0;

  this->pair = nullptr;

  this->sum = 0;
  this->plogp = 0;
  this->plogq = 0;

  this->logo = (float *)malloc(((period * NUMSYM) + 1) * sizeof(float));

  this->fillWithPseudoCount(pseudo_count);
}

void PairWindow::fillWithPseudoCount(float pseudo_count) {
  int default_value = pseudo_count / NUMSYM;
  this->sum = 0.0;
  for (int i = 0; i < period; ++i) {
    for (int j = 0; j < NUMSYM; ++j) {
      this->logo[(i * NUMSYM) + j] = default_value;
    }
    this->sum += pseudo_count;
  }
}

void PairWindow::calculateTerms() {
  assert(this->pair != nullptr); // I need to use more asserts...

  this->plogp = 0;
  this->plogq = 0;

  for (int i = 0; i < period * NUMSYM; ++i) {
    float p = this->logo[i] / this->sum;
    float q = this->pair->logo[i] / this->pair->sum;

    this->plogp += p * log2(p);
    this->plogq += p * log2(q);
  }
}

float PairWindow::KLD() {
  float kld = (this->plogp - this->plogq) / this->sum;
  kld += log2(this->pair->sum) - log2(this->sum);
  return kld;
}

float PairWindow::slowKLD() {
  float self_sum = 0.0;
  float partner_sum = 0.0;

  // Take initial sum
  for (int i = 0; i < this->period; ++i) {
    for (int j = 0; j < NUMSYM; ++j) {
      self_sum += logo[(i * NUMSYM) + j];
      partner_sum += this->pair->logo[(i * NUMSYM) + j];
    }
  }

  float kld = 0.0;

  // Calculate KLD
  for (int i = 0; i < this->period; ++i) {
    for (int j = 0; j < NUMSYM; ++j) {
      float p = logo[(i * NUMSYM) + j] / self_sum;
      float q = this->pair->logo[(i * NUMSYM) + j] / partner_sum;

      kld += p * log2(p / q);
    }
  }

  return kld;
}

float PairWindow::splitValue() {
  float kld = this->KLD();
  return 1.0 - exp2(-kld);
}

void PairWindow::editCell(int i, int j, float delta) {
  // float oldP = this->logo[(i * NUMSYM) + j] / this->sum;
  // float newP = (this->logo[(i * NUMSYM) + j] + delta) / (this->sum + delta);

  // this->sum += delta;
  printf("Not implemented\n");
  exit(-1);
}

// There is some room to squeeze a few cycles of performance
// Not worth it right now...
// Maybe in the future?
void PairWindow::selfishEditCell(int i, int j, float delta) {

  float q = this->pair->logo[(i * NUMSYM) + j];
  float oldP = this->logo[(i * NUMSYM) + j] / this->sum;

  this->sum += delta; // Adjust the sum appropriately
  this->logo[(i * NUMSYM) + j] += delta;

  float newP = this->logo[(i * NUMSYM) + j] / this->sum;

  // Remove the old value from our sums
  this->plogp -= oldP * log2(oldP);
  this->plogq -= oldP * log2(q);

  // Add the new value to our sums
  this->plogp += newP * log2(newP);
  this->plogq += newP * log2(q);
}

void PairWindow::dumbEditCell(int i, int j, float delta) {
  this->logo[(i * NUMSYM) + j] += delta;
  this->sum += delta;
}

void PairWindow::dumbMoveWindowForward(RepeatRegion *r) {
  int s = SymbolForChar(r->sequence[this->start]);
  int p = r->logoNumbers[this->start];

  if (s == N_UNK) {
    this->logo[(p * NUMSYM) + 0] -= 0.25;
    this->logo[(p * NUMSYM) + 1] -= 0.25;
    this->logo[(p * NUMSYM) + 2] -= 0.25;
    this->logo[(p * NUMSYM) + 3] -= 0.25;
    this->sum -= 1.0;
  } else {
    assert(s == N_A || s == N_T || s == N_C || s == N_G);
    s = s - 1; // N_A starts at 1, but we want indexing starting at 0

    this->logo[(p * NUMSYM) + s] -= 1.0;
    this->sum -= 1.0;
  }

  this->start += 1;
  s = SymbolForChar(r->sequence[this->end]);
  p = r->logoNumbers[this->end];

  if (s == N_UNK) {
    this->logo[(p * NUMSYM) + 0] -= 0.25;
    this->logo[(p * NUMSYM) + 1] -= 0.25;
    this->logo[(p * NUMSYM) + 2] -= 0.25;
    this->logo[(p * NUMSYM) + 3] -= 0.25;
    this->sum -= 1.0;
  } else {
    assert(s == N_A || s == N_T || s == N_C || s == N_G);
    s = s - 1; // N_A starts at 1, but we want indexing starting at 0

    this->logo[(p * NUMSYM) + s] -= 1.0;
    this->sum -= 1.0;
  }

  this->end += 1;
}

void PairWindow::moveWindowForward(RepeatRegion *r) {
  int s = SymbolForChar(r->sequence[this->start]);
  int p = r->logoNumbers[this->start];

  if (s == N_UNK) {
    this->editCell(p, 0, -0.25);
    this->editCell(p, 1, -0.25);
    this->editCell(p, 2, -0.25);
    this->editCell(p, 3, -0.25);
  } else {
    assert(s == N_A || s == N_T || s == N_C || s == N_G);
    s = s - 1; // N_A starts at 1, but we want indexing starting at 0
    this->editCell(p, s, -1.0);
  }
  this->start += 1;
  s = SymbolForChar(r->sequence[this->end]);
  p = r->logoNumbers[this->end];

  if (s == N_UNK) {
    this->editCell(p, 0, -0.25);
    this->editCell(p, 1, -0.25);
    this->editCell(p, 2, -0.25);
    this->editCell(p, 3, -0.25);
  } else {
    assert(s == N_A || s == N_T || s == N_C || s == N_G);
    s = s - 1; // N_A starts at 1, but we want indexing starting at 0
    this->editCell(p, s, -1.0);
  }

  this->end += 1;
}

void PairWindow::fillWindow(RepeatRegion *r, int start, int end) {
  for (int i = start; i < end; ++i) {
    int p = r->logoNumbers[i];
    if (p < 0)
      continue;

    int s = SymbolForChar(r->sequence[i]);
    if (s == N_UNK) {
      this->logo[(p * NUMSYM) + 0] += 0.25;
      this->logo[(p * NUMSYM) + 1] += 0.25;
      this->logo[(p * NUMSYM) + 2] += 0.25;
      this->logo[(p * NUMSYM) + 3] += 0.25;
      this->sum += 1.0;
    }

    else {
      assert(s == N_A || s == N_T || s == N_C || s == N_G);
      s = s - 1; // N_A starts at 1, but we want indexing starting at 0

      this->logo[(p * NUMSYM) + s] += 1.0;
      this->sum += 1.0;
    }
  }
  this->start = start;
  this->end = end;
}

std::vector<int> *SplitRepeat(RepeatRegion *r, float threshold, int windowUnits,
                              int minSize, int minLagtime) {

  int windowSize = windowUnits * r->repeatPeriod;
  if (windowSize < minSize) {
    windowSize = minSize + (minSize % r->repeatPeriod);
  }

  int lagtime = r->repeatPeriod;
  if (lagtime < minLagtime) {
    lagtime = minLagtime;
  }

  assert(windowSize * 2 < r->repeatLength);
  assert(r->logoNumbers == nullptr);
  assert(r->sequence.size() > 0);

  std::vector<int> *splits = new std::vector<int>();

  PairWindow *leftWindow = new PairWindow(r->repeatPeriod);
  PairWindow *rightWindow = new PairWindow(r->repeatPeriod);

  leftWindow->pair = rightWindow;
  rightWindow->pair = leftWindow;

  leftWindow->fillWindow(r, 0, windowSize);
  rightWindow->fillWindow(r, windowSize, 2 * windowSize);

  leftWindow->calculateTerms();
  rightWindow->calculateTerms();

  int splitPos = -1;
  float splitValue = leftWindow->splitValue();

  if (splitValue >= threshold) {
    splitPos = leftWindow->end;
  }

  int seqLen = r->sequence.size();

  while (rightWindow->end < seqLen) {

    // Move right window forward
    rightWindow->dumbMoveWindowForward(r);
    // Move left window forward
    leftWindow->moveWindowForward(r);

    // Determine and record if we  have a new split value peak
    float tmp = leftWindow->splitValue();
    if (tmp > threshold) {
      if (tmp > splitValue) {
        splitValue = tmp;
        splitPos = leftWindow->end;
      }
    }

    if (splitPos >= 0) {
      if (leftWindow->end - splitPos > lagtime) {

        splits->push_back(splitPos);
        splitPos = -1;
        splitValue = 0;
        // Move windows forward
        leftWindow->fillWithPseudoCount(1.0);
        rightWindow->fillWithPseudoCount(1.0);
        leftWindow->start = splitPos;
        leftWindow->end = splitPos + windowSize;
        rightWindow->start = splitPos + windowSize;
        rightWindow->end = splitPos + windowSize + windowSize;

        if (rightWindow->end >= seqLen) {
          break;
        }

        else {
          leftWindow->fillWindow(r, leftWindow->start, leftWindow->end);
          rightWindow->fillWindow(r, rightWindow->start, rightWindow->end);

          leftWindow->calculateTerms();
          rightWindow->calculateTerms();

          tmp = leftWindow->splitValue();
          if (tmp > threshold) {
            if (tmp > splitValue) {
              splitValue = tmp;
              splitPos = leftWindow->end;
            }
          }
        }
      }
    }
  }

  if (splitPos >= 0) {
    splits->push_back(splitPos);
  }

  return splits;
}

std::vector<std::string> *ConsensusForSplits(RepeatRegion *r,
                                             std::vector<int> *splits,
                                             float consensusThreshold) {

  std::vector<std::string> *consensi = new std::vector<std::string>();

  return consensi;
}

void FilterSplits(std::vector<int> *splits, std::vector<std::string> *consensus,
                  float threshold, float wildstar_weight) {}
