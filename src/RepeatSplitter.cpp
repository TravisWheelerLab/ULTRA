//
// Created by Olson, Daniel (NIH/NIAID) [E] on 4/28/22.
//

#include "RepeatSplitter.hpp"

#define score_func() (KLD_pq() + KLD_qp())

void SplitWindow::AllocateSplitWindow(int alphabet, int max_period) {
  alphabetSize = alphabet;
  maxPeriod = max_period;

  cp = (float **)malloc(sizeof(float *) * alphabetSize);
  cq = (float **)malloc(sizeof(float *) * alphabetSize);
  logo = (float **)malloc(sizeof(float *) * alphabetSize);

  logo_totals = (float *)malloc(sizeof(float) * max_period);
  p_mat = (float *)malloc(sizeof(float) * maxPeriod * alphabetSize);
  q_mat = (float *)malloc(sizeof(float) * maxPeriod * alphabetSize);
  logo_mat = (float *)malloc(sizeof(float) * alphabetSize * maxPeriod);

  for (int i = 0; i < alphabetSize; ++i) {
    logo[i] = &logo_mat[i * maxPeriod];
    cp[i] = &p_mat[i * maxPeriod];
    cq[i] = &q_mat[i * maxPeriod];
  }
}

void SplitWindow::DeallocSplitWindow() {
  free(cp);
  cp = nullptr;
  free(cq);
  cq = nullptr;
  free(logo);
  logo = nullptr;

  free(logo_totals);
  logo_totals = nullptr;
  free(p_mat);
  p_mat = nullptr;
  free(q_mat);
  q_mat = nullptr;
  free(logo_mat);
  logo_mat = nullptr;
}

void SplitWindow::StartSplitWindow(int p, float startingScore) {
  period = p;
  place = 0;

  float d = startingScore / ((float)alphabetSize * (float)period);

  plogp = d * log2(d) * (float)period * (float)alphabetSize;
  plogq = plogp;
  qlogq = plogp;
  qlogp = plogp;

  Np = startingScore;
  Nq = startingScore;

  logNp = log2(Np);
  logNq = log2(Nq);

  for (int a = 0; a < alphabetSize; ++a) {
    for (int i = 0; i < period; ++i) {
      cp[a][i] = d;
      cq[a][i] = d;
    }
  }
}

float SplitWindow::slow_KLD_pq() {
  float sum = 0;

  for (int a = 0; a < alphabetSize; ++a) {
    for (int i = 0; i < period; ++i) {
      float p = cp[a][i] / Np;
      float q = cq[a][i] / Nq;

      sum += p * log2(p / q);
    }
  }

  return sum;
}

float SplitWindow::slow_KLD_qp() {
  float sum = 0;

  for (int a = 0; a < alphabetSize; ++a) {
    for (int i = 0; i < period; ++i) {
      float p = cp[a][i] / Np;
      float q = cq[a][i] / Nq;

      sum += q * log2(q / p);
    }
  }
  return sum;
}

float SplitWindow::KLD_pq() {
  return (plogp / Np) - (plogq / Np) + logNq - logNp;
}

float SplitWindow::KLD_qp() {
  return (qlogq / Nq) - (qlogp / Np) + logNp - logNq;
}

float SplitWindow::AdjustMatrix(bool is_p, int i, int symbol, float delta) {
  float *aloga;
  float *alogb;
  float *bloga;
  float **ca;
  float **cb;

  if (i < 0)
    return 0;

  if (is_p) {
    aloga = &plogp;
    alogb = &plogq;
    bloga = &qlogp;
    ca = cp;
    cb = cq;
  }

  else {
    aloga = &qlogq;
    alogb = &qlogp;
    bloga = &plogq;
    ca = cq;
    cb = cp;
  }

  if (symbol < 0) {
    float d = delta / (float)alphabetSize;
    for (int a = 0; a < alphabetSize; ++a) {
      *aloga -= ca[a][i] * log2(ca[a][i]);
      *alogb -= ca[a][i] * log2(cb[a][i]);
      *bloga -= cb[a][i] * log2(ca[a][i]);

      ca[a][i] += d;

      *aloga += ca[a][i] * log2(ca[a][i]);
      *alogb += ca[a][i] * log2(cb[a][i]);
      *bloga += cb[a][i] * log2(ca[a][i]);
    }
  }

  else {
    *aloga -= ca[symbol][i] * log2(ca[symbol][i]);
    *alogb -= ca[symbol][i] * log2(cb[symbol][i]);
    *bloga -= cb[symbol][i] * log2(ca[symbol][i]);

    ca[symbol][i] += delta;

    *aloga += ca[symbol][i] * log2(ca[symbol][i]);
    *alogb += ca[symbol][i] * log2(cb[symbol][i]);
    *bloga += cb[symbol][i] * log2(ca[symbol][i]);
  }

  return delta;
}

void SplitWindow::FillWindow(RepeatRegion *r, int windowSize) {
  StartSplitWindow(r->repeatPeriod, 1.0);
  length = windowSize;

  // Fill left window
  for (int i = 0; i < windowSize; ++i) {
    int symbol = SymbolForChar(r->sequence[i]) - 1;
    int indx = r->logoNumbers[i];

    if (indx < 0)
      continue;

    if (symbol < 0) {
      float d = 1.0 / alphabetSize;

      for (int a = 0; a < alphabetSize; ++a) {
        cp[a][indx] += d;
      }
    }

    else {
      cp[symbol][indx] += 1.0;
    }
  }

  // Fill right window
  for (int i = length; i < 2 * length; ++i) {
    int symbol = SymbolForChar(r->sequence[i]) - 1;
    int indx = r->logoNumbers[i];

    if (indx < 0)
      continue;
    if (symbol < 0) {
      float d = 1.0 / alphabetSize;

      for (int a = 0; a < alphabetSize; ++a) {
        cq[a][indx] += d;
      }
    }

    else {
      cq[symbol][indx] += 1.0;
    }
  }

  // Count matrix
  plogp = 0;
  plogq = 0;
  qlogq = 0;
  qlogp = 0;
  Np = 0;
  Nq = 0;

  for (int i = 0; i < period; ++i) {
    for (int a = 0; a < alphabetSize; ++a) {
      plogp += cp[a][i] * log2(cp[a][i]);
      plogq += cp[a][i] * log2(cq[a][i]);
      qlogp += cq[a][i] * log2(cp[a][i]);
      qlogq += cq[a][i] * log2(cq[a][i]);
      Np += cp[a][i];
      Nq += cq[a][i];
    }
  }

  logNp = log2(Np);
  logNq = log2(Nq);
}

void SplitWindow::MoveWindowForward(RepeatRegion *r) {

  int old_left_symbol = SymbolForChar(r->sequence[place]) - 1;
  int old_right_symbol = SymbolForChar(r->sequence[place + length]) - 1;
  int new_left_symbol = old_right_symbol;
  int new_right_symbol =
      SymbolForChar(r->sequence[place + length + length]) - 1;

  int old_left_indx = r->logoNumbers[place];
  int old_right_indx = r->logoNumbers[place + length];
  int new_left_indx = old_right_indx;
  int new_right_indx = r->logoNumbers[place + length + length];

  Np += AdjustMatrix(true, old_left_indx, old_left_symbol, -1.0);
  Nq += AdjustMatrix(false, old_right_indx, old_right_symbol, -1.0);

  Np += AdjustMatrix(true, new_left_indx, new_left_symbol, 1.0);
  Nq += AdjustMatrix(false, new_right_indx, new_right_symbol, 1.0);

  logNp = log2(Np);
  logNq = log2(Nq);

  place += 1;
}

std::vector<float> *SplitWindow::RegionScores(RepeatRegion *r, int windowSize) {

  std::vector<float> *scores = new std::vector<float>();

  FillWindow(r, windowSize);

  for (int i = 0; i < length; ++i) {
    scores->push_back(0.0);
  }
  scores->push_back(score_func());

  // exit(0);

  for (int i = length; i < r->sequence.size() - length; ++i) {
    MoveWindowForward(r);
    scores->push_back(score_func());
  }

  for (int i = 0; i < length; ++i) {
    scores->push_back(0);
  }
  return scores;
}

std::vector<int> *SplitWindow::SplitsForRegion(RepeatRegion *r, int window_size,
                                               float threshold) {

  std::vector<int> *splits = new std::vector<int>();

  int split_pos = -1;
  float split_val = threshold;

  // printf("Filling window\n");
  FillWindow(r, window_size);

  float score = score_func();

  // printf("Splitting region %i %f\n", window_size, threshold);

  if (score > split_val) {
    split_val = score;
    split_pos = place + length;
  }

  // printf("Moving window forward\n");
  for (int i = 0; i < r->sequence.size() - length - length; ++i) {
    // printf("%i %i %i\n", i, place, length);
    MoveWindowForward(r);
    score = score_func();
    // printf("%i %f: ", place + length, score, threshold);
    if (score > split_val) {
      //   printf("split");
      split_val = score;
      split_pos = place + length;
    }

    else if (split_pos >= 0 && place - split_pos > 0) {
      // printf(" final %i", split_pos);
      splits->push_back(split_pos);
      split_pos = -1;
      split_val = threshold;
    }
    //  printf("\n");
  }
  // printf("finishing\n");

  if (split_pos >= 0) {
    // printf("split terminal\n");
    splits->push_back(split_pos);
    split_pos = -1;
    split_val = threshold;
  }

  return splits;
}

void SplitWindow::FillLogoForRegion(RepeatRegion *r, int start, int end) {
  float d = 1.0 / (float)alphabetSize;

  for (int i = 0; i < r->repeatPeriod; ++i) {
    for (int a = 0; a < alphabetSize; ++a) {
      logo[a][i] = d;
    }
    logo_totals[i] = 0;
  }

  for (int i = start; i < end; ++i) {
    int symbol = SymbolForChar(r->sequence[i]) - 1;
    int indx = r->logoNumbers[i];

    if (indx < 0 || symbol < 0)
      continue;

    logo[symbol][indx] += 1.0;
    logo_totals[indx] += 1.0;
  }
}

std::string SplitWindow::ConsensusForRegion(RepeatRegion *r, int start, int end,
                                            float consensus_threshold) {
  std::string consensus = std::string(r->repeatPeriod, '*');

  FillLogoForRegion(r, start, end);

  for (int i = 0; i < r->repeatPeriod; ++i) {
    float bestVal = consensus_threshold;
    for (int a = 0; a < alphabetSize; ++a) {
      float mass = logo[a][i] / logo_totals[i];
      if (mass > bestVal) {
        bestVal = mass;
        consensus[i] = CharForSymbol(a + 1);
      }
    }
  }

  return consensus;
}

std::vector<std::string> *
SplitWindow::ConsensiForSplit(RepeatRegion *r, std::vector<int> *splits,
                              float consensus_threshold) {
  std::vector<std::string> *consensi = new std::vector<std::string>();

  int start = 0;
  for (int i = 0; i < splits->size(); ++i) {
    consensi->push_back(
        ConsensusForRegion(r, start, splits->at(i), consensus_threshold));
    start = splits->at(i);
  }

  consensi->push_back(
      ConsensusForRegion(r, start, r->sequence.size(), consensus_threshold));

  return consensi;
}

float ConsensusSimilarity(std::string *c1, std::string *c2, int offset) {
  float similarity = 0.0;
  for (int i = 0; i < c1->size(); ++i) {
    int offi = (i + offset) % c2->size();
    if (c1->at(i) == c2->at(offi))
      similarity += 1.0;

    else if (c1->at(i) == '*' || c1->at(offi) == '*')
      similarity += 0.25;
  }

  return similarity / c1->size();
}

bool ShouldJoinConsensus(std::string *c1, std::string *c2,
                         float join_threshold) {
  for (int i = 0; i < c1->size(); ++i) {
    if (ConsensusSimilarity(c1, c2, i) > join_threshold)
      return true;
  }

  return false;
}

void ValidateSplits(std::vector<std::string> *consensi,
                    std::vector<int> *splits, float join_threshold) {

  for (int i = 1; i < consensi->size(); ++i) {
    if (ShouldJoinConsensus(&consensi->at(i - 1), &consensi->at(i),
                            join_threshold)) {
      consensi->at(i) = consensi->at(i - 1);
      splits->at(i - 1) = -1;
    }
  }
}