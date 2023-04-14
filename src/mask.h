

#ifndef ULTRA_MASK_H
#define ULTRA_MASK_H
#include <string>
#include <vector>
#include <cstdio>
#include <unordered_map>
#include <sstream>
#include <fstream>

typedef struct s_mask_region {
  unsigned long long start;
  unsigned long long end;

  bool operator()(s_mask_region m1, s_mask_region m2);
} mregion;

void OutputMaskedFASTA(const std::string &in_path,
                       FILE* out_file,
                       std::unordered_map<unsigned long long, std::vector<mregion> *> masks,
                       bool n_mask);

std::vector<mregion> *CleanedMasks(std::vector<mregion> *mask);



#endif // ULTRA_MASK_H
