//
// Created by Olson, Daniel (NIH/NIAID) [E] on 2/2/23.
//

#ifndef ULTRA_MASK_H
#define ULTRA_MASK_H
#include <string>
#include <vector>
#include <stdio.h>
#include <unordered_map>
#include <sstream>
#include <fstream>
typedef struct s_mask_region {
  unsigned long long start;
  unsigned long long end;
} mregion;

void OutputMaskedFASTA(std::string in_path,
                       FILE* out_file,
                       const std::unordered_map<std::string,
                                                std::vector<mregion> *> &masks);

#endif // ULTRA_MASK_H
