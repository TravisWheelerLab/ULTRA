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

  bool operator()(s_mask_region m1, s_mask_region m2);
} mregion;

void OutputCaseMaskedFASTA(std::string in_path,
                       FILE* out_file,
                        std::unordered_map<std::string,
                                                std::vector<mregion> *> masks);
void OutputNMaskedFASTA(std::string in_path,
                           FILE* out_file,
                            std::unordered_map<std::string,
                                                    std::vector<mregion> *> masks);
void OutputRemoveMaskedFASTA(std::string in_path,
                           FILE* out_file,
                            std::unordered_map<std::string,
                                                    std::vector<mregion> *> masks);
std::vector<mregion> *CleanedMasks(std::vector<mregion> *mask);



#endif // ULTRA_MASK_H
