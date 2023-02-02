//
// Created by Olson, Daniel (NIH/NIAID) [E] on 2/2/23.
//

#include "mask.h"

void OutputMaskedFASTA(std::string in_path,
                       FILE* out_file,
                       const std::unordered_map<std::string,
                                                std::vector<mregion> *> &masks) {

  std::ifstream in_file(in_path);
  if (!in_file.is_open()) {
    printf("Unable to open %s\n", in_path.c_str());
    return;
  }

  std::string line;

  std::string sequence_name

  while(std::getline(in_file, line)) {

  }

  in_file.close();
}