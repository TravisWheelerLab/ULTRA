//
// Created by Olson, Daniel (NIH/NIAID) [E] on 2/2/23.
//

#include "mask.h"

void OutputMaskedFASTA(std::string in_path,
                       FILE* out_file,
                        std::unordered_map<std::string,
                                                std::vector<mregion> *> masks) {

  std::ifstream in_file(in_path);
  if (!in_file.is_open()) {
    printf("Unable to open %s\n", in_path.c_str());
    return;
  }

  std::string line;

  std::string seq_name = "";

  unsigned long long seq_pos = 0;

  std::vector <mregion> *cmask = nullptr;

  while(std::getline(in_file, line)) {

    if (line[0] == '>') {
      seq_name = line.substr(1, std::string::npos);
      fprintf(out_file, "%s", line.c_str());

      if (masks.find(seq_name) != masks.end()) {
        delete cmask;
        cmask = CleanedMasks(masks[seq_name]);
      }
    }

    else {
      if (cmask == nullptr) {
        fprintf(out_file, "%s", line.c_str());
        continue;
      }
    }

  }

  in_file.close();
}

bool mregion::operator()(mregion m1, mregion m2) {
  if (m1.start != m2.start)
    return (m1.start > m2.start);

  else
    return m1.end > m2.end;
}

std::vector<mregion> *CleanedMasks(std::vector<mregion> *mask) {
  if (mask->empty())
    return nullptr;

  std::sort(mask->begin(), mask->end(), mregion());

  // combine overlapped regions
  std::vector<mregion> *cmask = new std::vector<mregion>();
  cmask->push_back(mask->at(0));
  for (int i = 1; i < mask->size(); ++i) {
    if (mask->at(i).start <= cmask->back().end) {
      if (mask->at(i).end > cmask->back().end)
        cmask->back().end = mask->at(i).end;
    }

    else {
      cmask->push_back(mask->at(i));
    }
  }

  return cmask;
}
