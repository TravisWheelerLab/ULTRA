
#include "mask.h"

void OutputMaskedFASTA(const std::string &in_path,
                       FILE* out_file,
                       std::unordered_map<std::string, std::vector<mregion> *> masks,
                       bool n_mask)
{
  std::ifstream in_file(in_path);
  if (!in_file.is_open()) {
    printf("Unable to open %s\n", in_path.c_str());
    return;
  }

  std::string line;
  std::string seq_name = "";
  unsigned long long seq_pos = 0;
  unsigned long long mask_i = 0;
  std::vector<mregion> *cmask = nullptr;

  while (std::getline(in_file, line)) {

    if (line[0] == '>') {
      seq_name = line.substr(1, std::string::npos);
      seq_pos = 0;
      fprintf(out_file, "%s", line.c_str());

      if (masks.find(seq_name) != masks.end()) {
        if (cmask != nullptr)
          delete cmask;
        cmask = CleanedMasks(masks[seq_name]);
        if (cmask->size() == 0) {
          cmask = nullptr;
        }
      }
    }

    else {
      // If there are no masks for this sequence, just skip over it
      if (cmask == nullptr) {
        fprintf(out_file, "%s", line.c_str());
        continue;
      }

      // There are masks for this sequence, check to see if we need to mask
      else {
        for (int i = 0; i < line.length(); ++i) {
          if (line[i] == 'a' || line[i] == 'A' ||
              line[i] == 'c' || line[i] == 'C' ||
              line[i] == 'g' || line[i] == 'G' ||
              line[i] == 't' || line[i] == 'T') {

            if (cmask) {
              if (cmask->at(mask_i).start <= seq_pos) {
                // We need to mask
                if (cmask->at(mask_i).end >= seq_pos) {
                  // We should mask
                  if (n_mask) {
                    line[i] = 'n';
                  }

                  else {
                    line[i] = std::tolower(line[i]);
                  }
                }
                // We are done with cmask_i
                else {
                  mask_i += 1;
                  if (mask_i >= cmask->size()) {
                    cmask = nullptr;
                  }
                }
              }
            }
            seq_pos += 1;
          }

          else if (line[i] == 'n' || line[i] == 'N') {
            seq_pos += 1;
          }
        }
        fprintf(out_file, "%s", line.c_str());
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
