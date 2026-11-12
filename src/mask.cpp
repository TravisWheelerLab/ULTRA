
#include "mask.hpp"
#include <cctype>

void OutputMaskedFASTA(
    const std::string &in_path, FILE *out_file,
    std::unordered_map<unsigned long long, std::vector<mregion> *> masks,
    bool n_mask) {
  std::ifstream in_file(in_path);
  if (!in_file.is_open()) {
    fprintf(stderr, "Unable to open %s\n", in_path.c_str());
    return;
  }

  std::string line;
  std::string seq_name = "";
  long long seq_id = 0;
  unsigned long long seq_pos = 0;
  unsigned long long mask_i = 0;
  std::vector<mregion> *cmask = nullptr;

  while (std::getline(in_file, line)) {
    if (line[0] == '>') {
      seq_name = line.substr(1, std::string::npos);
      seq_pos = 0;
      seq_id += 1;
      mask_i = 0;

      // Each record is masked with its own regions only
      if (cmask != nullptr) {
        delete cmask;
        cmask = nullptr;
      }
      if (masks.find(seq_id) != masks.end()) {
        cmask = CleanedMasks(masks[seq_id]);
      }
    }

    // Comment lines are not sequence; nothing to do once the masks run out
    else if (cmask != nullptr && line[0] != ';' && mask_i < cmask->size()) {
      for (int i = 0; i < line.length(); ++i) {
        // Count positions the way SequenceWindow::ReadLine does
        if (!isalnum((unsigned char)line[i]))
          continue;

        // Regions are sorted and disjoint
        while (mask_i < cmask->size() && cmask->at(mask_i).end <= seq_pos)
          mask_i += 1;

        if (mask_i < cmask->size() && cmask->at(mask_i).start <= seq_pos) {
          char u = std::toupper((unsigned char)line[i]);
          if (u == 'A' || u == 'C' || u == 'G' || u == 'T') {
            if (n_mask)
              line[i] = 'n';
            else
              line[i] = std::tolower((unsigned char)line[i]);
          }
        }
        seq_pos += 1;
      }
    }

    fprintf(out_file, "%s\n", line.c_str());
  }
  if (cmask != nullptr)
    delete cmask;
  in_file.close();
}

bool mregion::operator()(mregion m1, mregion m2) {
  if (m1.start != m2.start)
    return (m1.start < m2.start);

  else
    return m1.end < m2.end;
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
