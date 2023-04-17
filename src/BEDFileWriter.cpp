//
// Created by Daniel Olson on 3/10/22.
//

#include "BEDFileWriter.hpp"
#include "repeat.hpp"
#include "ultra.hpp"
#include <algorithm>
#include <iostream>
void BEDFileWriter::InitializeWriter(Ultra *ultra) { owner = ultra; }

void BEDFileWriter::WriteRepeat(RepeatRegion *repeat) {

  // We need a better behavior here - check with Travis
  std::string name = repeat->sequenceName;
  for (int i = 0; i < name.size(); ++i) {
    if ((name[i] >= 'a' && name[i] <= 'z') ||
        (name[i] >= 'A' && name[i] <= 'Z') ||
        (name[i] >= '0' && name[i] <= '9')) {
      continue;
    }

    else if (name[i] == ' ') {
      name[i] = '\0';
    }

    else {
      name[i] = '_';
    }
  }

  fprintf(owner->out, "%s\t%lu\t%lu", name.c_str(), repeat->sequenceStart,
          repeat->sequenceStart + repeat->repeatLength);
          
  // We need to decide what to do with the overall sequence

  fprintf(owner->out, "\t.\t%lu\t%lu\t0,0,0\t", repeat->sequenceStart,
          repeat->sequenceStart + repeat->repeatLength);

  if (owner->settings->max_split > 0) {

    std::string sizes = "";
    std::string starts = "0";
    int numberOfValidSplits = 0;

    int cstart = 0;
    if (repeat->splits != nullptr && !repeat->splits->empty()) {
      for (int i = 0; i < repeat->splits->size(); ++i) {
        int split_i = repeat->splits->at(i);
        if (split_i > 0) {
          if (numberOfValidSplits > 0)
            sizes.push_back(',');
          sizes += std::to_string(split_i - cstart);
          starts.push_back(',');
          starts += std::to_string(split_i);

          ++numberOfValidSplits;
          cstart = split_i;
        }
      }

      if (numberOfValidSplits > 0)
        sizes.push_back(',');
      sizes += std::to_string(repeat->repeatLength - cstart);

      fprintf(owner->out, "%i\t%s\t%s", numberOfValidSplits + 1, sizes.c_str(),
              starts.c_str());

    }

    else {
      fprintf(owner->out, "1\t%lu\t0", repeat->repeatLength);
    }
  }

  fprintf(owner->out, "\n");
}

void BEDFileWriter::EndWriter() {}