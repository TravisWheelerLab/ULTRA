//
// Created by Daniel Olson on 3/10/22.
//

#include "TabFileWriter.hpp"
#include "repeat.hpp"
#include "ultra.hpp"
#include <algorithm>
#include <iostream>
void TabFileWriter::InitializeWriter(Ultra *ultra) { owner = ultra; }

void TabFileWriter::WriteRepeat(RepeatRegion *repeat) {

  // We need a better behavior here - check with Travis
  std::string name = repeat->sequenceName;
  for (int i = 0; i < name.size(); ++i) {
    if ((name[i] >= 'a' && name[i] <= 'z') ||
        (name[i] >= 'A' && name[i] <= 'Z') ||
        (name[i] >= '0' && name[i] <= '9') ||
        name[i] == '-' || name[i] == '_' || name[i] == '.' ||
        name[i] == ':' || name[i] == '*' || name[i] == '#') {
      continue;
    }

    else {
      name = name.substr(0, i);
      break;
    }

   /* else {
      name[i] = '_';
    }*/
  }

  // Columns 1 (name) 2(start) 3 (end)
  fprintf(owner->out, "%s\t%lu\t%lu", name.c_str(), repeat->sequenceStart,
          repeat->sequenceStart + repeat->repeatLength);

  // We need to decide what to do with the overall sequence

  std::string rep_con = ".";
  if (!repeat->string_consensus.empty())
    rep_con = repeat->string_consensus;
  // Columns 4 (name) 5 (score) 6 (strand=.) 7 thickstart 8 thickend 9 rgb
  fprintf(owner->out, "\t%s\t%f\t.\t%lu\t%lu\t0,0,0\t", rep_con.c_str(), repeat->regionScore,
          repeat->sequenceStart, repeat->sequenceStart + repeat->repeatLength);

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

void TabFileWriter::EndWriter() {}