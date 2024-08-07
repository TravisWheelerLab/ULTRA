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

  // Columns 1 (name) 2(start) 3 (end)
  fprintf(owner->out, "%s\t%lu\t%lu", name.c_str(), repeat->sequenceStart,
          repeat->sequenceStart + repeat->repeatLength);

  // We need to decide what to do with the overall sequence

  std::string rep_con = ".";
  if (!repeat->string_consensus.empty())
    rep_con = repeat->string_consensus;
  // Columns 4 (name) 5 (score) 6 (strand=.) 7 thickstart 8 thickend 9 rgb
  fprintf(owner->out, "\t%s\t%f", rep_con.c_str(), repeat->regionScore);
  fprintf(owner->out, "\n");
}

void BEDFileWriter::EndWriter() {}