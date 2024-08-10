//
// Created by Daniel Olson on 3/10/22.
//

#include "BEDFileWriter.hpp"
#include "repeat.hpp"
#include "ultra.hpp"
#include <algorithm>
#include <iostream>
void BEDFileWriter::InitializeWriter(Ultra *ultra, FILE *out_file) { owner = ultra; out=out_file; }

void BEDFileWriter::WriteRepeat(RepeatRegion *repeat) {

  // We need a better behavior here - check with Travis
  std::string name = repeat->sequenceName;
  for (int i = 0; i < name.size(); ++i) {
    if ((name[i] >= 'a' && name[i] <= 'z') ||
        (name[i] >= 'A' && name[i] <= 'Z') ||
        (name[i] >= '0' && name[i] <= '9') || name[i] == '-' ||
        name[i] == '_' || name[i] == '.' || name[i] == ':' || name[i] == '*' ||
        name[i] == '#') {
      continue;
    }

    else {
      name = name.substr(0, i);
      break;
    }
  }

  // Columns 1 (name) 2(start) 3 (end)
  fprintf(out, "%s\t%lu\t%lu", name.c_str(), repeat->sequenceStart,
          repeat->sequenceStart + repeat->repeatLength);

  // We need to decide what to do with the overall sequence

  std::string rep_con = std::to_string(repeat->repeatPeriod);
  if (owner->settings->max_consensus_period >= repeat->repeatPeriod && !repeat->string_consensus.empty())
    rep_con = repeat->string_consensus;

  fprintf(out, "\t%s\n", rep_con.c_str());
}

void BEDFileWriter::EndWriter() {}