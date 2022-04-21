//
// Created by Daniel Olson on 3/10/22.
//

#include "BEDFileWriter.hpp"
#include "repeat.hpp"
#include "ultra.hpp"
#include <algorithm>

void BEDFileWriter::InitializeWriter(Ultra *ultra) { owner = ultra;

}

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

  fprintf(owner->out, "%s %lu %lu", name.c_str(), repeat->sequenceStart,
          repeat->sequenceStart + repeat->repeatLength);
  // We need to decide what to do with the overall sequence
  fprintf(owner->out, " incomplete %f", repeat->regionScore);

  fprintf(owner->out, " . %lu %lu 0", repeat->sequenceStart,
          repeat->sequenceStart + repeat->repeatLength);
  fprintf(owner->out, "\n");

  //TODO
  // Complete repeat split report
  /*
  if (repeat->splits->size() > 0) {
    fprintf(owner->out, "%lu", repeat->splits->size());
    int prevSplit = 0;

    for (int i = 0; i < repeat->splits->size()) {

    }
  }*/
}

void BEDFileWriter::EndWriter() {

}