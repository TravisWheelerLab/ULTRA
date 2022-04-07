//
// Created by Daniel Olson on 3/10/22.
//

#include "BEDFileWriter.hpp"
#include "repeat.hpp"
#include "ultra.hpp"
#include <algorithm>

void BEDFileWriter::InitializeWriter(Ultra *ultra) {
  owner = ultra;
}

void BEDFileWriter::WriteRepeat(RepeatRegion *repeat) {

  // We need a better behavior here - check with Travis
  std::string name = repeat->sequenceName;
  int j = 0;
  for (int i = 0; i < name.size(); ++i) {
    if ((name[i] >= 'a' && name[i] <= 'z') ||
        (name[i] >= 'A' && name[i] <= 'Z') ||
        (name[i] >= '0' && name[i] <= '9') ||
        name[i] == '_') {
        name[j] = name[i];
        ++j;
    }

    else {
      name[j] = '\0';
    }
  }
  name[j] = '\0';

  fprintf(owner->out,"%s", name.c_str(),
                            repeat->sequenceStart,
    repeat->sequenceStart + repeat->repeatLength);
  fprintf(owner->out, " incomplete %f", repeat->regionScore);

  // In the future only print this section if we are doing repeat
  //    splititng
  fprintf(owner->out, " . %lu %lu 0", repeat->sequenceStart,
              repeat->sequenceStart + repeat->repeatLength);
  fprintf(owner->out, "\n");

}

void BEDFileWriter::EndWriter() {

}