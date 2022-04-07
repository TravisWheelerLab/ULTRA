//
// Created by Daniel Olson on 3/10/22.
//

#include "BEDFileWriter.hpp"
#include "repeat.hpp"
#include "ultra.hpp"

void BEDFileWriter::InitializeWriter(Ultra *ultra) {
  owner = ultra;
}

void BEDFileWriter::WriteRepeat(RepeatRegion *repeat) {
 // printf("test\n");
  //printf("Test it %llu\n", (unsigned long long)repeat);
  return;

  printf("%llu\n", (unsigned long long)repeat);
  printf("A %llu\n", (unsigned long long)&repeat->sequenceName);
  printf("%s\n", repeat->sequenceName.c_str());
  printf("%llu\n", (unsigned long long)owner->out);
  printf("We are reaching BED out\n"); // remove this eventually
  printf("A %llu\n", (unsigned long long)&repeat->sequenceName);
  printf(/*owner->out, */"test: \"%s\"", repeat->sequenceName.c_str());

      //                              repeat->sequenceStart,
       //     repeat->sequenceStart + repeat->repeatLength);
  fprintf(owner->out, " incomplete %f", repeat->regionScore);

  // In the future only print this section if we are doing repeat
  //    splititng
  fprintf(owner->out, " . %lu %lu 0", repeat->sequenceStart,
              repeat->sequenceStart + repeat->repeatLength);

}

void BEDFileWriter::EndWriter() {

}