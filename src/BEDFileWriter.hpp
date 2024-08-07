//
// Created by Daniel Olson on 3/10/22.
//

#ifndef ULTRA_BEDFILEWRITER_HPP
#define ULTRA_BEDFILEWRITER_HPP
#include "RepeatFileWriter.hpp"

class RepeatRegion;
class Ultra;
class BEDFileWriter : virtual public RepeatFileWriter {
public:
  Ultra *owner;
  FILE *out;
  bool outputMultilineSplitRepeats = false;

  void InitializeWriter(Ultra *ultra, FILE *out);
  void WriteRepeat(RepeatRegion *repeat);
  void EndWriter();
};

#endif // ULTRA_BEDFILEWRITER_HPP
