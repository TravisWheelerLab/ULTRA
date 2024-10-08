//
// Created by Daniel Olson on 08/07/24.
//

#ifndef ULTRA_TABFILEWRITER_HPP
#define ULTRA_TABFILEWRITER_HPP
#include "RepeatFileWriter.hpp"

class RepeatRegion;
class Ultra;
class TabFileWriter : virtual public RepeatFileWriter {
public:
  Ultra *owner;
  FILE *out;

  void InitializeWriter(Ultra *ultra, FILE *out);
  void WriteRepeat(RepeatRegion *repeat);
  void EndWriter();
};

#endif // ULTRA_TABFILEWRITER_HPP
