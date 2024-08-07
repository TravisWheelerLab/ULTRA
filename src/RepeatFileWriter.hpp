//
// Created by Daniel Olson
//

#ifndef ULTRA_REPEATFILEWRITER_HPP
#define ULTRA_REPEATFILEWRITER_HPP
#include <stdio.h>
class Ultra;
class RepeatRegion;

class RepeatFileWriter {
public:
  virtual void InitializeWriter(Ultra *ultra, FILE *out) {}
  virtual void WriteRepeat(RepeatRegion *repeat) {}
  virtual void EndWriter() {}
};

#endif // ULTRA_REPEATFILEWRITER_HPP
