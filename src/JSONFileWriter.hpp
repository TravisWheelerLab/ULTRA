//
// Created by Daniel Olson on 3/11/22.
//

#ifndef ULTRA_JSONFILEWRITER_HPP
#define ULTRA_JSONFILEWRITER_HPP
#include "RepeatFileWriter.hpp"

class RepeatRegion;
class Ultra;
class JSONFileWriter : virtual public RepeatFileWriter
{
public:
  Ultra *owner;
  void InitializeWriter(Ultra *ultra);
  void WriteRepeat(RepeatRegion *repeat);
  void EndWriter();
}

#endif // ULTRA_JSONFILEWRITER_HPP
