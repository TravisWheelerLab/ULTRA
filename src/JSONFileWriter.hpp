//
// Created by Daniel Olson on 3/11/22.
//

#ifndef ULTRA_JSONFILEWRITER_HPP
#define ULTRA_JSONFILEWRITER_HPP
#include "RepeatFileWriter.hpp"
#include <string>

class RepeatRegion;
class Ultra;
class JSONFileWriter : virtual public RepeatFileWriter {
  int repeatsOutput = 0;

  void OutputJSONKeyValue(std::string key, std::string value);

public:
  Ultra *owner;
  void InitializeWriter(Ultra *ultra);
  void WriteRepeat(RepeatRegion *repeat);
  void EndWriter();
};

#endif // ULTRA_JSONFILEWRITER_HPP
