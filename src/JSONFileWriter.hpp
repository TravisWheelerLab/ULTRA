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

  void OutputJSONKeyValue(std::string key, std::string value, bool quotes=false);

public:
  Ultra *owner;
  void InitializeWriter(Ultra *ultra);
  void WriteRepeat(RepeatRegion *repeat);
  void EndWriter();
  std::string StringForSubRepeat(RepeatRegion *r, int split_index, int start_pos);
  std::string SubRepeatsString(RepeatRegion *r);
};

#endif // ULTRA_JSONFILEWRITER_HPP
