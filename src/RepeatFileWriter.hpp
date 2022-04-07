//
// Created by Daniel Olson
//

#ifndef ULTRA_REPEATFILEWRITER_HPP
#define ULTRA_REPEATFILEWRITER_HPP

class Ultra;
class RepeatRegion;

class RepeatFileWriter {
public:

  virtual void InitializeWriter(Ultra *ultra) {}
  virtual void WriteRepeat(RepeatRegion *repeat) {}
  virtual void EndWriter() {}

};


#endif // ULTRA_REPEATFILEWRITER_HPP
