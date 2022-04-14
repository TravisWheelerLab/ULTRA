//
// Created by Daniel Olson on 3/11/22.
//

#include "JSONFileWriter.hpp"
#include "repeat.hpp"
#include "ultra.hpp"

void JSONFileWriter::OutputJSONKeyValue(std::string key, std::string value) {
  fprintf(owner->out, ",\n\"%s\": %s", key.c_str(), value.c_str());
}

void JSONFileWriter::InitializeWriter(Ultra *ultra) {
  owner = ultra;
  fprintf(owner->out, "{\"Repeats\": [");
}

void JSONFileWriter::WriteRepeat(RepeatRegion *repeat) {

  if (this->repeatsOutput > 0) {
    fprintf(owner->out, ",\n\n");
  }

  ++this->repeatsOutput;

  fprintf(owner->out, "{\"SequenceName\": \"%s\"",
          repeat->sequenceName.c_str());

  this->OutputJSONKeyValue("Start", std::to_string(repeat->sequenceStart));
  this->OutputJSONKeyValue("Length", std::to_string(repeat->repeatLength));
  this->OutputJSONKeyValue("Period", std::to_string(repeat->repeatPeriod));
  this->OutputJSONKeyValue("Score", std::to_string(repeat->regionScore));
  if (owner->settings->v_calculateLogPVal) {
    this->OutputJSONKeyValue("Log2PVal", std::to_string(repeat->logPVal));
  }

  this->OutputJSONKeyValue("Substitutions", std::to_string(repeat->mismatches));
  this->OutputJSONKeyValue("Insertions", std::to_string(repeat->insertions));
  this->OutputJSONKeyValue("Deletions", std::to_string(repeat->deletions));
  this->OutputJSONKeyValue("Consensus", repeat->GetConsensus());

  if (owner->outputReadID) {
    this->OutputJSONKeyValue("ReadID", std::to_string(repeat->readID));
  }

  if (owner->outputRepeatSequence) {
    this->OutputJSONKeyValue("Sequence", repeat->sequence);
  }

  if (owner->settings->v_showTraceback) {
    this->OutputJSONKeyValue("Traceback", "\"" + repeat->traceback + "\"");
  }

  if (owner->settings->v_showLogoNumbers) {
    std::string logoNumbers = "\"" + std::to_string(repeat->logoNumbers[0]);
    for (int i = 0; i < repeat->repeatLength; ++i) {
      logoNumbers.append("," + std::to_string(repeat->logoNumbers[i]));
    }

    logoNumbers += "\"";

    this->OutputJSONKeyValue("LogoNumbers", logoNumbers);
  }

  // TODO
  // POSITION SCORE DELTAS HERE

  // TODO
  // split repeats here

  if (owner->settings->v_splitRepeats) {
  }
}

void JSONFileWriter::EndWriter() { fprintf(owner->out, "]\n}"); }