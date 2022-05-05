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

std::string JSONFileWriter::StringForSubRepeat(RepeatRegion *r, int split_index, int start_pos) {
  std::string repeatString = "{";

  int start = start_pos + r->sequenceStart;
  int end = r->repeatLength;

  int consensusPosition = split_index;

  if (split_index >= 0) {
    end = r->splits->at(split_index);
    consensusPosition = r->splits->size();
  }

  float subScore = r->regionScore / (float)(end - start);

  repeatString += "Start: ";
  repeatString += std::to_string(start);
  repeatString += ",\nEnd: ";
  repeatString += std::to_string(end);
  repeatString += ",\nScore: ";
  repeatString += std::to_string(subScore);
  repeatString += "Consensus: \"";
  repeatString += r->consensi->at(consensusPosition);


  repeatString += "\"}";


  return repeatString;
}

std::string JSONFileWriter::SubRepeatsString(RepeatRegion *r) {
  std::string subRepeats = "[";

  int start = 0;
  for (int s = 0; s < r->splits->size(); ++s) {
    int split_i = r->splits->at(s);
    if (split_i >= 0) {

      if (start > 0)
        subRepeats += ",\n";
      subRepeats += StringForSubRepeat(r, split_i, start);

      start = split_i;
    }
  }

  if (start > 0)
    subRepeats += ",\n";
  subRepeats += StringForSubRepeat(r, -1, start);

  subRepeats.push_back(']');
  return subRepeats;
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


  if (owner->settings->v_showScores) {
    std::string positionScoreDeltas = "[";

    for (int i = 0; i < repeat->repeatLength; ++i) {
      if (i > 0)
        positionScoreDeltas += ",";
      positionScoreDeltas += std::to_string(repeat->scores[i]);
    }

    this->OutputJSONKeyValue("LogoNumbers", positionScoreDeltas);
    positionScoreDeltas.push_back(']');
  }


  int numberOfValidSplits = 0;
  if (repeat->splits != nullptr)
    for (int i = 0; i < repeat->splits->size(); ++i) {
      if (repeat->splits->at(i) > 0)
        ++numberOfValidSplits;
    }

  if (numberOfValidSplits > 0) {
    std::string subRepeats = SubRepeatsString(repeat);
    this->OutputJSONKeyValue("SubRepeats", subRepeats);
  }

}

void JSONFileWriter::EndWriter() { fprintf(owner->out, "]\n}"); }