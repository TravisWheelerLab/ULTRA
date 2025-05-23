//
// Created by Daniel Olson on 3/11/22.
//

#include "JSONFileWriter.hpp"
#include "repeat.hpp"
#include "ultra.hpp"

void JSONFileWriter::OutputJSONKeyValue(std::string key, std::string value,
                                        bool quotes) {
  if (!quotes)
    fprintf(out, ",\n\"%s\": %s", key.c_str(), value.c_str());
  else
    fprintf(out, ",\n\"%s\": \"%s\"", key.c_str(), value.c_str());
}

void JSONFileWriter::InitializeWriter(Ultra *ultra, FILE *out_f) {
  owner = ultra;
  out = out_f;
  fprintf(out, "{\"Repeats\": [\n");
}

std::string JSONFileWriter::StringForSubRepeat(RepeatRegion *r, int split_index,
                                               int start_pos) {
  std::string repeatString = "{";

  int start = start_pos;
  int end = r->repeatLength;

  int consensusPosition = r->splits->size();

  if (split_index >= 0) {
    end = r->splits->at(split_index);
    consensusPosition = split_index;
  }

  float frac_rep = (float)r->repeatLength / (float)(end - start_pos);
  float subScore = r->regionScore / frac_rep;

  repeatString += "\"Start\": ";
  repeatString += std::to_string(start);
  repeatString += ",\n\"Consensus\": \"";
  // Converting to a c string is important.
  // ...Dunno why, but here we are
  repeatString += r->consensi->at(consensusPosition).c_str();

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
      subRepeats += StringForSubRepeat(r, s, start);
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
    fprintf(out, ",\n\n");
  }

  ++this->repeatsOutput;

  fprintf(out, "{\"SequenceName\": \"%s\"", repeat->sequenceName.c_str());

  this->OutputJSONKeyValue("Start", std::to_string(repeat->sequenceStart));
  this->OutputJSONKeyValue("Length", std::to_string(repeat->repeatLength));
  this->OutputJSONKeyValue("Period", std::to_string(repeat->repeatPeriod));
  this->OutputJSONKeyValue("Score", std::to_string(repeat->regionScore));
  if (owner->settings->pval) {
    double pval = owner->PvalForScore(repeat->regionScore);
    this->OutputJSONKeyValue("PVal", std::to_string(pval));
  }

  if (owner->settings->show_counts) {
    auto copies =
        (repeat->repeatLength - repeat->insertions + repeat->deletions) /
        repeat->repeatPeriod;
    this->OutputJSONKeyValue("Copies", std::to_string(copies));
    this->OutputJSONKeyValue("Substitutions",
                             std::to_string(repeat->mismatches));
    this->OutputJSONKeyValue("Insertions", std::to_string(repeat->insertions));
    this->OutputJSONKeyValue("Deletions", std::to_string(repeat->deletions));
  }

  this->OutputJSONKeyValue("Consensus", repeat->string_consensus, true);

  if (owner->settings->show_seq) {
    this->OutputJSONKeyValue("Sequence", repeat->sequence, true);
  }

  if (owner->settings->show_trace) {
    this->OutputJSONKeyValue("Traceback", repeat->traceback, true);
  }

  if (owner->settings->show_logo_nums && repeat->logoNumbers != nullptr) {
    std::string logoNumbers = "\"" + std::to_string(repeat->logoNumbers[0]);
    for (int i = 0; i < repeat->repeatLength; ++i) {
      logoNumbers.append("," + std::to_string(repeat->logoNumbers[i]));
    }

    logoNumbers += "\"";

    this->OutputJSONKeyValue("LogoNumbers", logoNumbers);
  }

  if (owner->settings->show_deltas && repeat->scores != nullptr) {
    std::string positionScoreDeltas = "[";

    for (int i = 0; i < repeat->repeatLength; ++i) {
      if (i > 0)
        positionScoreDeltas += ",";
      positionScoreDeltas += std::to_string(repeat->scores[i]);
    }

    positionScoreDeltas.push_back(']');
    this->OutputJSONKeyValue("PositionScoreDeltas", positionScoreDeltas);
  }

  int numberOfValidSplits = 0;
  if (repeat->splits != nullptr) {
    for (int i = 0; i < repeat->splits->size(); ++i) {
      if (repeat->splits->at(i) > 0)
        ++numberOfValidSplits;
    }
  }

  if (numberOfValidSplits > 0) {
    std::string subRepeats = SubRepeatsString(repeat);
    this->OutputJSONKeyValue("SubRepeats", subRepeats);
  }
  fprintf(out, "}");
}

void JSONFileWriter::EndWriter() { fprintf(out, "]\n}"); }