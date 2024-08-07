//
// Created by Daniel Olson on 3/10/22.
//

#include "TabFileWriter.hpp"
#include "repeat.hpp"
#include "ultra.hpp"
#include <algorithm>
#include <iostream>
void TabFileWriter::InitializeWriter(Ultra *ultra) { owner = ultra; }

void TabFileWriter::WriteRepeat(RepeatRegion *repeat) {

  // We need a better behavior here - check with Travis
  std::string name = repeat->sequenceName;
  for (int i = 0; i < name.size(); ++i) {
    if ((name[i] >= 'a' && name[i] <= 'z') ||
        (name[i] >= 'A' && name[i] <= 'Z') ||
        (name[i] >= '0' && name[i] <= '9') ||
        name[i] == '-' || name[i] == '_' || name[i] == '.' ||
        name[i] == ':' || name[i] == '*' || name[i] == '#') {
      continue;
    }

    else {
      name = name.substr(0, i);
      break;
    }

   /* else {
      name[i] = '_';
    }*/
  }

  // Columns 1 (name) 2 (start) 3 (end) 4 (score)
  fprintf(owner->out, "%s\t%lu\t%lu\t%i\t%f", name.c_str(), repeat->sequenceStart,
          repeat->sequenceStart + repeat->repeatLength, repeat->repeatPeriod, repeat->regionScore);
  if (owner->settings->pval) {
    fprintf(owner->out, ",%f", owner->PvalForScore(repeat->regionScore));
  }

  // We need to decide what to do with the overall sequence
  std::string rep_con = ".";
  if (owner->settings->max_consensus_period != 0) {
    if (owner->settings->max_consensus_period >= repeat->repeatPeriod &&
        !repeat->string_consensus.empty())
      rep_con = repeat->string_consensus;

    printf("\t%s", rep_con.c_str());
  }

  if (owner->settings->max_split > 0) {
    std::string sizes = "";
    std::string starts = "0";
    std::string consensi = "";
    int numberOfValidSplits = 0;

    int cstart = 0;
    if (repeat->splits != nullptr && !repeat->splits->empty()) {
      for (int i = 0; i < repeat->splits->size(); ++i) {
        int split_i = repeat->splits->at(i);
        if (split_i > 0) {
          if (numberOfValidSplits > 0) {
            sizes.push_back(',');
            consensi.push_back(',');
          }

          sizes += std::to_string(split_i - cstart);
          starts.push_back(',');
          starts += std::to_string(split_i);

          if (owner->settings->max_consensus_period != 0) {
            std::string con = ".";
            if (owner->settings->max_consensus_period >= repeat->repeatPeriod) {
              if (repeat->consensi->size() > i) {
                con = repeat->consensi->at(i);
              }
            }
            consensi += con;
          }

          ++numberOfValidSplits;
          cstart = split_i;
        }
      }

      if (numberOfValidSplits > 0)
        sizes.push_back(',');
      sizes += std::to_string(repeat->repeatLength - cstart);

      fprintf(owner->out, "%i\t%s\t%s", numberOfValidSplits + 1, sizes.c_str(),
              starts.c_str());
      if (owner->settings->max_consensus_period != 0) {
        fprintf(owner->out, "\t%s", consensi.c_str());
      }

    }

    else {
      fprintf(owner->out, "1\t%lu\t0", repeat->repeatLength);
      if (owner->settings->max_consensus_period != 0) {
        fprintf(owner->out, "\t%s", rep_con.c_str());
      }
    }
  }

  fprintf(owner->out, "\n");
}

void TabFileWriter::EndWriter() {}