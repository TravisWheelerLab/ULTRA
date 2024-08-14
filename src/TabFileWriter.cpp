//
// Created by Daniel Olson on 08/07/24.
//

#include "TabFileWriter.hpp"
#include "repeat.hpp"
#include "ultra.hpp"
#include <algorithm>
#include <iostream>

void TabFileWriter::InitializeWriter(Ultra *ultra, FILE *out_f) {
  owner = ultra;
  out = out_f;
  fprintf(out, "SeqID");
  fprintf(out, "\tStart");
  fprintf(out, "\tEnd");
  fprintf(out, "\tPeriod");
  fprintf(out, "\tScore");
  if (owner->settings->pval)
    fprintf(out, ",PValue");
  if (owner->settings->max_consensus_period >= 0)
    fprintf(out, "\tConsensus");
  if (owner->settings->max_split >= 0) {
    fprintf(out, "\t#Subrepeats");
    fprintf(out, "\tSubrepeatStarts");
    if (owner->settings->max_consensus_period >= 0)
      fprintf(out, "\tSubrepeatConsensi");
  }

  if (owner->settings->show_seq) {
    fprintf(out, "\tSequence");
  }

  fprintf(out, "\n");

}

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
  fprintf(out, "%s\t%lu\t%lu\t%i\t%f", name.c_str(), repeat->sequenceStart,
          repeat->sequenceStart + repeat->repeatLength, repeat->repeatPeriod, repeat->regionScore);
  if (owner->settings->pval) {
    fprintf(out, ",%g", owner->PvalForScore(repeat->regionScore));
  }

  // We need to decide what to do with the overall sequence
  std::string rep_con = ".";
  if (owner->settings->max_consensus_period != 0) {
    if (owner->settings->max_consensus_period >= repeat->repeatPeriod &&
        !repeat->string_consensus.empty())
      rep_con = repeat->string_consensus;

    fprintf(out, "\t%s", rep_con.c_str());
  }

  if (owner->settings->max_split >= 0) {
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
          }

          sizes += std::to_string(split_i - cstart);
          starts.push_back(',');
          starts += std::to_string(split_i);

          ++numberOfValidSplits;
          cstart = split_i;
        }
      }

      if (owner->settings->max_consensus_period >= 0) {
        for (int i = 0; i < repeat->consensi->size(); ++i) {
          std::string con = ".";
          if (owner->settings->max_consensus_period >= repeat->repeatPeriod) {
            if (repeat->consensi != nullptr && repeat->consensi->size() > i) {
              if (i > 0 && repeat->consensi->at(i) == repeat->consensi->at(i - 1))
                continue;
              con = repeat->consensi->at(i);
            }
          }
          if (i > 0)
            consensi.push_back(',');
          consensi += con;
        }
      }

      if (numberOfValidSplits > 0)
        sizes.push_back(',');
      sizes += std::to_string(repeat->repeatLength - cstart);

      fprintf(out, "\t%i\t%s", numberOfValidSplits + 1,
              starts.c_str());
      if (owner->settings->max_consensus_period >= 0) {
        fprintf(out, "\t%s", consensi.c_str());
      }

    }

    else {
      fprintf(out, "\t1\t0");
      if (owner->settings->max_consensus_period >= 0) {
        fprintf(out, "\t%s", rep_con.c_str());
      }
    }
  }

  if (owner->settings->show_seq) {
    fprintf(out, "\t%s", repeat->sequence.c_str());
  }
  fprintf(out, "\n");
}

void TabFileWriter::EndWriter() {}