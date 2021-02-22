//
//  JSONRepeat.cpp
//  ultrax
//
//  Created by Daniel Olson on 9/25/19.
//  Copyright © 2019 Daniel Olson. All rights reserved.
//

#include "JSONRepeat.hpp"

bool JSONRepeat::InterpretRepeat(json11::Json repeat) {
  if (!repeat.is_object())
    return false;

  if (repeat.object_items().find("PassID") == repeat.object_items().end())
    return false;
  passID = repeat["PassID"].int_value();

  if (repeat.object_items().find("Start") == repeat.object_items().end())
    return false;
  start = repeat["Start"].number_value();

  if (repeat.object_items().find("Length") == repeat.object_items().end())
    return false;
  length = repeat["Length"].number_value();

  if (repeat.object_items().find("Period") == repeat.object_items().end())
    return false;
  period = repeat["Period"].number_value();

  if (repeat.object_items().find("Score") == repeat.object_items().end())
    return false;
  score = repeat["Score"].number_value();

  substitutions = repeat["Substitutions"].int_value();
  insertions = repeat["Insertions"].number_value();
  deletions = repeat["Deletions"].number_value();
  consensus = repeat["Consensus"].string_value();
  sequence = repeat["Sequence"].string_value();
  sequenceName = repeat["SequenceName"].string_value();
  traceback = repeat["Traceback"].string_value();
  overlap_code = repeat["OC"].int_value();

  return true;
}

void OutputJSONKey(FILE *out, std::string key) {
  fprintf(out, "\"%s\": ", key.c_str());
}

void JSONRepeat::OutputRepeat(FILE *out, float threshold, bool first) {

  if (score > threshold) {
    if (!first)
      fprintf(out, ",\n");
    fprintf(out, "{");

    OutputJSONKey(out, "PassID");
    fprintf(out, "%i", passID);

    fprintf(out, ",\n");
    OutputJSONKey(out, "SequenceName");
    fprintf(out, "\"%s\"", sequenceName.c_str());

    fprintf(out, ",\n");
    OutputJSONKey(out, "Start");
    fprintf(out, "%lu", start);

    fprintf(out, ",\n");
    OutputJSONKey(out, "Length");
    fprintf(out, "%lu", length);

    fprintf(out, ",\n");
    OutputJSONKey(out, "Period");
    fprintf(out, "%i", period);

    fprintf(out, ",\n");
    OutputJSONKey(out, "Score");
    fprintf(out, "%f", score);

    fprintf(out, ",\n");
    OutputJSONKey(out, "Substitutions");
    fprintf(out, "%i", substitutions);

    fprintf(out, ",\n");
    OutputJSONKey(out, "Insertions");
    fprintf(out, "%i", insertions);

    fprintf(out, ",\n");
    OutputJSONKey(out, "Deletions");
    fprintf(out, "%i", deletions);

    fprintf(out, ",\n");
    OutputJSONKey(out, "Consensus");
    fprintf(out, "\"%s\"", consensus.c_str());

    if (sequence.length() > 0) {
      fprintf(out, ",\n");
      OutputJSONKey(out, "Sequence");
      fprintf(out, "\"%s\"", sequence.c_str());
    }

    if (traceback.length() > 0) {
      fprintf(out, ",\n");
      OutputJSONKey(out, "Traceback");
      fprintf(out, "\"%s\"", traceback.c_str());
    }

    if (subrepeats.size() > 0) {
      fprintf(out, ",\n");
      OutputJSONKey(out, "Subrepeats");
      fprintf(out, "[");
    }

    // correct the subrepeats array before OutRepeating them
    for (int i = 0; i < subrepeats.size(); ++i) {
      subrepeats[i]->OutputRepeat(out, threshold, i == 0);
    }

    if (subrepeats.size() > 0) {
      fprintf(out, "]\n");
    }

    fprintf(out, "}");
  }
}

void JSONRepeat::FixSubrepeatOverlap() {
  if (subrepeats.size() < 2)
    return;

  JSONRepeat *n = subrepeats.back();
  subrepeats.pop_back();
  JSONRepeat *c = subrepeats.back();
  subrepeats.pop_back();

  if (c->start + c->length > n->start) {
    subrepeats.push_back(c);
    subrepeats.push_back(n);
  }
}
