//
// Created by Daniel Olson on 4/14/23.
//

#include "cli.hpp"

Settings_t::Settings_t() {

  // *************
  // Input options
  // *************
  app.add_option("input_file",
                 this->in_file,
                 "DNA FASTA input file")
      ->required(true)
      ->group("Input");

  app.add_flag("-r,--readall",
               this->read_all,
               "Read entire input file into memory"
               "(disables streaming input)")->group("Input");


  // *************
  // Output options
  // *************

  app.add_flag("--pval",
               this->pval,
               "Use p-values instead of scores in BED output")
      ->group("Output");

  app.add_flag("-j,--json",
               this->json,
               "Use JSON outuput format instead of BED")
      ->group("Output");

  app.add_flag("--hideseq",
               this->hide_seq,
               "Hide sequence descriptor in JSON output")
      ->group("Output");

  app.add_flag("--showdelta",
               this->show_deltas,
               "Show positional score deltas in JSON output")
      ->group("Output");

  app.add_flag("--showtrace",
               this->show_trace,
               "Show Viterbi trace in JSON output")
      ->group("Output");

  app.add_flag("--showwid",
               this->show_wid,
               "Show sequence window IDs in JSON output")
      ->group("Output");

  app.add_flag("--showlogo",
               this->show_logo_nums,
               "Show logo numbers in JSON output")
      ->group("Output");

  app.add_flag("--suppress",
               this->suppress_out,
               "Do not output BED or JSON annotation")
      ->group("Output");

  // *************
  // Mask options
  // *************

  app.add_option("-o,--out",
               this->out_file,
               "Output file path")
      ->group("Masking");

  app.add_option("--mask",
                 this->mask_file,
                 "File path to save a masked FASTA")
      ->group("Masking");

  app.add_flag("--nmask",
               mask_with_n,
               "Use n masking instead of case masking")
      ->group("Masking");

  // *************
  // System options
  // *************
  app.add_option("-t,--threads",
                 this->threads,
                 "Number of threads to use")
      ->default_str("1")
      ->group("System");

  app.add_option("--winsize",
                 this->window_size,
                 "Manually set sequence window size")
  ->group("System");

  app.add_option("--overlap",
                 this->overlap,
                 "Manually set sequence window overlap size")
      ->group("System");

  app.add_option("--windows",
                 this->windows,
                 "Number of sequence windows to store at once")
      ->default_val(this->windows)
      ->group("System");

  app.add_flag("--mem",
                 this->show_memory,
                 "Display memory requirements for current settings")
      ->group("System");

  // *************
  // Filter options
  // *************
  app.add_option("-s, --minscore",
                 this->min_score,
                 "Minimum reportable repeat score")
      ->default_val(this->min_score)
      ->group("Filter");

  app.add_option("--minlen",
                 this->min_length,
                 "Minimum reportable repeat length")
      ->default_val(this->min_length)
      ->group("Filter");

  app.add_option("--minunit",
                 this->min_units,
                 "Minimum reportable number of repeat units")
      ->default_val(this->min_units)
      ->group("Filter");

  // *************
  // Model options
  // *************
  app.add_option("-p, --period",
                 this->max_period,
                 "Maximum detectable repeat period")
      ->default_val(this->max_period)
      ->group("Model");

  app.add_option("-i, --inserts",
                 this->max_insert,
                 "Maximum number of insertion states")
      ->default_val(this->max_insert)
      ->group("Model");

  app.add_option("-d, --deletes",
                 this->max_delete,
                 "Maximum number of deletion states")
      ->default_val(this->max_delete)
      ->group("Model");

  // *************
  // Probability options
  // *************
  app.add_option("--at",
                 this->at,
                 "Expected AT content")
      ->default_val(this->at)
      ->group("Probabilities");

  app.add_option("--acgt",
                 this->acgt,
                 "Expected frequency of A C G T")
      ->default_val(this->acgt)
      ->expected(4)
      ->group("Probabilities");

  app.add_option("-m,--match",
                 this->match_probability,
                 "Expected pattern conservation in repeats")
      ->default_val(this->match_probability)
      ->group("Probabilities");

  app.add_option("--decay",
                 this->period_decay,
                 "Decay penalty applied to repetitive states")
      ->default_val(this->period_decay)
      ->group("Probabilities");

  app.add_option("--nr",
                 this->transition_nr,
                 "Probability of a repeat starting")
      ->default_val(this->transition_nr)
      ->group("Probabilities");

  app.add_option("--rn",
                 this->transition_rn,
                 "Probability of a repeat stopping")
      ->default_val(this->transition_rn)
      ->group("Probabilities");

  app.add_option("--ri",
                 this->transition_ri,
                 "Probability of an insertion occurring")
      ->default_val(this->transition_ri)
      ->group("Probabilities");

  app.add_option("--rd",
                 this->transition_rd,
                 "Probability of a deletion occurring")
      ->default_val(this->transition_rd)
      ->group("Probabilities");

  app.add_option("--ii",
                 this->transition_ii,
                 "Probability of consecutive insertions occurring")
      ->default_val(this->transition_ii)
      ->group("Probabilities");

  app.add_option("--dd",
                 this->transition_dd,
                 "Probability of consecutive deletions occurring")
      ->default_val(this->transition_dd)
      ->group("Probabilities");

  // *************
  // Split options
  // *************

  app.add_flag("--nosplit",
               this->no_split,
               "Do not perform repeat splitting")
  ->group("Splitting");

  app.add_option("--maxsplit",
                 this->max_split,
                 "The maximum repeat period to perform repeat splitting")
      ->default_val(this->max_split)
      ->group("Splitting");

  app.add_option("--maxsplit",
                 this->split_threshold,
                 "Split threshold value")
      ->default_val(this->split_threshold)
      ->group("Splitting");

  app.add_option("--splitdepth",
                 this->split_depth,
                 "Number of repeat units to use in repeat splitting")
      ->default_val(this->split_depth)
      ->group("Splitting");

  app.add_option("--minsplitwin",
                 this->min_split_window,
                 "Minimum repeat split window size")
      ->default_val(this->min_split_window)
      ->group("Splitting");


}