//
// Created by Daniel Olson on 4/14/23.
//

#include "cli.hpp"
#include <cstring>
void Settings::prepare_settings() {

  app.footer("For additional information see README\n"
             "Supported by: NIH NIGMS P20GM103546 and NIH NHGRI U24HG010136\n");

  // *************
  // Input options
  // *************
  app.add_option("input_file", this->in_file, "DNA FASTA input file")
      ->required(false)
      ->group("Input");

  app.add_flag("-r,--readall", this->read_all,
               "Read entire input file into memory"
               " (disables streaming input)")
      ->group("Input");

  // *************
  // Output options
  // *************

  app.add_option("-o,--out", this->out_file, "Output file path")
      ->group("Output");
  app.add_flag("--disable_streaming_out", this->disable_streaming_out,
               "Disables streaming output; no output will be created until all analysis has been completed")
      ->group("Output");

  app.add_flag("--pval", this->pval,
               "Use p-values instead of scores in BED output")
      ->group("Output");

  app.add_flag("-j,--json", this->json,
               "Use JSON outuput format instead of BED")
      ->group("Output");

  app.add_flag("--hideseq", this->hide_seq,
               "Hide sequence descriptor in JSON output")
      ->group("Output");

  app.add_flag("--showdelta", this->show_deltas,
               "Show positional score deltas in JSON output")
      ->group("Output");

  app.add_flag("--showtrace", this->show_trace,
               "Show Viterbi trace in JSON output")
      ->group("Output");

  app.add_flag("--showwid", this->show_wid,
               "Show sequence window IDs in JSON output")
      ->group("Output");

  app.add_flag("--showlogo", this->show_logo_nums,
               "Show logo numbers in JSON output")
      ->group("Output");

  app.add_flag("--hs, --hidesettings", this->hide_settings,
               "Do not output settings")
      ->group("Output");

  app.add_flag("--suppress", this->suppress_out,
               "Do not output BED or JSON annotation")
      ->group("Output");

  // *************
  // Mask options
  // *************

  app.add_option("--mask", this->mask_file, "File path to save a masked FASTA")
      ->group("Masking");

  app.add_flag("--nmask", mask_with_n, "Use n masking instead of case masking")
      ->group("Masking");

  // *************
  // System options
  // *************
  app.add_option("-t,--threads", this->threads, "Number of threads to use")
      ->default_str("1")
      ->group("System");

  app.add_option("--winsize", this->window_size,
                 "Manually set sequence window size")
      ->group("System");

  app.add_option("--overlap", this->overlap,
                 "Manually set sequence window overlap size")
      ->group("System");

  app.add_option("--windows", this->windows,
                 "Number of sequence windows to store at once")
      ->default_val(this->windows)
      ->group("System");

  app.add_flag("--mem", this->show_memory,
               "Display memory requirements for current settings and exit")
      ->group("System");

  // *************
  // Filter options
  // *************
  app.add_option("-s, --score", this->min_score,
                 "Minimum reportable repeat score")
      ->default_str("-100.0")
      ->group("Filter");

  app.add_option("--minlen", this->min_length,
                 "Minimum reportable repeat length")
      ->default_val(this->min_length)
      ->group("Filter");

  app.add_option("--minunit", this->min_units,
                 "Minimum reportable number of repeat units")
      ->default_val(this->min_units)
      ->group("Filter");

  // *************
  // Model options
  // *************
  app.add_option("-p, --period", this->max_period,
                 "Maximum detectable repeat period")
      ->default_val(this->max_period)
      ->group("Model");

  app.add_option("-i, --inserts", this->max_insert,
                 "Maximum number of insertion states")
      ->default_val(this->max_insert)
      ->group("Model");

  app.add_option("-d, --deletes", this->max_delete,
                 "Maximum number of deletion states")
      ->default_val(this->max_delete)
      ->group("Model");

  // *************
  // Probability options
  // *************
  app.add_option("--at", this->at, "Expected AT content")
      ->default_val(this->at)
      ->group("Probabilities");

  app.add_option("--acgt", this->acgt, "Expected frequency of A C G T")
      ->expected(4)
      ->group("Probabilities");

  app.add_option("-m,--match", this->match_probability,
                 "Expected conservation in repeats")
      ->default_val(this->match_probability)
      ->group("Probabilities");

  app.add_option("--decay", this->period_decay,
                 "Decay penalty applied to repetitive states")
      ->default_val(this->period_decay)
      ->group("Probabilities");

  app.add_option("--nr", this->transition_nr,
                 "Probability of a repeat starting")
      ->default_val(this->transition_nr)
      ->group("Probabilities");

  app.add_option("--rn", this->transition_rn,
                 "Probability of a repeat stopping")
      ->default_val(this->transition_rn)
      ->group("Probabilities");

  app.add_option("--ri", this->transition_ri,
                 "Probability of an insertion occurring")
      ->default_val(this->transition_ri)
      ->group("Probabilities");

  app.add_option("--rd", this->transition_rd,
                 "Probability of a deletion occurring")
      ->default_val(this->transition_rd)
      ->group("Probabilities");

  app.add_option("--ii", this->transition_ii,
                 "Probability of consecutive insertions occurring")
      ->default_val(this->transition_ii)
      ->group("Probabilities");

  app.add_option("--dd", this->transition_dd,
                 "Probability of consecutive deletions occurring")
      ->default_val(this->transition_dd)
      ->group("Probabilities");

  // *************
  // Split options
  // *************

  app.add_flag("--nosplit", this->no_split, "Do not perform repeat splitting")
      ->group("Splitting and Naming");

  app.add_option("--maxsplit", this->max_split,
                 "The maximum repeat period to perform repeat splitting")
      ->default_val(this->max_split)
      ->group("Splitting and Naming");

  app.add_option("--splitval", this->split_threshold, "Split threshold value")
      ->default_val(this->split_threshold)
      ->group("Splitting and Naming");

  app.add_option("--splitdepth", this->split_depth,
                 "Number of repeat units to use in repeat splitting")
      ->default_val(this->split_depth)
      ->group("Splitting and Naming");

  app.add_option("--minsplitwin", this->min_split_window,
                 "Minimum repeat split window size")
      ->default_val(this->min_split_window)
      ->group("Splitting and Naming");
}

bool Settings::parse_input(int argc, const char **argv) {

  for (int i = 0; i < argc; ++i) {
    if (i > 0) {
      this->args += " ";
    }
    this->args += argv[i];

    // Argument preprocessor
    if (strlen(argv[i]) >= 3) {
      if (argv[i][0] == '-') {
        if (isalpha(argv[i][1])) {
          printf("Argument '%s' is not allowed (long arguments begin with --, "
                 "filenames may not begin with -)\n",
                 argv[i]);
          return false;
        }
      }
    }
  }

  CLI11_PARSE(this->app, argc, argv);
  bool passed = true;
  if (this->in_file.empty() && !this->show_memory) {
    printf("Input file required.\n");
    passed = false;
  }

  if (!this->json) {
    if (this->hide_seq) {
      printf("--hideseq is only available with --json\n");
      passed = false;
    }

    if (this->show_deltas) {
      printf("--showdelta is only available with --json\n");
      passed = false;
    }

    if (this->show_trace) {
      printf("--showtrace is only available with --json\n");
      passed = false;
    }

    if (this->show_logo_nums) {
      printf("--showlogo is only available with --json\n");
      passed = false;
    }
  }

  if (this->json && this->suppress_out) {
    printf("--suppress is incompatible with --json\n"
           "--suppress disables both JSON and BED output\n");
    passed = false;
  }

  if (this->mask_file.empty() && this->mask_with_n) {
    printf("--maskn requires an output file path provided to --mask\n");
    passed = false;
  }

  if (this->threads < 1) {
    printf("--threads must be at least 1.\n");
    passed = false;
  }

  if (this->window_size != -1 && this->window_size < this->max_period) {
    printf("--winsize cannot be smaller than --period\n");
    passed = false;
  }

  if (this->windows < this->threads) {
    printf("--windows must be at least as large as --threads\n");
    passed = false;
  }

  if (this->max_period < 1) {
    printf("--period must be at least 1\n");
    passed = false;
  }

  if (this->at < 0.0 || this->at > 1.0) {
    printf("--at must be > 0.0 and < 1.0\n");
    passed = false;
  }

  if (this->acgt.size() == 4) {
    float sum = 0.0;
    for (int i = 0; i < 4; ++i) {
      sum += this->acgt[i];
      if (this->acgt[i] < 0 || this->acgt[i] > 1.0) {
        printf("--acgt argument [%i] must be >= 0.0 and <= 1.0\n", i);
        passed = false;
      }
    }

    if (abs(1.0 - sum) > 0.001) {
      printf("--acgt arguments must sum to 1.0 ± 0.001\n");
      passed = false;
    }
  }

  if (this->match_probability < 0.0 || this->match_probability > 1.0) {
    printf("--match must be > 0.0 and < 1.0\n");
    passed = false;
  }

  if (this->period_decay <= 0.0 || this->period_decay > 1.0) {
    printf("--decay must be > 0.0 and <= 1.0\n");
    passed = false;
  }

  if (this->transition_nr < 0.0 || this->transition_nr > 1.0) {
    printf("--nr must be > 0.0 and < 1.0\n");
    passed = false;
  }

  if (this->transition_rn < 0.0 || this->transition_rn > 1.0) {
    printf("--rn must be > 0.0 and < 1.0\n");
    passed = false;
  }

  if (this->transition_ri < 0.0 || this->transition_ri > 1.0) {
    printf("--ri must be > 0.0 and < 1.0\n");
    passed = false;
  }

  if (this->transition_rd < 0.0 || this->transition_rd > 1.0) {
    printf("--rd must be > 0.0 and < 1.0\n");
    passed = false;
  }

  if (this->transition_ii < 0.0 || this->transition_ii > 1.0) {
    printf("--ii must be > 0.0 and < 1.0\n");
    passed = false;
  }

  if (this->transition_dd < 0.0 || this->transition_dd > 1.0) {
    printf("--dd must be > 0.0 and < 1.0\n");
    passed = false;
  }

  if (this->split_threshold <= 0) {
    printf("--splitval must be > 0.0\n");
    passed = false;
  }

  return passed;
}

int Settings::calculate_num_states() {
  int num_of_states = 0;

  // Calculate the number of cells needed
  for (int i = 0; i < this->max_period; ++i) {
    ++num_of_states;

    if (i <= this->max_insert) {
      if (i > 1) {
        num_of_states += i - 1;
      }
    }

    else {
      num_of_states += this->max_insert;
    }

    if (i <= this->max_delete) {
      if (i > 1) {
        num_of_states += i - 1;
      }
    }

    else {
      num_of_states += this->max_delete;
    }
  }

  if (this->read_all) {
    this->windows = -1;
  }

  return num_of_states;
}

void Settings::assign_settings() {
  if (!this->mask_file.empty()) {
    this->produce_mask = true;
  }

  // Assign default system setting parameters

  if (this->overlap == -1) {
    this->overlap = this->max_period;
  }

  if (this->window_size == -1) {
    int num_states = this->calculate_num_states();
    unsigned long cmem = num_states * this->max_period;

    // itty bitty models use ~40 mb per thread
    if (cmem <= 1000) {
      this->window_size = 10000 * this->max_period;
    }

    // tiny models use 80 mb per thread
    else if (cmem <= 10000) {
      this->window_size = 2000 * this->max_period;
    }

    // Small models use less than 160 mb per thread
    else if (cmem <= 200000) {
      this->window_size = 200 * this->max_period;
    }

    // medium models use less than 1 GB per thread
    else if (cmem <= 5000000) {
      this->window_size = 50 * this->max_period;
    }

    // Large models use less than 4 GB per thread
    else if (cmem <= 80000000) {
      this->window_size = 10 * this->max_period;
    }

    // Massive models use a lot of memory
    else {
      this->window_size = 2 * this->max_period;
    }
  }

  this->a_freq = this->at / 2.0;
  this->t_freq = this->at / 2.0;
  this->c_freq = (1.0 - this->at) / 2.0;
  this->g_freq = (1.0 - this->at) / 2.0;

  if (this->acgt.size() == 4) {
    float sum = 0.0;
    for (int i = 0; i < 4; ++i) {
      sum += this->acgt[i];
    }

    this->a_freq = this->acgt[0] / sum;
    this->c_freq = this->acgt[1] / sum;
    this->g_freq = this->acgt[2] / sum;
    this->t_freq = this->acgt[3] / sum;
  }
}

void Settings::print_memory_usage() {

  int num_states = this->calculate_num_states();
  unsigned long long dp_size = (unsigned long long)num_states;
  dp_size *= (unsigned long long)(this->window_size + (2 * this->overlap));
  printf("----------------------------\n");
  printf("Maximum repeat period: %llu\n", this->max_period);
  printf("Number of states: %i\n", num_states);
  printf("Total sequence window size: %lli\n",
         (this->window_size + 2 * this->overlap));
  printf("DP matrix cells: %llu\n", dp_size);
  printf("Threads: %i\n", this->threads);
  dp_size *= (unsigned long long)this->threads;
  dp_size *= 4;
  printf("Sequence window queue length: %llu\n", this->windows);
  dp_size += (this->window_size + (2 * this->overlap)) * this->windows;
  float gb_size = (float)dp_size / (1024.0 * 1024.0 * 1024.0);
  if (gb_size > 0.09) {
    printf("Total size: %.2f GB (%llu bytes)\n", gb_size, dp_size);
  } else {
    printf("Total size: (%llu bytes)\n", dp_size);
  }
  printf("----------------------------\n");

  printf("*Actual memory usage will be slightly greater due to "
         "output repeat queue and repeat split matrices\n");
}

// This could be done in less code with templates. Yikcy wicky yucky wucky.
std::string json_var(std::string name, std::string value) {
  return "\"" + name + "\": " + "\"" + value + "\",\n";
}

std::string json_var(std::string name, bool value) {
  if (value)
    return "\"" + name + "\": " + "true" + ",\n";
  else
    return "\"" + name + "\": " + "false" + ",\n";
}

std::string json_var(std::string name, int value) {
  return "\"" + name + "\": " + std::to_string(value) + ",\n";
}

std::string json_var(std::string name, float value) {
  return "\"" + name + "\": " + std::to_string(value) + ",\n";
}

std::string json_var(std::string name, unsigned long long value) {
  return "\"" + name + "\": " + std::to_string(value) + ",\n";
}

std::string json_var(std::string name, long long value) {
  return "\"" + name + "\": " + std::to_string(value) + ",\n";
}

#define JSONMACRO(NAME) json_string += json_var(#NAME, NAME)

std::string Settings::json_string() {
  std::string json_string = "";
  JSONMACRO(args);

  JSONMACRO(in_file);
  JSONMACRO(read_all);

  JSONMACRO(out_file);
  JSONMACRO(pval);
  JSONMACRO(pval_exponent_loc_m);
  JSONMACRO(pval_exponent_loc_b);
  JSONMACRO(pval_exponent_scale_m);
  JSONMACRO(pval_exponent_scale_b);

  JSONMACRO(json);
  JSONMACRO(hide_seq);
  JSONMACRO(show_deltas);
  JSONMACRO(show_trace);
  JSONMACRO(show_wid);
  JSONMACRO(show_logo_nums);
  JSONMACRO(suppress_out);

  JSONMACRO(produce_mask);
  JSONMACRO(mask_file);
  JSONMACRO(mask_with_n);

  JSONMACRO(threads);
  JSONMACRO(window_size);
  JSONMACRO(overlap);
  JSONMACRO(windows);

  JSONMACRO(min_score);
  JSONMACRO(min_length);
  JSONMACRO(min_units);

  JSONMACRO(max_period);
  JSONMACRO(max_insert);
  JSONMACRO(max_delete);

  JSONMACRO(a_freq);
  JSONMACRO(c_freq);
  JSONMACRO(g_freq);
  JSONMACRO(t_freq);
  JSONMACRO(match_probability);

  JSONMACRO(period_decay);
  JSONMACRO(transition_nr);
  JSONMACRO(transition_rn);
  JSONMACRO(transition_ri);
  JSONMACRO(transition_rd);
  JSONMACRO(transition_ii);
  JSONMACRO(transition_dd);

  JSONMACRO(no_split);
  JSONMACRO(max_split);
  JSONMACRO(split_threshold);
  JSONMACRO(split_depth);
  JSONMACRO(min_split_window);
  json_string.pop_back();
  json_string.pop_back();
  json_string += "\n";

  return json_string;
}

#undef JSONMACRO