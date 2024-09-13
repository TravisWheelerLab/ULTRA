//
// Created by Daniel Olson on 4/14/23.
//

#include "cli.hpp"
#include <cstring>
#include <iostream>
#include <vector>

void Settings::prepare_settings() {

  app.footer("For additional information see README\n"
             "Supported by: NIH NIGMS P20GM103546 and NIH NHGRI U24HG010136\n");

  // *************
  // Input options
  // *************
  app.add_option("input_file", this->in_file, "DNA FASTA input file")
      ->required(false)
      ->group("Input");

  app.add_flag("-r,--read_all", this->read_all,
               "Read entire input file into memory"
               " (disables streaming input)")
      ->group("Input");
  app.add_flag("--readall", this->read_all,
               "Read entire input file into memory"
               " (disables streaming input)")
      ->group("");

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

  app.add_flag("--pval_loc", this->p_value_loc,
               "The exponential location used for converting scores to p-values.")
      ->group("Output");

  app.add_flag("--pval_scale", this->p_value_scale,
               "The exponential scale used for converting scores to p-values")
      ->group("Output");

  app.add_flag("--tsv", this->ultra_out,
               "Use TSV output format")
      ->group("Output");

  app.add_flag("--json", this->json_out,
               "Use JSON output format")
      ->group("Output");

  app.add_flag("--bed", this->bed_out,
               "Use BED output format")
      ->group("Output");

  app.add_option("--max_consensus", this->max_consensus_period,
                 "The maximum length of consensus pattern to include in output")
      ->default_val(this->max_consensus_period)
      ->group("Output");

  app.add_flag("--show_seq", this->show_seq,
               "Output repetitive region")
      ->group("Output");

  app.add_flag("--show_delta", this->show_deltas,
               "Show positional score deltas in JSON output")
      ->group("Output");
  app.add_flag("--showdelta", this->show_deltas,
               "Show positional score deltas in JSON output")
      ->group("");

  app.add_flag("--show_trace", this->show_trace,
               "Show Viterbi trace in JSON output")
      ->group("Output");
  app.add_flag("--showtrace", this->show_trace,
               "Show Viterbi trace in JSON output")
      ->group("");

  app.add_flag("--show_wid", this->show_wid,
               "Show sequence window IDs in JSON output")
      ->group("Output");
  app.add_flag("--showwid", this->show_wid,
               "Show sequence window IDs in JSON output")
      ->group("");

  app.add_flag("--show_logo", this->show_logo_nums,
               "Show logo numbers in JSON output")
      ->group("Output");
  app.add_flag("--showlogo", this->show_logo_nums,
               "Show logo numbers in JSON output")
      ->group("");

  app.add_flag("--hs, --hide_settings", this->hide_settings,
               "Do not output settings")
      ->group("Output");
  app.add_flag("--hidesettings", this->hide_settings,
               "Do not output settings")
      ->group("");

  app.add_flag("--suppress", this->suppress_out,
               "Do not output BED or JSON annotation")
      ->group("Output");

  app.add_flag("--fdr", this->estimate_fdr,
               "Estimate the False Discovery rate (runtime will be twice as long)")
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

  app.add_option("--win_size", this->window_size,
                 "Manually set sequence window size")
      ->group("System");
  app.add_option("--winsize", this->window_size)
      ->group("");

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

  app.add_option("--min_length", this->min_length,
                 "Minimum reportable repeat length")
      ->default_val(this->min_length)
      ->group("Filter");
  app.add_option("--minlen", this->min_length,
                 "Minimum reportable repeat length")
      ->default_val(this->min_length)
      ->group("");

  app.add_option("--min_unit", this->min_units,
                 "Minimum reportable number of repeat units")
      ->default_val(this->min_units)
      ->group("Filter");
  app.add_option("--minunit", this->min_units,
                 "Minimum reportable number of repeat units")
      ->default_val(this->min_units)
      ->group("");

  // *************
  // Tune options
  // *************


  app.add_flag("--tune", this->tune,
               "Tune parameters using a small search grid before running (see README)")
      ->group("Parameter Tuning");

  app.add_flag("--tune_medium", this->tune_medium,
               "Tune parameters before running (see README)")
      ->group("Parameter Tuning");

  app.add_flag("--tune_large", this->tune_large,
               "Tune parameters using a larger search grid before running (see README)")
      ->group("Parameter Tuning");

  app.add_option("--tune_file", this->tune_param_path,
                 "Use custom parameter search during tuning (see README)")
      ->default_val("")
      ->group("Parameter Tuning");

  app.add_flag("--tune_indel", this->tune_indels,
               "Enable indels while tuning")
      ->group("Parameter Tuning");
  app.add_flag("--tune_indels", this->tune_indels,
               "Enable indels while tuning")
      ->group("");

  app.add_option("--tune_fdr", this->tune_fdr,
                 "FDR to be tuned against (see README)")
      ->group("Parameter Tuning");

  app.add_flag("--tune_only", this->tune_only,
               "Tune parameters and don't run (see README)")
      ->group("Parameter Tuning");

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

  app.add_flag("--no_split", this->no_split, "Do not perform repeat splitting")
      ->group("Splitting and Naming");
  app.add_flag("--nosplit", this->no_split, "Do not perform repeat splitting")
      ->group("");

  app.add_option("--max_split", this->max_split,
                 "The maximum repeat period to perform repeat splitting")
      ->default_val(this->max_split)
      ->group("Splitting and Naming");
  app.add_option("--maxsplit", this->max_split,
                 "The maximum repeat period to perform repeat splitting")
      ->default_val(this->max_split)
      ->group("");

  app.add_option("--split_threshold", this->split_threshold, "Split threshold value")
      ->default_val(this->split_threshold)
      ->group("Splitting and Naming");
  app.add_option("--splitval", this->split_threshold, "Split threshold value")
      ->default_val(this->split_threshold)
      ->group("");

  app.add_option("--split_depth", this->split_depth,
                 "Number of repeat units to use in repeat splitting")
      ->default_val(this->split_depth)
      ->group("Splitting and Naming");
  app.add_option("--splitdepth", this->split_depth,
                 "Number of repeat units to use in repeat splitting")
      ->default_val(this->split_depth)
      ->group("");

  app.add_option("--min_split_window", this->min_split_window,
                 "Minimum repeat split window size")
      ->default_val(this->min_split_window)
      ->group("Splitting and Naming");
  app.add_option("--minsplitwin", this->min_split_window,
                 "Minimum repeat split window size")
      ->default_val(this->min_split_window)
      ->group("");

  app.add_flag("--cite", this->cite)->group("");
}

void Settings::set_multi_option() {
  for (auto &opt : app.get_options()) {
    opt->multi_option_policy(CLI::MultiOptionPolicy::TakeLast);
  }
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
          fprintf(stderr, "Argument '%s' is not allowed (long arguments begin with --, "
                          "filenames may not begin with -)\n",
                  argv[i]);
          return false;
        }
      }
    }
  }

  try {
    this->app.parse(argc, argv);
  } catch (const CLI::ExtrasError &e) {
    std::cerr << "Unrecognized flag or argument: " << e.what() << std::endl;
    exit(0); // or any other error handling
  } catch (const CLI::CallForHelp &e) {
    std::cout << this->app.help();
    exit(0); // or any other error handling
  }

  if (this->cite) {
    printf("BibTeX: \n"
           "@article {Olson2024ultra,\n"
           "  author = {Olson, Daniel R. and Wheeler, Travis J.},\n"
           "  title = {ULTRA-Effective Labeling of Repetitive Genomic Sequence},\n"
           "  elocation-id = {2024.06.03.597269},\n"
           "  year = {2024},\n"
           "  doi = {10.1101/2024.06.03.597269},\n"
           "  publisher = {Cold Spring Harbor Laboratory},\n"
           "  URL = {https://www.biorxiv.org/content/early/2024/06/04/2024.06.03.597269},\n"
           "  eprint = {https://www.biorxiv.org/content/early/2024/06/04/2024.06.03.597269.full.pdf},\n"
           "  journal = {bioRxiv}\n"
           "}\n");
    exit(0);
  }

  bool passed = true;
  if (this->in_file.empty() && !this->show_memory) {
    fprintf(stderr, "Input file required.\n");
    passed = false;
  }

  if (!this->json_out) {

    if (this->show_deltas) {
      fprintf(stderr, "--show_delta is only available with --json\n");
      passed = false;
    }

    if (this->show_trace) {
      fprintf(stderr, "--show_trace is only available with --json\n");
      passed = false;
    }

    if (this->show_logo_nums) {
      fprintf(stderr, "--show_logo is only available with --json\n");
      passed = false;
    }
  }

  if ((this->ultra_out || this->json_out || this->bed_out) && this->suppress_out) {
    fprintf(stderr, "--suppress is incompatible with --tsv, --json, and --bed\n");
    passed = false;
  }

  if (this->mask_file.empty() && this->mask_with_n) {
    fprintf(stderr, "--maskn requires an output file path provided to --mask\n");
    passed = false;
  }

  if (this->threads < 1) {
    fprintf(stderr, "--threads must be at least 1.\n");
    passed = false;
  }

  if (this->window_size != -1 && this->window_size < this->max_period) {
    fprintf(stderr, "--winsize cannot be smaller than --period\n");
    passed = false;
  }

  if (this->windows < this->threads) {
    fprintf(stderr, "--windows must be at least as large as --threads\n");
    passed = false;
  }

  if (this->max_period < 1) {
    fprintf(stderr, "--period must be at least 1\n");
    passed = false;
  }

  if (this->at < 0.0 || this->at > 1.0) {
    fprintf(stderr, "--at must be >= 0.0 and <= 1.0\n");
    passed = false;
  }

  if (this->acgt.size() == 4) {
    float sum = 0.0;
    for (int i = 0; i < 4; ++i) {
      sum += this->acgt[i];
      if (this->acgt[i] < 0 || this->acgt[i] > 1.0) {
        fprintf(stderr, "--acgt argument [%i] must be >= 0.0 and <= 1.0\n", i);
        passed = false;
      }
    }

    if (abs(1.0 - sum) > 0.001) {
      fprintf(stderr, "--acgt arguments must sum to 1.0 ± 0.001\n");
      passed = false;
    }
  }

  if (this->match_probability <= 0.0 || this->match_probability > 1.0) {
    fprintf(stderr, "--match must be > 0.0 and <= 1.0\n");
    passed = false;
  }

  if (this->period_decay <= 0.0 || this->period_decay > 1.0) {
    fprintf(stderr, "--decay must be > 0.0 and <= 1.0\n");
    passed = false;
  }

  if (this->transition_nr < 0.0 || this->transition_nr > 1.0) {
    fprintf(stderr, "--nr must be >= 0.0 and <= 1.0\n");
    passed = false;
  }

  if (this->transition_rn < 0.0 || this->transition_rn > 1.0) {
    fprintf(stderr, "--rn must be >= 0.0 and <= 1.0\n");
    passed = false;
  }

  if (this->transition_ri < 0.0 || this->transition_ri > 1.0) {
    fprintf(stderr, "--ri must be >= 0.0 and <= 1.0\n");
    passed = false;
  }

  if (this->transition_rd < 0.0 || this->transition_rd > 1.0) {
    fprintf(stderr, "--rd must be >= 0.0 and <= 1.0\n");
    passed = false;
  }

  if (this->transition_ii < 0.0 || this->transition_ii > 1.0) {
    fprintf(stderr, "--ii must be >= 0.0 and <= 1.0\n");
    passed = false;
  }

  if (this->transition_dd < 0.0 || this->transition_dd > 1.0) {
    fprintf(stderr, "--dd must be >= 0.0 and <= 1.0\n");
    passed = false;
  }

  if (this->split_threshold <= 0) {
    fprintf(stderr, "--splitval must be > 0.0\n");
    passed = false;
  }

  if (tune_fdr < 0.0 || tune_fdr > 1.0) {
    fprintf(stderr, "--tune_fdr must be >= 0 and <= 1.0\n");
    passed = false;
  }

  if (this->tune_only || this->tune_medium || this->tune_large || this->tune_indels) {
    this->tune = true;
  }

  if (this->tune_medium && this->tune_large) {
    fprintf(stderr, "Cannot use both --tune_medium and --tune_large\n");
    passed = false;
  }

  if (!this->tune_param_path.empty() && (this->tune_medium || this->tune_large)) {
    fprintf(stderr, "Cannot use both --tune_file and (--tune_small or --tune_large).\n");
    passed = false;
  }

  int c = 0;
  if (this->ultra_out) c++;
  if (this->json_out) c++;
  if (this->bed_out) c++;

  if (c > 1) {
    if (this->out_file.empty()) {
      fprintf(stderr, "Output file path must be provided when using multiple output formats .\n");
      passed = false;
    }
  }

  if (!(this->ultra_out || this->json_out || this->bed_out)) {
    this->ultra_out = true;
  }

  return passed;
}

bool Settings::parse_multi_input(int argc, const char **argv, std::string arg_str) {
  // Create combined arguments
  int new_argc;
  char **new_argv;
  string_to_args(arg_str, new_argc, new_argv);
  auto pair = combine_args(argc, argv, new_argc, new_argv);
  auto combined_argc = pair.first;
  auto combined_argv = pair.second;
  // Parse the combined arguments
  bool result = parse_input(combined_argc, (const char**)combined_argv);

  // Free the argument memory
  for (int i = 0; i < new_argc; ++i) {
    delete[] new_argv[i];
  }
  delete[] new_argv;

  for (int i = 0; i < combined_argc; ++i) {
    delete[] combined_argv[i];
  }
  delete[] combined_argv;

  return result;
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
    unsigned long period_memory = num_states * this->max_period;

    // itty bitty models use ~40 mb per thread
    if (period_memory <= 1000) {
      this->window_size = 10000 * this->max_period;
    }

    // tiny models use 80 mb per thread
    else if (period_memory <= 10000) {
      this->window_size = 2000 * this->max_period;
    }

    // Small models use less than 160 mb per thread
    else if (period_memory <= 200000) {
      this->window_size = 1000 * this->max_period;
    }

    // medium models use less than 1 GB per thread
    else if (period_memory <= 2000000) {
      this->window_size = 500 * this->max_period;
    }

    // medium models use less than 1 GB per thread
    else if (period_memory <= 20000000) {
      this->window_size = 100 * this->max_period;
    }

    else if (period_memory <= 50000000) {
      this->window_size = 50 * this->max_period;
    }

    // Large models use less than 4 GB per thread
    else  {
      this->window_size = 25 * this->max_period;
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

  if (this->no_split) {
    this->max_split = 0;
  }
}

void Settings::print_memory_usage() {

  int num_states = this->calculate_num_states();
  unsigned long long dp_size = (unsigned long long)num_states;
  dp_size *= (unsigned long long)(this->window_size + (2 * this->overlap));
  printf("----------------------------\n");
  printf("Maximum repeat period: %llu\n", this->max_period);
  printf("Number of model states: %i\n", num_states);
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
  printf("============================\n");
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

// TODO: Update json_string to include correct tuning info
std::string Settings::json_string() {
  std::string json_string = "";
  JSONMACRO(args);

  JSONMACRO(in_file);
  JSONMACRO(read_all);

  JSONMACRO(out_file);
  JSONMACRO(pval);
  JSONMACRO(p_value_loc);
  JSONMACRO(p_value_scale);
  JSONMACRO(p_value_freq);

  JSONMACRO(ultra_out);
  JSONMACRO(json_out);
  JSONMACRO(bed_out);
  JSONMACRO(max_consensus_period);
  JSONMACRO(show_seq);
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

std::vector<std::string> small_tune_settings()
{

  std::vector<std::string> settings;

  std::vector<float> match_settings = std::vector<float>{0.6, 0.75, 0.9};
  std::vector<float> at_settings = std::vector<float>{0.4, 0.5, 0.6};
  std::vector<float> repeat_start = std::vector<float>{0.001, 0.01};
  std::vector<float> repeat_stop = std::vector<float>{0.005, 0.05};

  for (auto match : match_settings) {
    for (auto at : at_settings) {
      for (int i = 0; i < 2; ++i) {
        std::string param_name = "-m " + std::to_string(match) + " ";
        param_name += "--at " + std::to_string(at) + " ";
        param_name += "--rn " + std::to_string(repeat_start[i]) + " ";
        param_name += "--nr " + std::to_string(repeat_stop[i]);
        settings.push_back(param_name);
      }
    }
  }

  return settings;
}

std::vector<std::string> medium_tune_settings()
{

  std::vector<std::string> settings;

  std::vector<float> match_settings = std::vector<float>{0.6, 0.7, 0.8, 0.9};
  std::vector<float> at_settings = std::vector<float>{0.3, 0.4, 0.5, 0.6, 0.7};
  std::vector<float> repeat_start = std::vector<float>{0.001, 0.01};
  std::vector<float> repeat_stop = std::vector<float>{0.005, 0.05};

  for (auto match : match_settings) {
    for (auto at : at_settings) {
      for (int i = 0; i < 2; ++i) {

        std::string param_name = "-m " + std::to_string(match) + " ";
        param_name += "--at " + std::to_string(at) + " ";
        param_name += "--rn " + std::to_string(repeat_start[i]) + " ";
        param_name += "--nr " + std::to_string(repeat_stop[i]);
        settings.push_back(param_name);
      }
    }
  }

  return settings;
}

std::vector<std::string> large_tune_settings()
{

  std::vector<std::string> settings;

  std::vector<float> match_settings = std::vector<float>{0.6, 0.7, 0.8, 0.9};
  std::vector<float> at_settings = std::vector<float>{0.3, 0.35, 0.4, 0.5, 0.6, 0.65, 0.7};
  std::vector<float> repeat_start = std::vector<float>{0.001, 0.005, 0.01};
  std::vector<float> repeat_stop = std::vector<float>{0.005, 0.01, 0.05};

  for (auto match : match_settings) {
    for (auto at : at_settings) {
      for (auto rep_start : repeat_start) {
        for (auto rep_stop : repeat_stop) {
          std::string param_name = "-m " + std::to_string(match) + " ";
          param_name += "--at " + std::to_string(at) + " ";
          param_name += "--rn " + std::to_string(rep_start) + " ";
          param_name += "--nr " + std::to_string(rep_stop);
          settings.push_back(param_name);
        }
      }
    }
  }

  return settings;
}

std::vector<std::string> tune_settings_for_path(std::string path) {
  std::vector<std::string> settings;
  std::ifstream file(path);

  // Check if the file is opened successfully
  if (!file.is_open()) {
    // Handle the error, for example, by logging or throwing an exception
    std::cerr << "Failed to open settings file: " << path << std::endl;
    exit(0);
  }

  std::string line;
  int line_num = 0;
  while (std::getline(file, line)) {
    ++line_num; // we 1 index line numbers
    // Check if the line is not empty
    if (!line.empty()) {
      int argc;
      char **argv;
      string_to_args(line, argc, argv);
      Settings *test_settings = new Settings();
      test_settings->prepare_settings();
      if (!test_settings->parse_input(argc, (const char**)argv)) {
        std::cerr << "Invalid arguments on line " << line_num << " in tune file. \"" << line << "\"" << std::endl;
        exit(0);
      }

      for (int i = 0; i < argc; ++i) {
        delete[] argv[i];
      }
      delete[] argv;
      delete test_settings;

      settings.push_back(line);
    }

  }

  file.close();
  return settings;
}



void string_to_args(const std::string& str, int& argc, char**& argv) {
  std::istringstream iss(str);
  std::vector<std::string> tokens;
  std::string token;

  // Tokenize the string
  while (iss >> token) {
    tokens.push_back(token);
  }

  // Set argc
  argc = tokens.size();

  // Allocate argv
  argv = new char*[argc + 1];

  // Copy tokens to argv
  for (int i = 0; i < argc; ++i) {
    argv[i] = new char[tokens[i].length() + 1];
    std::strcpy(argv[i], tokens[i].c_str());
  }

  // Null-terminate argv
  argv[argc] = nullptr;
}

std::pair<int, char**> combine_args(int argc1, const char** argv1, int argc2, char** argv2) {
  int combinedArgc = argc1 + argc2;
  char** combinedArgv = new char*[combinedArgc + 1];

  for (int i = 0; i < argc1; ++i) {
    combinedArgv[i] = new char[std::strlen(argv1[i]) + 1];
    std::strcpy(combinedArgv[i], argv1[i]);
  }

  for (int i = 0; i < argc2; ++i) {
    combinedArgv[argc1 + i] = new char[std::strlen(argv2[i]) + 1];
    std::strcpy(combinedArgv[argc1 + i], argv2[i]);
  }

  combinedArgv[combinedArgc] = nullptr;
  return std::make_pair(combinedArgc, combinedArgv);
}