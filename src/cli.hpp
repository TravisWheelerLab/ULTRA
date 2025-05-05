//
// Created by Daniel Olson on 4/14/23.
//

#ifndef ULTRA_CLI_HPP
#define ULTRA_CLI_HPP

#define ULTRA_VERSION_STRING "1.1.0"
#define DEBUG_STRING ""
#ifdef DEBUG_PRAGMA
#undef DEBUG_STRING
#define DEBUG_STRING                                                           \
  " **WARNING** BUILT WITHOUT RELEASE OPTIMIZATION **WARNING** \n"
#endif

#include "../lib/CLI11.hpp"
#include <string>
#include <vector>

struct Settings {
  std::string args = "";

  // Are we displaying the citation text?
  bool cite = false;

  // Input settings
  std::string in_file = "";
  bool read_all = false;

  // Output settings
  std::string out_file = "";
  bool pval = false;
  bool disable_streaming_out = false;
  float p_value_loc = 4.27294921875;
  float p_value_scale = 1.8913602828979492;
  float p_value_freq = 0.000330256;
  bool estimate_fdr = false;

  int max_consensus_period = 1000000;
  bool ultra_out = false;
  bool bed_out = false;
  bool json_out = false;

  bool show_counts = false;
  bool show_seq = false;
  bool show_deltas = false;
  bool show_trace = false;
  bool show_wid = false;
  bool show_logo_nums = false;
  bool hide_settings = false;
  bool suppress_out = false;

  // Masking settings
  bool produce_mask = false;
  std::string mask_file = "";
  bool mask_with_n = false;

  // System settings
  int threads = 1;
  long long window_size = -1;
  long long overlap = -1;
  long long windows = 1024;
  bool show_memory = false;

  // Filter settings
  float min_score = -100.0;
  unsigned long long min_length = 10;
  unsigned long long min_units = 2;

  // Tuning parameters
  double tune_fdr = 0.05;
  bool tune = false;
  bool tune_medium = false;
  bool tune_large = false;
  bool tune_only = false;
  bool tune_indels = false;
  std::string tune_param_path;

  // Model parameters
  unsigned long long max_period = 100;
  unsigned long long max_insert = 10;
  unsigned long long max_delete = 10;

  // Probability parameters
  float a_freq = 0.25;
  float c_freq = 0.25;
  float g_freq = 0.25;
  float t_freq = 0.25;
  std::vector<float> acgt;
  float at = 0.5;
  float match_probability = 0.7;

  float period_decay = 0.85;
  float transition_nr = 0.01;
  float transition_rn = 0.05;
  float transition_ri = 0.02;
  float transition_rd = 0.02;
  float transition_ii = 0.02;
  float transition_dd = 0.02;

  // Split and naming parameters
  bool no_split = false;
  int max_split = 10;
  float split_threshold = 3.5;
  unsigned long long split_depth = 5;
  unsigned long long min_split_window = 15;
  // TODO
  //  implement these
  unsigned long long max_namable_period = 50;
  unsigned long long max_highfi_naming = 20;

  bool run_without_reader = false;

  CLI::App app{"\n"
               "=================================================\n"
               "(U)ltra (L)ocates (T)andemly (R)epetitive (A)reas\n"
               "     Daniel R. Olson and Travis J. Wheeler\n"
               "                 Version " ULTRA_VERSION_STRING
               "\n" DEBUG_STRING "     Use '--cite' for citation instructions\n"
               "=================================================\n"};

  void prepare_settings();
  void set_multi_option();
  bool parse_input(int argc, const char **argv);
  bool parse_multi_input(int argc, const char **argv, std::string arg_str);
  int calculate_num_states();
  void assign_settings();
  void print_memory_usage();
  std::string json_string();
};

std::vector<std::string> small_tune_settings();
std::vector<std::string> medium_tune_settings();
std::vector<std::string> large_tune_settings();
std::vector<std::string> tune_settings_for_path(std::string path);

void string_to_args(const std::string &str, int &argc, char **&argv);
std::pair<int, char **> combine_args(int argc1, const char **argv1, int argc2,
                                     char **argv2);

#endif // ULTRA_CLI_HPP
