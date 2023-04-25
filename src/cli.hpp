//
// Created by Daniel Olson on 4/14/23.
//

#ifndef ULTRA_CLI_HPP
#define ULTRA_CLI_HPP
#define ULTRA_VERSION_STRING "1.0.0 (beta 2)"
#include "../lib/CLI11.hpp"
#include <string>
#include <vector>

struct Settings {
  std::string args = "";

  // Input settings
  std::string in_file = "";
  bool read_all = false;

  // Output settings
  std::string out_file = "";
  bool pval = false;
  float pval_exponent_loc_m = -0.13;
  float pval_exponent_loc_b = 2.82;
  float pval_exponent_scale_m = 0.081;
  float pval_exponent_scale_b = 1.28;

  bool json = false;
  bool hide_seq = false;
  bool show_deltas = false;
  bool show_trace = false;
  bool show_wid = false;
  bool show_logo_nums = false;
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

  // Model parameters
  unsigned long long max_period = 100;
  unsigned long long max_insert = 10;
  unsigned long long max_delete = 10;

  // Probability parameters
  float a_freq = 0.3;
  float c_freq = 0.2;
  float g_freq = 0.2;
  float t_freq = 0.3;
  std::vector<int> acgt;
  float at = 0.6;
  float match_probability = 0.8;

  float period_decay = 0.85;
  float transition_nr = 0.01;
  float transition_rn = 0.05;
  float transition_ri = 0.02;
  float transition_rd = 0.02;
  float transition_ii = 0.02;
  float transition_dd = 0.02;

  // Split and naming parameters
  bool no_split = false;
  unsigned long long max_split = 10;
  float split_threshold = 3.5;
  unsigned long long split_depth = 5;
  unsigned long long min_split_window = 10;
  // TODO
  //  implement these
  unsigned long long max_namable_period = 50;
  unsigned long long max_highfi_naming = 20;

  CLI::App app{"\n"
               "=================================================\n"
               "(U)ltra (L)ocates (T)andemly (R)epetitive (A)reas\n"
               "     Daniel R. Olson and Travis J. Wheeler\n"
               "            Version " ULTRA_VERSION_STRING "\n"
               "=================================================\n"};

  void prepare_settings();
  bool parse_input(int argc, const char **argv);
  int calculate_num_states();
  void assign_settings();
  void print_memory_usage();
  std::string json_string();
};

#endif // ULTRA_CLI_HPP
