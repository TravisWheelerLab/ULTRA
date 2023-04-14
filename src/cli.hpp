//
// Created by Daniel Olson on 4/14/23.
//

#ifndef ULTRA_CLI_HPP
#define ULTRA_CLI_HPP

#include <string>
#include <vector>
#include "../lib/CLI11.hpp"

struct Settings_t {
  std::string in_file = "";
  bool read_all = false;

  float min_score = -1000.0;
  unsigned long long min_length = 10;
  unsigned long long min_units = 2;

  std::string out_file = "";

  bool produce_mask = false;
  std::string mask_file = "";
  bool mask_with_n = false;

  bool suppress_out = false;

  bool pval = false;

  bool json = false;
  bool hide_seq = false;
  bool show_deltas = false;
  bool show_trace = false;
  bool show_wid = false;
  bool show_logo_nums = false;

  int threads = 1;
  long long window_size = -1;
  long long overlap = -1;
  unsigned long long windows = 1024;
  bool show_memory = false;

  unsigned long long max_period = 100;
  unsigned long long max_insert = 10;
  unsigned long long max_delete = 10;

  float a_freq = 0.3;
  float c_freq = 0.2;
  float g_freq = 0.2;
  float t_freq = 0.3;
  float match_probability = 0.8;

  float period_decay = 0.85;
  float transition_nr =0.01;
  float transition_rn = 0.05;
  float transition_ri = 0.02;
  float transition_rd = 0.02;
  float transition_ii = 0.02;
  float transition_dd = 0.02;

  bool no_split = false;
  unsigned long long max_split = 10;
  float split_threshold = 3.5;
  unsigned long long split_depth = 5;
  unsigned long long min_split_window = 10;

  CLI::App app{"(U)ltra (L)ocates (T)andemly (R)epetitive (A)reas\n"};

  Settings_t();

};

#endif // ULTRA_CLI_HPP
