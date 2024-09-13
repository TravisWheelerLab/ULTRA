#include "cli.hpp"
#include "mask.hpp"
#include "ultra.hpp"
#include <string>
#include <new>  // for std::bad_alloc


int main_wrapper(int argc, const char * argv[]) {
  // Prepare settings
  Settings *settings = new Settings();
  settings->prepare_settings();
  if (!settings->parse_input(argc, argv)) {
    exit(0);
  }

  settings->assign_settings();
  settings->print_memory_usage();
  if (settings->show_memory) {
    exit(0);
  }

  // Tuning stuff.
  std::vector<unsigned long long> coverage;
  std::vector<unsigned long long> shuffled_coverage;
  std::vector<std::string> param_strings;
  unsigned long long seq_length = 0;
  if (settings->tune) {

    // Get the parameters we are going to run with
    if (settings->tune_param_path.empty()) {
      if (settings->tune_medium)
        param_strings = medium_tune_settings();
      else if (settings->tune_large)
        param_strings = large_tune_settings();
      else
        param_strings = small_tune_settings();
    }else {
      param_strings = tune_settings_for_path(settings->tune_param_path);
    }

    // Run ULTRA on each parameter set
    printf("Performing parameter search...\n");
    printf("Coverage, False Discovery Rate, Parameters\n");

    int best_coverage_index = -1;
    double best_coverage = 0.0;
    for (auto arg_string : param_strings) {
      unsigned long cvg_tmp;
      double real_coverage = 0;
      double false_coverage = 0;

      Settings search_setting;
      search_setting.prepare_settings();
      search_setting.set_multi_option();
      if (!search_setting.parse_multi_input(argc, argv, arg_string)) {
        exit(0);
      }
      search_setting.assign_settings();

      if (!search_setting.tune_indels) {
        search_setting.max_insert = 0;
        search_setting.max_delete = 0;
      }
      search_setting.suppress_out = true;
      search_setting.hide_settings = true;
      search_setting.produce_mask = true;
      search_setting.no_split = true;
      search_setting.max_split = 0;
      search_setting.ultra_out = false;
      search_setting.json_out = false;
      search_setting.bed_out = false;

      // We first get the normal coverage
      auto *ultra = new Ultra(&search_setting);
      ultra->AnalyzeFile();
      ultra->OutputRepeats(true);
      seq_length = ultra->reader->fastaReader->total_seq_length;
      cvg_tmp = ultra->Coverage();
      real_coverage = (double)cvg_tmp / (double)seq_length;
      coverage.push_back(cvg_tmp);
      delete ultra;

      // Next we get our shuffled coverage
      ultra = new Ultra(&search_setting);
      ultra->shuffleSequence = true;
      ultra->AnalyzeFile();
      ultra->OutputRepeats(true);
      cvg_tmp = ultra->Coverage();
      false_coverage = (double)cvg_tmp / (double)seq_length;
      shuffled_coverage.push_back(cvg_tmp);
      delete ultra;

      double fdr = 1.0;
      if (real_coverage > 0.0) {
        fdr = false_coverage / real_coverage;
      }

      if (fdr <= settings->tune_fdr) {
        if (real_coverage > best_coverage) {
          best_coverage = real_coverage;
          best_coverage_index = coverage.size() - 1;
        }
      }

      printf("(%zu/%zu): %g, %g, %s\n",(coverage.size()),
             param_strings.size(),
             real_coverage, fdr, arg_string.c_str());
    }

    printf("-----------\n");
    if (best_coverage_index >= 0) {
      double real_coverage = (double)coverage[best_coverage_index] / (double)seq_length;
      double false_coverage = (double)shuffled_coverage[best_coverage_index] / (double)seq_length;
      double fdr = false_coverage / real_coverage;

      printf("Best coverage within FDR limit: %g, %g, %s\n", real_coverage,
             fdr,
             param_strings[best_coverage_index].c_str());

      delete settings;
      settings = new Settings();
      settings->prepare_settings();
      settings->set_multi_option();
      if (!settings->parse_multi_input(argc, argv, param_strings[best_coverage_index])) {
        exit(0);
      }
      settings->assign_settings();
    } else {
      printf("No parameters found within FDR limit.\n");
      exit(0);
    }

    if (settings->tune_only) {
      exit(0);
    }
  }

  // Perform the actual run


  auto *ultra = new Ultra(settings);
  ultra->AnalyzeFile();
  ultra->OutputRepeats(true);
  unsigned long long true_coverage = ultra->Coverage();
  unsigned long long shuff_coverage = 0;
  // Check to see if we are making a masked file
  if (settings->produce_mask) {
    FILE *f = fopen(settings->mask_file.c_str(), "w");
    OutputMaskedFASTA(settings->in_file, f, ultra->masks_for_seq,
                      settings->mask_with_n);
    fclose(f);
  }

  delete ultra;
  if (settings->estimate_fdr) {
    settings->suppress_out = true;
    settings->hide_settings = true;
    settings->produce_mask = true;
    settings->no_split = true;
    settings->max_split = 0;
    settings->ultra_out = false;
    settings->json_out = false;
    settings->bed_out = false;
    settings->out_file = "";
    ultra = new Ultra(settings);
    ultra->shuffleSequence = true;
    ultra->AnalyzeFile();
    ultra->OutputRepeats(true);
    shuff_coverage = ultra->Coverage();
    float fdr = (float)shuff_coverage / (float)true_coverage;
    printf("Estimated false discovery rate: %g\n", fdr);
  }

  return 0;
}


int main(int argc, const char *argv[]) {

  char *reserve_memory = (char *)malloc(65536);
  try {
    int r = main_wrapper(argc, argv);
    return r;
  }
  catch (const std::bad_alloc& e) {
    // This block is executed if memory allocation fails
    free(reserve_memory); // May be necessary in order to print
    std::cerr << "Memory allocation failed: " << e.what() << std::endl;
    std::cerr << "Your model may be too large to fit in memory"  << std::endl;
    std::cerr << "Try running: ultra --mem <your arguments> to see expected memory usage" << std::endl;
  }
  catch (const std::exception& e) {
    // This block is executed for any other standard exceptions
    std::cerr << "Standard exception caught: " << e.what() << std::endl;
  }
  catch (...) {
    // This block catches any other non-standard exceptions
    std::cerr << "Unknown exception caught!" << std::endl;
  }

  return -1;
}
