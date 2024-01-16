#include "cli.hpp"
#include "mask.hpp"
#include "ultra.hpp"
#include <string>

int main(int argc, const char *argv[]) {

  // Prepare settings
  Settings settings;
  settings.prepare_settings();
  if (!settings.parse_input(argc, argv)) {
    exit(0);
  }

  settings.assign_settings();
  if (settings.show_memory) {
    settings.print_memory_usage();
    exit(0);
  }

  // Tuning stuff.
  std::vector<unsigned long long> coverage;
  std::vector<unsigned long long> shuffled_coverage;
  std::vector<std::string> param_strings;
  std::vector<std::tuple<std::string,Settings *>> search_params;
  unsigned long long seq_length = 0;
  if (settings.tune) {

    // Get the parameters we are going to run with
    if (settings.tune_param_path.empty()) {
      if (settings.tune_small)
        search_params = small_tune_settings(argc, argv);
      else if (settings.tune_large)
        search_params = large_tune_settings(argc, argv);
      else
        search_params = default_tune_settings(argc, argv);
    }else {
      search_params = tune_settings_for_path(settings.tune_param_path);
    }

    // Run ULTRA on each parameter set
    printf("Performing parameter search...\n");
    printf("Coverage, False Discovery Rate, Parameters\n");

    int best_coverage_index = -1;
    double best_coverage = 0.0;
    for (auto [params, setting] : search_params) {
      unsigned long cvg_tmp;
      double real_coverage = 0;
      double false_coverage = 0;

      // Push back our params
      param_strings.push_back(params);

      // We first get the normal coverage
      auto *ultra = new Ultra(setting);
      ultra->AnalyzeFile();
      ultra->OutputRepeats(true);
      seq_length = ultra->reader->fastaReader->total_seq_length;
      cvg_tmp = ultra->Coverage();
      real_coverage = (double)cvg_tmp / (double)seq_length;
      coverage.push_back(cvg_tmp);
      delete ultra;

      // Next we get our shuffled coverage
      ultra = new Ultra(setting);
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

      if (fdr <= settings.tune_fdr) {
        if (real_coverage > best_coverage) {
          best_coverage = real_coverage;
          best_coverage_index = param_strings.size() - 1;
        }
      }

      printf("(%zu/%zu): %g, %g, %s\n",(param_strings.size()),
             search_params.size(),
             real_coverage, fdr, params.c_str());
    }

    printf("-----------\n");
    if (best_coverage_index >= 0) {
      double real_coverage = (double)coverage[best_coverage_index] / (double)seq_length;
      double false_coverage = (double)shuffled_coverage[best_coverage_index] / (double)seq_length;
      double fdr = false_coverage / real_coverage;

      printf("Best coverage within FDR limit: %g, %g, %s\n", real_coverage,
             fdr,
             param_strings[best_coverage_index].c_str());
    } else {
      printf("No parameters found within FDR limit.\n");
      exit(0);
    }

    if (settings.tune_only) {
      exit(0);
    }
  }

  // Perform the actual run
  auto *ultra = new Ultra(&settings);
  ultra->AnalyzeFile();
  ultra->OutputRepeats(true);

  // Check to see if we are making a masked file
  if (settings.produce_mask) {
    FILE *f = fopen(settings.mask_file.c_str(), "w");
    OutputMaskedFASTA(settings.in_file, f, ultra->masks_for_seq,
                      settings.mask_with_n);
    fclose(f);
  }

  delete ultra;

  return 0;
}
