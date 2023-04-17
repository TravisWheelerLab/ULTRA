#include "cli.hpp"
#include "mask.hpp"
#include "ultra.hpp"
#include <string>

int main(int argc, const char *argv[]) {

  auto settings = Settings();
  if (!settings.parse_input(argc, argv)) {
    exit(0);
  }

  settings.assign_settings();
  if (settings.show_memory) {
    settings.print_memory_usage();
    exit(0);
  }

  auto *ultra = new Ultra(&settings, 0);
  ultra->AnalyzeFile();
  ultra->OutputRepeats(true);

  if (!settings.out_file.empty()) {
    fclose(ultra->out);
  }

  if (settings.produce_mask) {
    FILE *f = fopen(settings.mask_file.c_str(), "w");
    OutputMaskedFASTA(settings.in_file, f, ultra->masks_for_seq,
                      settings.mask_with_n);
    fclose(f);
  }
  return 0;
}
