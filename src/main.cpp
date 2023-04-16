#include "settings.hpp"
#include "ultra.hpp"
#include "mask.h"
#include <string>
#include "cli.hpp"

int main(int argc, const char *argv[]) {

  auto settings = Settings_t();
  if (!settings.parse_input(argc, argv)) {
    exit(0);
  }

  settings.assign_settings();
  if (settings.show_memory) {
    settings.print_memory_usage();
    exit(0);
  }

/*
  auto settings = Settings(argc, argv);
  auto *ultra = new Ultra(&settings, 0);
  ultra->AnalyzeFile();
  ultra->OutputRepeats(true);

  if (!ultra->settings->v_outFilePath.empty()) {
    fclose(ultra->out);
  }

  FILE* f = fopen("./mask.fa", "w");
  if (settings.v_mask) {
    OutputMaskedFASTA(settings.v_filePath,
                      f,
                      ultra->masks_for_seq,
                      settings.v_maskWithN);
  }
  fclose(f);*/
  return 0;
}
