#include "settings.hpp"
#include "ultra.hpp"

#include <string>

int main(int argc, const char *argv[]) {
  auto settings = Settings(argc, argv);
  auto *ultra = new Ultra(&settings, 0);
  ultra->AnalyzeFile();
  ultra->OutputRepeats(true);

  if (!ultra->settings->v_outFilePath.empty()) {
    fclose(ultra->out);
  }

  return 0;
}
