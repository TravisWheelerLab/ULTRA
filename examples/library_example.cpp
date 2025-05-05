//
// Created by Daniel Olson on 5/5/25.
//
#include <cstdio>
#include <ultra.hpp>
#include <cli.hpp>

int main() {
  printf("Note that the Makefile is currently linking against ULTRA/build/libultra_core.a\n");
  printf("   and is using ULTRA/src/ as the source of ULTRA header files.\n");
  printf("The location of the ULTRA lib and header files will need to be adjusted in real use-cases.\n\n\n");

  printf("Starting test of library.\n");


  // These are the settings we will use with ULTRA
  Settings *settings = new Settings();
  settings->run_without_reader = true; // This must be set to true when using as lib
  settings->window_size = 1000000; // window size should be equal to the largest sequence you will analyze
  settings->overlap = 0; // overlap should be turned off


  // We can also pass in an argument string to settings like so:
  int argc;
  char **argv;
  std::string arg_string = "-p 25 -i 3 -d 3";
  string_to_args(arg_string, argc, argv);

  settings->prepare_settings(); // This should be called before running settings->parse_input()
  if (!settings->parse_input(argc, (const char**)argv)) {
    exit(0);
  }

  // Finally, we call settings->assign_settings()
  settings->assign_settings();

  // We create a reusable ULTRA object
  auto ultra = new Ultra(settings);
  // We can find repeats in a C++ string using "FindRepeatsInString(std::string)"

  printf("Finding repeats in string ``aaaaaaaaaaaaaaaaaaaaaaaaaa''\n");
  auto repeats = ultra->FindRepeatsInString("aaaaaaaaaaaaaaaaaaaaaaaaaa");

  printf("Start Length Pattern\n");
  for (int i= 0; i < repeats->size(); ++i) {
    auto r = repeats->at(i);
    printf("%lu %lu %i %s\n", r->windowStart, r->repeatLength, r->repeatPeriod, r->GetConsensus().c_str());
    delete r; // We are in charge of the memory management of the repeats returned by FindRepeatsInString
  }

  delete repeats; // We are also in charge of the memory for the repeat array itself
  printf("---------\n");


  printf("Finding repeats in string ``aggtaaggtaaggtaaggtaaggtaagcggtataacatacagatctgactactactactactactactactactac''\n");
  repeats = ultra->FindRepeatsInString("aggtaaggtaaggtaaggtaaggtaagcggtataacatacagatctgactactactactactactactactactac");
  printf("Start Length Pattern\n");
  for (int i= 0; i < repeats->size(); ++i) {
    auto r = repeats->at(i);
    printf("%lu %lu %i %s\n", r->windowStart, r->repeatLength, r->repeatPeriod, r->GetConsensus().c_str());
    delete r;
  }

  delete repeats;
  printf("---------\n");

  printf("All positions reported are 0 indexed.\n");

  return 0;
}