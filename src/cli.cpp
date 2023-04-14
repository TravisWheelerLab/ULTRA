//
// Created by Daniel Olson on 4/14/23.
//

#include "cli.hpp"

Settings_t::Settings_t() {
  app.add_option("input_file",
                 this->in_file,
                 "DNA FASTA input file");
  app.add_flag("-r,--readall",
               this->read_all,
               "Read entire input file into memory (disables streaming input)");
}