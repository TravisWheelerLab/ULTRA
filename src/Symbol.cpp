//
//  Symbol.cpp
//  ultraP
//

#include "Symbol.hpp"

symbol SymbolForChar(char c) {

  switch (c) {
  case 'a':
  case 'A':
    return N_A;

  case 't':
  case 'T':
    return N_T;

  case 'c':
  case 'C':
    return N_C;

  case 'g':
  case 'G':
    return N_G;

  default:
    return N_UNK;
  }
}

char CharForSymbol(symbol s) {

  switch (s) {
  case N_A:
    return 'A';
  case N_T:
    return 'T';
  case N_G:
    return 'G';
  case N_C:
    return 'C';

  default:
    return 'N';
  }
}
