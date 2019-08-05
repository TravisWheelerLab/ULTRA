//
//  Symbol.hpp
//  ultraP
//


#ifndef Symbol_hpp
#define Symbol_hpp

#include <stdio.h>

typedef unsigned char symbol;


#define NUM_SYMBOLS     27
#define N_T             1
#define N_C             2
#define N_A             3
#define N_G             4
#define N_R             5
#define N_Y             6
#define N_K             7
#define N_M             8
#define N_S             9
#define N_W             10
#define N_B             11
#define N_D             12
#define N_H             13
#define N_V             14
#define N_N             15
#define N_GAP           16

#define N_UNK    0

symbol SymbolForChar(char c);
char CharForSymbol(symbol s);


#endif /* Symbol_hpp */
