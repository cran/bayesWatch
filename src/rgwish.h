// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
//     Copyright (C) 2012 - 2020  Reza Mohammadi                                                   |
//                                                                                                 |
//     This file is part of BDgraph package.                                                       |
//                                                                                                 |
//     BDgraph is free software: you can redistribute it and/or modify it under                    |
//     the terms of the GNU General Public License as published by the Free                        |
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
//                                                                                                 |
//     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                             |
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
  
#ifndef rgwish_H
#define rgwish_H

#include "matrix.h"
#include <stdio.h>
#include <string.h>

extern "C" {
  void rmvn_c( double rand_values[], double mus[], double K[], int p );
  
  void rwish_c( double Ts[], double K[], int *b, int *p );

  void rgwish_c( double G[], double D[], double K[], int *b, int *p, double *threshold, int *failed);

  void log_exp_mc( int G[], int nu[], int *b, double H[], int *check_H, int *mc, int *p, double f_T[] );
}
#endif
