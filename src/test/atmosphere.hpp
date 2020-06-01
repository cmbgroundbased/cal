/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <string>

#ifndef ATM_TEST_HPP
#define ATM_TEST_HPP


double azmin = 0.0;
double azmax = M_PI/2;
double elmin = M_PI/4;
double elmax = M_PI/2.5;

double tmin = 0;
double tmax = 60;

double lmin_center = 0.01;
double lmin_sigma = 0.001;

double lmax_center = 10;
double lmax_sigma = 10;

double w_center = 0;
double w_sigma = 10;

double wdir_center = 0;
double wdir_sigma = 100;

double z0_center = 2000;
double z0_sigma = 0;

double T0_center = 280;
double T0_sigma = 10;

double zatm = 40000;
double zmax = 2000;

double xstep = 100;
double ystep = 100;
double zstep = 100;

long nelem_sim_max = 100;

int verbosity = 0;

uint64_t key1 = 0;
uint64_t key2 = 11;
uint64_t counterval1 = 0;
uint64_t counterval2 = 0;
std::string cachedir = std::string("/tmp/");

double rmin = 0;
double rmax = 10000;

#endif // ATM_TEST_HPP
