/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <string>

#ifndef ATM_TEST_HPP
#define ATM_TEST_HPP

double azmin = (0/180)*M_PI;
double azmax = (359/180)*M_PI;
double elmin = (60/180)*M_PI;
double elmax = (80/180)*M_PI;

double tmin = 0;
double tmax_sim = 3600.0;
double tmax_tod = tmax_sim/2;

double lmin_center = 0.01;
double lmin_sigma = 0.001;

double lmax_center = 30;
double lmax_sigma = 10;

double w_center = 2;
double w_sigma = 0;

double wdir_center = -1.78;
double wdir_sigma = 0;

double z0_center = 2000;
double z0_sigma = 0;

double T0_center = 280;
double T0_sigma = 0;

double zatm = 40000;
double zmax = 2000;

double xstep = 20;
double ystep = 20;
double zstep = 20;

long nelem_sim_max = 100000;

int verbosity = 1;

uint64_t key1 = 0;
uint64_t key2 = 2<<32;
uint64_t counterval1 = 0;
uint64_t counterval2 = 0;
std::string cachedir = std::string("/tmp/");

double rmin = 0;
double rmax = 10000;

#endif // ATM_TEST_HPP
