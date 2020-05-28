/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <iostream>
#include <CAL_MPI_AtmSim.hpp>
#include <atmosphere.hpp>
#include <mpi.h>
using namespace std;

int main(){

    MPI_Comm *comm;

    cal::mpi_atm_sim *atm_strip = new cal::mpi_atm_sim(azmin, azmax, elmin, elmax, tmin, tmax, lmin_center, lmin_sigma, lmax_center, lmax_sigma, w_center, w_sigma, wdir_center, wdir_sigma, z0_center, z0_sigma, T0_center, T0_sigma, zatm, zmax, xstep, ystep, zstep, nelem_sim_max, verbosity, MPI_COMM_WORLD, key1, key2, counterval1, counterval2, cachedir, rmin, rmax);


    // atm_strip->simulate(true);
    // atm_strop->observe() prova


    return 0;
}
