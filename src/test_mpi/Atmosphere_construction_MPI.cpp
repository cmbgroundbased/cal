/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <iostream>
#include <CAL_MPI_AtmSim.hpp>
#include <mpi_init.hpp>
#include <atmosphere.hpp>
#include <mpi.h>
#define DEBUG

int main(int argc, char * argv[]){

    std::cout << "Simulation start" << std::endl;
    double * t;
    double * az;
    double * el;
    double * tod;
    long nsamp;

    double fs_hz = 25; // Hz
    nsamp = long((tmax - tmin) * fs_hz);
    double ss = 0.017; // deg/sec
    bool direction_lr = true; // scanning direction

    cal::mpi_init(argc, argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);

    if (rank == 0){
        cal::Environment *env;
        env = &(env->get());
        env->print();

        t=(double *)malloc(nsamp*sizeof(double));
        az=(double *)malloc(nsamp*sizeof(double));
        el=(double *)malloc(nsamp*sizeof(double));
        tod=(double *)malloc(nsamp*sizeof(double));

        // Fill the elevation vector (el) with a CES at
        // el=cost=M_PI/3 // ok
        for(int i=0; i<nsamp; i++){
            el[i] = M_PI/3.0;
        }

        for(int i=0; i<nsamp; i++){
            tod[i] = 0.0;
        }

        // 15 deg of atmosphere slab.
        double d_az = ss / fs_hz;
        double az_now = azmin;
        for(int i=0; i<nsamp; i++){
            if(direction_lr){
                az_now += d_az;
                az[i] = az_now;
                if(az_now > 0.017*15-2*d_az) direction_lr = !direction_lr;
            } else {
                az_now -= d_az;
                az[i] = az_now;
                if(az_now < 0+2*d_az) direction_lr = !direction_lr;
            }
        }

        double dt = 1/fs_hz;
        for(int i=0; i<nsamp; i++){
            t[i] = tmin + dt*i;
        }

        std::cout << "Number of samples: " << nsamp << std::endl;
    }

    cal::mpi_atm_sim *atm_strip = new cal::mpi_atm_sim(azmin, azmax, elmin, elmax, tmin, tmax, lmin_center, lmin_sigma, lmax_center, lmax_sigma, w_center, w_sigma, wdir_center, wdir_sigma, z0_center, z0_sigma, T0_center, T0_sigma, zatm, zmax, xstep, ystep, zstep, nelem_sim_max, verbosity, comm, key1, key2, counterval1, counterval2, cachedir, rmin, rmax);


    atm_strip->simulate(true);
    atm_strip->observe(&(*t), &(*az), &(*el), &(*tod), nsamp);

    // if(rank==0){
    //     for(int i=0; i<nsamp; i++){
    //         std::cout << t[i] << " " << tod[i] << std::endl;
    //     }
    // }

    MPI_Finalize();
    std::cout << "#Time\tElevation\tAzimuth\tTOD" << std::endl;
    for(int i=0; i<nsamp; i++){
        std::cout << t[i] << " " << el[i] << " " << az[i] << " " <<  tod[i]/tod[1] << std::endl;
    }

    return 0;
}
