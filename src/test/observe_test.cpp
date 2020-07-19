/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <iostream>
#include <CALAtmSim.hpp>
#include <atmosphere.hpp>

int main(int argc, char * argv[]){

    std::cout << "Simulation start" << std::endl;
    double * t;
    double * az;
    double * el;
    double * tod;
    long nsamp;

    double fs_hz = 25; // Hz
    nsamp = long((tmax_tod - tmin) * fs_hz);
    double ss = 0.017; // deg/sec
    bool direction_lr = true; // scanning direction


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
        el[i] = (M_PI/3.0)+(M_PI/24.0);
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
        if(az_now > 0.017*15-d_az) direction_lr = !direction_lr;
    } else {
        az_now -= d_az;
        az[i] = az_now;
        if(az_now < 0+d_az) direction_lr = !direction_lr;
        }
    }

    double dt = 1/fs_hz;
    for(int i=0; i<nsamp; i++){
    t[i] = tmin + dt*i;
    }

    std::cout << "Number of samples: " << nsamp << std::endl;

    cal::atm_sim *atm_strip = new cal::atm_sim(azmin, azmax, elmin, elmax, tmin, tmax_sim, lmin_center, lmin_sigma, lmax_center, lmax_sigma, w_center, w_sigma, wdir_center, wdir_sigma, z0_center, z0_sigma, T0_center, T0_sigma, zatm, zmax, xstep, ystep, zstep, nelem_sim_max, verbosity, key1, key2, counterval1, counterval2, cachedir, rmin, rmax);


    atm_strip->simulate(false);


    atm_strip->observe(&(*t), &(*az), &(*el), &(*tod), nsamp, -1.0);

    return 0;
}
