/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <CALAtmSim.hpp>

/**
* Contructror for atm_sim class. This is the version NON-MPI.
* The class get from the Environment class the MAX_THRS_NUM and
* split the parallelization from them. This version of the library
* is not indicated from HPC-cluster but it's for personal use on your
* laptop.
*
* For a more cluster-oriented implementation see the class mpi_atm_sim.
*/
cal::atm_sim::atm_sim(double azmin, double azmax, double elmin, double elmax,
                        double tmin, double tmax,
                        double lmin_center, double lmin_sigma,
                        double lmax_center, double lmax_sigma,
                        double w_center, double w_sigma,
                        double wdir_center, double wdir_sigma,
                        double z0_center, double z0_sigma,
                        double T0_center, double T0_sigma,
                        double zatm, double zmax,
                        double xstep, double ystep, double zstep,
                        long nelem_sim_max,
                        int verbosity,
                        uint64_t key1, uint64_t key2,
                        uint64_t counterval1, uint64_t counterval2,
                        std::string cachedir,
                        double rmin, double rmax)
                        :
                        cachedir(cachedir),
                        verbosity(verbosity),
                        key1(key1), key2(key2),
                        counter1start(counterval1), counter2start(counterval2),
                        azmin(azmin), azmax(azmax),
                        elmin(elmin), elmax(elmax), tmin(tmin), tmax(tmax),
                        lmin_center(lmin_center), lmin_sigma(lmin_sigma),
                        lmax_center(lmax_center), lmax_sigma(lmax_sigma),
                        w_center(w_center), w_sigma(w_sigma),
                        wdir_center(wdir_center), wdir_sigma(wdir_sigma),
                        z0_center(z0_center), z0_sigma(z0_sigma),
                        T0_center(T0_center), T0_sigma(T0_sigma),
                        zatm(zatm), zmax(zmax),
                        xstep(xstep), ystep(ystep), zstep(zstep),
                        nelem_sim_max(nelem_sim_max),
                        rmin(rmin), rmax(rmax)
{

    std::cout << "CTOR atm_sim class" << std::endl;
    counter1 = counter1start;
    counter2 = counter2start;

    corrlim = 1e-3;

    ntask = 1;
    rank = 0;

    auto &  env = cal::Environment::get();
    nthread = env.max_threads();
    if ((rank == 0) && (verbosity > 0))
    {
        std::cerr << "atmsim constructed with " << ntask << " process, " << nthread << " threads per process." << std::endl;
    }

    /**
    * Check the sky patch boundary that you are going to observe
    * Not it works fine only for sky-patch. A lot of efforts are needed to
    * improve a continuos scanning strategy like QUIJOTE / STRIP / and others (LAT?)
    */
    if (azmin >= azmax) throw std::runtime_error("CTOR_atmsim ERROR: azmin >= azmax.");

    if (elmin < 0) throw std::runtime_error("CTOR_atmsim ERROR: elmin < 0, you are going to observe the ground ... ");

    if (elmax > M_PI_2) throw std::runtime_error("CTOR_atmsim ERROR: elmax > pi/2.");

    if (elmin > elmax) throw std::runtime_error("CTOR_atmsim ERROR: elmin > elmax.");

    if (tmin > tmax) throw std::runtime_error("CTOR_atmsim ERROR: tmin > tmax.");

    if (lmin_center > lmax_center) throw std::runtime_error("CTOR_atmsim ERROR: lmin_center > lmax_center");

    xstepinv = 1 / xstep;
    ystepinv = 1 / ystep;
    zstepinv = 1 / zstep;

    // Angular span observed
    delta_az = (azmax - azmin);
    delta_el = (elmax - elmin);
    delta_t = (tmax - tmin);

    // Starting point
    az0 = azmin + delta_az / 2;
    el0 = elmin + delta_el / 2;
    sinel0 = sin(el0);
    cosel0 = cos(el0);

    // Rotate the coordinate system. Align \hat{z} axis with \hat{r}.
    xxstep = xstep * cosel0 - zstep * sinel0;
    yystep = ystep;
    zzstep = xstep * sinel0 + zstep * cosel0;

    // speed up the in-cone calculation
    double tol = 0.1 * M_PI / 180; // 0.1 degree tolerance
    tanmin = tan(-0.5 * delta_az - tol);
    tanmax = tan(0.5 * delta_az + tol);

    // Some prints
    if ((rank == 0) && (verbosity > 0)) {
        std::cerr << std::endl;
        std::cerr << "Input parameters:" << std::endl;
        std::cerr << "             az = [" << azmin * 180. / M_PI << " - "
                  << azmax * 180. / M_PI << "] (" << delta_az * 180. / M_PI
                  << " degrees)" << std::endl;
        std::cerr << "             el = [" << elmin * 180. / M_PI << " - "
                  << elmax * 180. / M_PI << "] (" << delta_el * 180 / M_PI
                  << " degrees)" << std::endl;
        std::cerr << "              t = [" << tmin << " - " << tmax
                  << "] (" << delta_t << " s)" << std::endl;
        std::cerr << "           lmin = " << lmin_center << " +- " << lmin_sigma
                  << " m" << std::endl;
        std::cerr << "           lmax = " << lmax_center << " +- " << lmax_sigma
                  << " m" << std::endl;
        std::cerr << "              w = " << w_center << " +- " << w_sigma
                  << " m" << std::endl;
        std::cerr << "           wdir = " << wdir_center * 180. / M_PI << " +- "
                  << wdir_sigma * 180. / M_PI << " degrees " << std::endl;
        std::cerr << "             z0 = " << z0_center << " +- " << z0_sigma
                  << " m" << std::endl;
        std::cerr << "             T0 = " << T0_center << " +- " << T0_sigma
                  << " K" << std::endl;
        std::cerr << "           zatm = " << zatm << " m" << std::endl;
        std::cerr << "           zmax = " << zmax << " m" << std::endl;
        std::cerr << "       scan frame: " << std::endl;
        std::cerr << "          xstep = " << xstep << " m" << std::endl;
        std::cerr << "          ystep = " << ystep << " m" << std::endl;
        std::cerr << "          zstep = " << zstep << " m" << std::endl;
        std::cerr << " horizontal frame: " << std::endl;
        std::cerr << "         xxstep = " << xxstep << " m" << std::endl;
        std::cerr << "         yystep = " << yystep << " m" << std::endl;
        std::cerr << "         zzstep = " << zzstep << " m" << std::endl;
        std::cerr << "  nelem_sim_max = " << nelem_sim_max << std::endl;
        std::cerr << "        corrlim = " << corrlim << std::endl;
        std::cerr << "      verbosity = " << verbosity << std::endl;
        std::cerr << "           rmin = " << rmin << " m" << std::endl;
        std::cerr << "           rmax = " << rmax << " m" << std::endl;
    }

    // Cholesky decomposition
    chcommon = &cholcommon;
    cholmod_start(chcommon);
    if (verbosity > 1){
        chcommon->print = 3;
    }
    else{
        chcommon->print = 1;
    }
    chcommon->itype = CHOLMOD_INT;
    chcommon->dtype = CHOLMOD_DOUBLE;

    // The factorization is LL' no LDL'
    chcommon->final_ll = 1;
}

/**
* We don't neet a DTOR carefully definition. The unique_ptr(s) are authomatically
* free once they are derefenced.
*/
cal::atm_sim::~atm_sim()
{
    std::cout << "DTOR atm_sim class" << std::endl;
    compressed_index.reset();
    full_index.reset();
    realization.reset();
    cholmod_finish(chcommon);
}
