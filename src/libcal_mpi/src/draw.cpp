/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <CAL_MPI_AtmSim.hpp>
#include <math_rng.hpp>

/**
* Draw 10 000 gaussian variates to use in the drawing the
* simulation parameters
*
* Precalculate the ratio for covariance
*
* The wind is parallel to the surface, and it's indipendent
* to z. It's a good practice transform the coordinate
* in order to rotate the frame where the scan is across the
* X-axis.
*
* We take into account the opposite of the wind speed v = -v,
* in this way we can apply this to the telescope position.
* By example, S. Church 1995 (the paragraph of the wind)
*/
void cal::mpi_atm_sim::draw()
{
    // Draw 10 000 gaussian variates to use in the drawing the simulation parameters
    const uint64_t nrand = 10000;
    double randn[nrand];
    cal::rng_dist_normal(nrand, key1, key2, counter1, counter2, randn);
    counter2 += nrand;
    double * prand = randn;
    uint64_t irand = 0;

    if (rank == 0){
        lmin = 0;
        lmax = 0;
        w = -1;
        wdir = 0;
        z0 = 0;
        T0 = 0;

        while (lmin >= lmax){
            lmin = 0;
            lmax = 0;
            while (lmin <= 0 && irand < nrand -1){
                lmin = lmin_center + randn[irand++] * lmin_sigma;
            }
            while (lmax <= 0 && irand < nrand -1){
                lmax = lmax_center + randn[irand++] * lmax_sigma;
            }
        }

        while (w < 0 && irand < nrand - 1){
            w = w_center + randn[irand++] * w_sigma;
        }
        wdir = fmod(wdir_center + randn[irand++] * wdir_sigma, M_PI);
        while (z0 <= 0 && irand < nrand - 1){
            z0 = z0_center + randn[irand++] * z0_sigma;
        }
        while (T0 <= 0 && irand < nrand - 1){
            T0 = T0_center + randn[irand++] * T0_sigma;
        }

        if (irand == nrand) throw std::runtime_error("Failed to draw parameters in order to satisfy the boundary condictions");
    }

    if (MPI_Bcast(&lmin, 1, MPI_DOUBLE, 0, comm)) throw std::runtime_error(
                  "Failed to bcast lmin");

    if (MPI_Bcast(&lmax, 1, MPI_DOUBLE, 0, comm)) throw std::runtime_error(
                  "Failed to bcast lmax");

    if (MPI_Bcast(&w, 1, MPI_DOUBLE, 0, comm)) throw std::runtime_error(
                  "Failed to bcast w");

    if (MPI_Bcast(&wdir, 1, MPI_DOUBLE, 0, comm)) throw std::runtime_error(
                  "Failed to bcast wdir");

    if (MPI_Bcast(&z0, 1, MPI_DOUBLE, 0, comm)) throw std::runtime_error(
                  "Failed to bcast z0");

    if (MPI_Bcast(&T0, 1, MPI_DOUBLE, 0, comm)) throw std::runtime_error(
                  "Failed to bcast T0");

    // Precalculate the ratio for covariance
    z0inv = 1. / (2. * z0);

    // The wind is parallel to the surface. Here we rotate a frame where the scan is across the X-axis.

    double eastward_wind = w * cos(wdir);
    double northward_wind = w * sin(wdir);

    double angle = az0 - M_PI / 2;
    double wx_h = eastward_wind * cos(angle) - northward_wind * sin(angle);

    wy = eastward_wind * sin(angle) + northward_wind *cos(angle);

    wx = wx_h * cosel0;
    wz = -wx_h * sinel0;

    // Inverse the wind direction so we can apply it to the telescope position.

    wx = -wx;
    wy = -wy;

    if ((rank == 0) && (verbosity > 0)) {
        std::cerr << std::endl;
        std::cerr << "Atmospheric realization parameters:" << std::endl;
        std::cerr << " lmin = " << lmin << " m" << std::endl;
        std::cerr << " lmax = " << lmax << " m" << std::endl;
        std::cerr << "    w = " << w << " m/s" << std::endl;
        std::cerr << " easthward wind = " << eastward_wind << " m/s" << std::endl;
        std::cerr << " northward wind = " << northward_wind << " m/s" << std::endl;
        std::cerr << "  az0 = " << az0 * 180. / M_PI << " degrees" << std::endl;
        std::cerr << "  el0 = " << el0 * 180. / M_PI << " degrees" << std::endl;
        std::cerr << "   wx = " << wx << " m/s" << std::endl;
        std::cerr << "   wy = " << wy << " m/s" << std::endl;
        std::cerr << "   wz = " << wz << " m/s" << std::endl;
        std::cerr << " wdir = " << wdir * 180. / M_PI << " degrees" << std::endl;
        std::cerr << "   z0 = " << z0 << " m" << std::endl;
        std::cerr << "   T0 = " << T0 << " K" << std::endl;
    }

    return;

}
