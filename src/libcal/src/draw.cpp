#include <CALAtmSim.hpp>
#include <sys_utils.hpp>
#include <sys_env.hpp>
#include <math_rng.hpp>
// #inluce <qualcosa per PRNG>

#include <sstream>
#include <iostream>
#include <fstream>
#include <cstring>
#include <random>    // Ha un sacco di generatori
#include <functional>
#include <cmath>
#include <algorithm> // per fare il std::sort

void cal::atm_sim::draw()
{
    // Draw 10 000 gaussian variates to use in the drawing the simulation parameters
    const size_t nrand = 10000;
    double randn[nrand];
    cal::rng_dist_normal(nrand, key1, key2, counter1, counter2, randn);
    counter2 += nrand;
    double * prand = randn;
    long irand = 0;

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

    return;

}
