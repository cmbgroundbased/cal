/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

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


double cal::mpi_atm_sim::interp(double x, double y, double z, std::vector <long> & last_ind,
              std::vector <double> & last_nodes)
{
    // Trilinear interpolation

    long ix = (x - xstart) * xstepinv;
    long iy = (y - ystart) * ystepinv;
    long iz = (z - zstart) * zstepinv;

    double dx = (x - (xstart + (double)ix * xstep)) * xstepinv;
    double dy = (y - (ystart + (double)iy * ystep)) * ystepinv;
    double dz = (z - (zstart + (double)iz * zstep)) * zstepinv;

# ifdef DEBUG
    if ((dx < 0) || (dx > 1) || (dy < 0) || (dy > 1) || (dz < 0) || (dz > 1)) {
        std::ostringstream o;
        o.precision(16);
        o << "atmsim::interp : bad fractional step: " << std::endl
          << "x = " << x << std::endl
          << "y = " << y << std::endl
          << "z = " << z << std::endl
          << "dx = " << dx << std::endl
          << "dy = " << dy << std::endl
          << "dz = " << dz << std::endl;
        std::cerr << o.str() << std::endl;
        throw std::runtime_error(o.str().c_str());
    }
# endif //ifdef DEBUG

    double c000, c001, c010, c011, c100, c101, c110, c111;

    if ((ix != last_ind[0]) || (iy != last_ind[1]) || (iz != last_ind[2])) {
# ifdef DEBUG
        if ((ix < 0) || (ix > nx - 2) || (iy < 0) || (iy > ny - 2)
            || (iz < 0) || (iz > nz - 2)) {
            std::ostringstream o;
            o.precision(16);
            o << "atmsim::interp : full index out of bounds at"
              << std::endl << "("
              << x << ", " << y << ", " << z << ") = ("
              << ix << "/" << nx << ", "
              << iy << "/" << ny << ", "
              << iz << "/" << nz << ")";
            std::cerr << o.str() << std::endl;
            throw std::runtime_error(o.str().c_str());
        }
# endif // ifdef DEBUG

        size_t offset = ix * xstride + iy * ystride + iz * zstride;

        size_t ifull000 = offset;
        size_t ifull001 = offset + zstride;
        size_t ifull010 = offset + ystride;
        size_t ifull011 = ifull010 + zstride;
        size_t ifull100 = offset + xstride;
        size_t ifull101 = ifull100 + zstride;
        size_t ifull110 = ifull100 + ystride;
        size_t ifull111 = ifull110 + zstride;

# ifdef DEBUG
        long ifullmax = compressed_index->size() - 1;
        if (
            (ifull000 < 0) || (ifull000 > ifullmax) ||
            (ifull001 < 0) || (ifull001 > ifullmax) ||
            (ifull010 < 0) || (ifull010 > ifullmax) ||
            (ifull011 < 0) || (ifull011 > ifullmax) ||
            (ifull100 < 0) || (ifull100 > ifullmax) ||
            (ifull101 < 0) || (ifull101 > ifullmax) ||
            (ifull110 < 0) || (ifull110 > ifullmax) ||
            (ifull111 < 0) || (ifull111 > ifullmax)) {
            std::ostringstream o;
            o.precision(16);
            o << "atmsim::observe : bad full index. "
              << "ifullmax = " << ifullmax << std::endl
              << "ifull000 = " << ifull000 << std::endl
              << "ifull001 = " << ifull001 << std::endl
              << "ifull010 = " << ifull010 << std::endl
              << "ifull011 = " << ifull011 << std::endl
              << "ifull100 = " << ifull100 << std::endl
              << "ifull101 = " << ifull101 << std::endl
              << "ifull110 = " << ifull110 << std::endl
              << "ifull111 = " << ifull111 << std::endl;
            std::cerr << o.str() << std::endl;
            throw std::runtime_error(o.str().c_str());
        }
# endif // ifdef DEBUG

        long i000 = (*compressed_index)[ifull000];
        long i001 = (*compressed_index)[ifull001];
        long i010 = (*compressed_index)[ifull010];
        long i011 = (*compressed_index)[ifull011];
        long i100 = (*compressed_index)[ifull100];
        long i101 = (*compressed_index)[ifull101];
        long i110 = (*compressed_index)[ifull110];
        long i111 = (*compressed_index)[ifull111];

# ifdef DEBUG
        long imax = realization->size() - 1;
        if (
            (i000 < 0) || (i000 > imax) ||
            (i001 < 0) || (i001 > imax) ||
            (i010 < 0) || (i010 > imax) ||
            (i011 < 0) || (i011 > imax) ||
            (i100 < 0) || (i100 > imax) ||
            (i101 < 0) || (i101 > imax) ||
            (i110 < 0) || (i110 > imax) ||
            (i111 < 0) || (i111 > imax)) {
            std::ostringstream o;
            o.precision(16);
            o << "atmsim::observe : bad compressed index. "
              << "imax = " << imax << std::endl
              << "i000 = " << i000 << std::endl
              << "i001 = " << i001 << std::endl
              << "i010 = " << i010 << std::endl
              << "i011 = " << i011 << std::endl
              << "i100 = " << i100 << std::endl
              << "i101 = " << i101 << std::endl
              << "i110 = " << i110 << std::endl
              << "i111 = " << i111 << std::endl
              << "(x, y, z) = " << x << ", " << y << ", " << z << ")"
              << std::endl
              << "in_cone(x, y, z) = " << in_cone(x, y, z)
              << std::endl;
            std::cerr << o.str() << std::endl;
            throw std::runtime_error(o.str().c_str());
        }
# endif // ifdef DEBUG

        c000 = (*realization)[i000];
        c001 = (*realization)[i001];
        c010 = (*realization)[i010];
        c011 = (*realization)[i011];
        c100 = (*realization)[i100];
        c101 = (*realization)[i101];
        c110 = (*realization)[i110];
        c111 = (*realization)[i111];

        last_ind[0] = ix;
        last_ind[1] = iy;
        last_ind[2] = iz;

        last_nodes[0] = c000;
        last_nodes[1] = c001;
        last_nodes[2] = c010;
        last_nodes[3] = c011;
        last_nodes[4] = c100;
        last_nodes[5] = c101;
        last_nodes[6] = c110;
        last_nodes[7] = c111;
    } else {
        c000 = last_nodes[0];
        c001 = last_nodes[1];
        c010 = last_nodes[2];
        c011 = last_nodes[3];
        c100 = last_nodes[4];
        c101 = last_nodes[5];
        c110 = last_nodes[6];
        c111 = last_nodes[7];
    }

    double c00 = c000 + (c100 - c000) * dx;
    double c01 = c001 + (c101 - c001) * dx;
    double c10 = c010 + (c110 - c010) * dx;
    double c11 = c011 + (c111 - c011) * dx;

    double c0 = c00 + (c10 - c00) * dy;
    double c1 = c01 + (c11 - c01) * dy;

    double c = c0 + (c1 - c0) * dz;

    return c;
}
