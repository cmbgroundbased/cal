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

/**
* Translate a compressed index into xyz-coordinates
* in the horizontal frame
* x, y, z are the coordinates in the scan frame,
* coord[0,1,2] are the coordinates in the horizontal frame
*/
void cal::mpi_atm_sim::ind2coord(long i, double * coord)
{
    long ifull = (*full_index)[i];

    long ix = ifull * xstrideinv;
    long iy = (ifull - ix * xstride) * ystrideinv;
    long iz = ifull - ix * xstride - iy * ystride;

    // coordinates in the scan frame
    double x = xstart + ix * xstep;
    double y = ystart + iy * ystep;
    double z = zstart + iz * zstep;

    // Into the horizontal frame
    coord[0] = x * cosel0 - z * sinel0;
    coord[1] = y;
    coord[2] = x * sinel0 + z * cosel0;
}

/**
* Translate scan frame xyz-coordinates into a compressed index
*/
long cal::mpi_atm_sim::coord2ind(double x, double y, double z)
{

    long ix = (x - xstart) * xstepinv;
    long iy = (y - ystart) * ystepinv;
    long iz = (z - zstart) * zstepinv;

# ifdef DEBUG
    if ((ix < 0) || (ix > nx - 1) || (iy < 0) || (iy > ny - 1) || (iz < 0) ||
        (iz > nz - 1)) {
        std::ostringstream o;
        o.precision(16);
        o << "atmsim::coord2ind : full index out of bounds at ("
          << x << ", " << y << ", " << z << ") = ("
          << ix << " /  " << nx << ", " << iy << " / " << ny << ", "
          << iz << ", " << nz << ")";
        throw std::runtime_error(o.str().c_str());
    }
# endif // ifdef DEBUG

    size_t ifull = ix * xstride + iy * ystride + iz * zstride;

    return (*compressed_index)[ifull];
}
