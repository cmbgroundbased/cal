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

double median(std::vector <double> vec) {
    if (vec.size() == 0) return 0;

    std::sort(vec.begin(), vec.end());
    int half1 = (vec.size() - 1) * .5;
    int half2 = vec.size() * .5;

    return .5 * (vec[half1] + vec[half2]);
}

double mean(std::vector <double> vec) {
    if (vec.size() == 0) return 0;

    double sum = 0;
    for (auto & val : vec) sum += val;

    return sum / vec.size();
}

/**
* Replace each vertex with a mean of its first-N-neighbors.
* The smoothing length is driven by the w parameter.
* In this implementation w = 3
*/
void cal::mpi_atm_sim::smooth()
{
    double t1 = MPI_Wtime();

    double coord[3];

    std::vector <double> smoothed_realization(realization->size());

    for (size_t i = 0; i < full_index->size(); ++i) {
        ind2coord(i, coord);
        long ix = coord[0] * xstepinv;
        long iy = coord[1] * ystepinv;
        long iz = coord[2] * zstepinv;

        long offset = ix * xstride + iy * ystride + iz * zstride;

        long w = 3;
        long ifullmax = compressed_index->size();

        std::vector <double> vals;

        for (int xoff = 0; xoff <= 0; ++xoff) {
            if (ix + xoff < 0) continue;
            if (ix + xoff >= nx) break;

            for (int yoff = -w; yoff <= w; ++yoff) {
                if (iy + yoff < 0) continue;
                if (iy + yoff >= ny) break;

                for (int zoff = -w; zoff <= w; ++zoff) {
                    if (iz + zoff < 0) continue;
                    if (iz + zoff >= nz) break;

                    long ifull = offset + xoff * xstride + yoff * ystride
                                 + zoff * zstride;

                    if ((ifull < 0) || (ifull >= ifullmax)) throw std::runtime_error(
                                  "Index out of range in smoothing.");

                    long ii = (*compressed_index)[ifull];

                    if (ii >= 0) {
                        vals.push_back((*realization)[ii]);
                    }
                }
            }
        }

        // Get the smoothed value

        smoothed_realization[i] = mean(vals);
    }

    if (realization->rank() == 0) {
        for (int i = 0; i < realization->size(); ++i) {
            (*realization)[i] = smoothed_realization[i];
        }
    }

    double t2 = MPI_Wtime();

    if ((rank == 0) && (verbosity > 0)) {
        std::cerr << "Realization smoothed in " << t2 - t1 << " sec." << std::endl;
    }
    return;
}
