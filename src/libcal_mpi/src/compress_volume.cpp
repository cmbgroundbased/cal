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

void cal::atm_sim::compress_volume()
{
    // Establish a mapping between full volume indices and observed
    // volume indices
    cal::Timer tm;
    tm.start();

    if ((rank == 0) && (verbosity > 0)) {
        std::cerr << "Compressing volume, N = " << nn << std::endl;
    }

    std::vector <unsigned char> hit;
    try {
        compressed_index.reset(new AlignedVector <long> (nn));
        std::fill(compressed_index->begin(), compressed_index->end(), -1);

        full_index.reset(new AlignedVector <long> (nn));
        std::fill(full_index->begin(), full_index->end(), -1);

        hit.resize(nn, false);
    } catch (...) {
        std::cerr << rank
                  << " : Failed to allocate element indices. nn = "
                  << nn << std::endl;
        throw;
    }
    // Start by flagging all elements that are hit

    for (long ix = 0; ix < nx - 1; ++ix) {
        if (ix % ntask != rank) continue;
        double x = xstart + ix * xstep;

        # pragma omp parallel for schedule(static, 10)
        for (long iy = 0; iy < ny - 1; ++iy) {
            double y = ystart + iy * ystep;

            for (long iz = 0; iz < nz - 1; ++iz) {
                double z = zstart + iz * zstep;
                if (in_cone(x, y, z)) {
# ifdef DEBUG
                    hit.at(ix * xstride + iy * ystride + iz * zstride) = true;
# else // ifdef DEBUG
                    hit[ix * xstride + iy * ystride + iz * zstride] = true;
# endif // ifdef DEBUG
                }
            }
        }
    }

    if ((rank == 0) && (verbosity > 0)) {
        std::cerr << "Flagged hits, flagging neighbors" << std::endl;
    }

    // For extra margin, flag all the neighbors of the hit elements

    std::vector <unsigned char> hit2 = hit;

    for (long ix = 1; ix < nx - 1; ++ix) {
        if (ix % ntask != rank) continue;

        # pragma omp parallel for schedule(static, 10)
        for (long iy = 1; iy < ny - 1; ++iy) {
            for (long iz = 1; iz < nz - 1; ++iz) {
                long offset = ix * xstride + iy * ystride + iz * zstride;

                if (hit2[offset]) {
                    // Flag this element but also its neighbours to facilitate
                    // interpolation

                    for (double xmul = -2; xmul < 4; ++xmul) {
                        if ((ix + xmul < 0) || (ix + xmul > nx - 1)) continue;

                        for (double ymul = -2; ymul < 4; ++ymul) {
                            if ((iy + ymul < 0) || (iy + ymul > ny - 1)) continue;

                            for (double zmul = -2; zmul < 4; ++zmul) {
                                if ((iz + zmul < 0) || (iz + zmul > nz - 1)) continue;

# ifdef DEBUG
                                hit.at(offset + xmul * xstride
                                       + ymul * ystride + zmul * zstride) = true;
# else // ifdef DEBUG
                                hit[offset + xmul * xstride
                                    + ymul * ystride + zmul * zstride] = true;
# endif // ifdef DEBUG
                            }
                        }
                    }
                }
            }
        }
    }

    hit2.resize(0);

    if ((rank == 0) && (verbosity > 0)) {
        std::cerr << "Creating compression table" << std::endl;
    }

    // Then create the mappings between the compressed and full indices

    long i = 0;
    for (long ifull = 0; ifull < nn; ++ifull) {
        if (hit[ifull]) {
            (*full_index)[i] = ifull;
            (*compressed_index)[ifull] = i;
            ++i;
        }
    }

    hit.resize(0);
    nelem = i;

    full_index->resize(nelem);

    tm.stop();

    if (rank == 0) {
        // if ( verbosity > 0 ) {
        tm.report("Volume compressed in");
        std::cout << i << " / " << nn << "(" << i * 100. / nn << " %)"
                  << " volume elements are needed for the simulation"
                  << std::endl
                  << "nx = " << nx << " ny = " << ny << " nz = " << nz
                  << std::endl
                  << "wx = " << wx << " wy = " << wy << " wz = " << wz
                  << std::endl;

        // }
    }

    if (nelem == 0) throw std::runtime_error("No elements in the observation cone.");
}
