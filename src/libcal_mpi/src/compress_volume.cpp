/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <CAL_MPI_AtmSim.hpp>
#include <vector>
/**
* Establish a mapping between full volume indices and observed
* volume indices.
*/
void cal::mpi_atm_sim::compress_volume()
{
    double t1 = MPI_Wtime();

    if ((rank == 0) && (verbosity > 0)) {
        std::cerr << "Compressing volume, N = " << nn << std::endl;
    }

    std::vector <unsigned char> hit;
    try {
        compressed_index = new mpi_shmem_long(nn, comm);
        compressed_index->set(-1);

        full_index = new mpi_shmem_long(nn, comm);
        full_index->set(-1);

        hit.resize(nn, false);
    } catch (...) {
        std::cerr << rank
                  << " : Failed to allocate element indices. nn = "
                  << nn << std::endl;
        throw;
    }

    // Start by flagging all elements that are hit
    for (int64_t ix = 0; ix < nx - 1; ++ix) {
        if (ix % ntask != rank){
            continue;
        }
        double x = xstart + ix * xstep;

        # pragma omp parallel for schedule(static, 10)
        for (int64_t iy = 0; iy < ny - 1; ++iy) {
            double y = ystart + iy * ystep;

            for (int64_t iz = 0; iz < nz - 1; ++iz) {
                double z = zstart + iz * zstep;
                if (in_cone(x, y, z)) {
                    hit[ix * xstride + iy * ystride + iz * zstride] = true;
                }
            }
        }
    }

    if ((rank == 0) && (verbosity > 0)) {
        std::cerr << "Flagged hits, flagging neighbors" << std::endl;
    }

    // For extra margin, flag all the neighbors of the
    // hit elements
    if (MPI_Allreduce(MPI_IN_PLACE, hit.data(), (int)nn,
                      MPI_UNSIGNED_CHAR, MPI_LOR, comm)) throw std::runtime_error(
                  "Failed to gather hits");

    std::vector <unsigned char> hit2 = hit;

    for (int64_t ix = 1; ix < nx - 1; ++ix) {
        if (ix % ntask != rank) continue;

        # pragma omp parallel for schedule(static, 10)
        for (int64_t iy = 1; iy < ny - 1; ++iy) {
            for (int64_t iz = 1; iz < nz - 1; ++iz) {
                int64_t offset = ix * xstride + iy * ystride + iz * zstride;

                if (hit2[offset]) {
                    // Flag this element but also its
                    // neighbours to facilitate
                    // interpolation

                    for (int64_t xmul = -2; xmul < 4; ++xmul) {
                        if ((ix + xmul < 0) || (ix + xmul > nx - 1)) continue;

                        for (int64_t ymul = -2; ymul < 4; ++ymul) {
                            if ((iy + ymul < 0) || (iy + ymul > ny - 1)) continue;

                            for (int64_t zmul = -2; zmul < 4; ++zmul) {
                                if ((iz + zmul < 0) || (iz + zmul > nz - 1)) continue;
                                hit[offset + xmul * xstride + ymul * ystride + zmul * zstride] = true;
                            }
                        }
                    }
                }
            }
        }
    }

    hit2.resize(0);

    if (MPI_Allreduce(MPI_IN_PLACE, hit.data(), (int)nn,
                      MPI_UNSIGNED_CHAR, MPI_LOR, comm)) throw std::runtime_error(
                  "Failed to gather hits");

    if ((rank == 0) && (verbosity > 0)) {
        std::cerr << "Creating compression table" << std::endl;
    }

    // Then create the mappings between the compressed and
    // full indices
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

    double t2 = MPI_Wtime();

    if (rank == 0) {
        std::cout << "Volume compressed in " << t2 - t1
          << " s." << std::endl
          << i << " / " << nn
          << "(" << i * 100. / nn << " %)"
          << " volume elements are needed for the simulation"
          << std::endl
          << "nx = " << nx << " ny = " << ny << " nz = " << nz
          << std::endl
          << "wx = " << wx << " wy = " << wy << " wz = " << wz
          << std::endl;
    }

    if (nelem == 0)
        throw std::runtime_error("No elements in the observation cone.");
}
