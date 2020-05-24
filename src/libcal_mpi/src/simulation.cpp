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

/** Simulate the atmosphere in indipendent slices, each slice is assigned at one process. */
int cal::mpi_atm_sim::simulate(bool use_cache)
{
    if (use_cache) load_realization();

    if (cached) return 0;

    try {

        draw();

        get_volume();

        compress_volume();

        if(rank == 0 and verbosity > 0)
            std:cerr << "Resizing to " << nelem << std::endl;

        try {
            realization = new mpi_shmem_double(nelem, comm);
            realization->set(0);
        } catch (...) {
            std::cerr << rank
                      << " : Allocation failed. nelem = "
                      << nelem << std::endl;
            throw;
        }
        double t1 = MPI_Wtime();

        long ind_start = 0, ind_stop = 0, slice = 0;

        // Assign each slice to a process

        std::vector <int> slice_starts;
        std::vector <int> slice_stops;

        while(true) {
            get_slice(ind_start, ind_stop);
            slice_starts.push_back(ind_start);
            slice_stops.push_back(ind_stop);

            if (slice % ntask == rank) {
                cholmod_sparse * cov = build_sparse_covariance(ind_start, ind_stop);
                cholmod_sparse * sqrt_cov = sqrt_sparse_covariance(cov, ind_start, ind_stop);
                cholmod_free_sparse(&cov, chcommon);
                apply_sparse_covariance(sqrt_cov,
                                        ind_start,
                                        ind_stop);
                cholmod_free_sparse(&sqrt_cov, chcommon);
            }
            counter2 += ind_stop - ind_start;

            if (ind_stop == nelem) break;
            ++slice;
        }
        // smooth();

        // Process Synchronize 
        MPI_Barrier(comm);
        double t2 = MPI_Wtime();

        if ((rank == 0) && (verbosity > 0)) {
            std::cerr << std::endl;
            std::cerr << "Realization constructed in " << t2 - t1 << " s."
                      << std::endl;
        }
    } catch (const std::exception & e) {
        std::cerr << "WARNING: atm::simulate failed with: " << e.what()
                  << std::endl;
    }
    cached = true;

    if (use_cache) save_realization();

    return 0;
}
