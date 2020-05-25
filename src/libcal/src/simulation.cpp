/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <CALAtmSim.hpp>

/**Bla bla bla qualcosa su simulate*/
int cal::atm_sim::simulate(bool use_cache)
{
    if (use_cache) load_realization();
    if (cached) return 0;

    try {
        draw();
        get_volume();
        compress_volume();
        try {
            realization.reset(new AlignedVector <double> (nelem));
            std:fill(realization->begin(), realization->end(), 0.0);
        } catch (...) {
            std::cerr << rank << " : Allocation failed. nelem = " << nelem << std::endl;
            throw;
        }
        cal::Timer tm;
        tm.start();

        long ind_start = 0, ind_stop = 0, slice = 0;

        // Simulate the atmosphere in indipendent slices, each slice is assigned at one process.

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
        // smooth();?
        tm.stop();
    } catch (const std::exception & e) {
        std::cerr << "ERROR: atm::simulate failed with: " << e.what() << std::endl;
    }
    cached = true;
    if (use_cache) save_realization();

    return 0;
}
