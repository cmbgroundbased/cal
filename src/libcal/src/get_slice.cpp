/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <CALAtmSim.hpp>

void cal::atm_sim::get_slice(long & ind_start, long & ind_stop)
{
    // Identify a manageable slice of compressed indices to simulate next

    // Move element counter to the end of the most recent simulated slice
    ind_start = ind_stop;

    long ix_start = (*full_index)[ind_start] * xstrideinv;
    long ix1 = ix_start;
    long ix2;

    while (true) {
        // Advance element counter by one layer of elements
        ix2 = ix1;
        while (ix1 == ix2) {
            ++ind_stop;
            if (ind_stop == nelem) break;
            ix2 = (*full_index)[ind_stop] * xstrideinv;
        }

        // Check if there are no more elements
        if (ind_stop == nelem) break;

        // Check if we have enough to meet the minimum number of elements
        if (ind_stop - ind_start >= nelem_sim_max) break;

        // Check if we have enough layers
        // const int nlayer_sim_max = 10;
        // if ( ix2 - ix_start >= nlayer_sim_max ) break;
        ix1 = ix2;
    }

    if ((rank == 0) && (verbosity > 0)) {
        std::cerr << "X-slice: " << ix_start * xstep << " -- " << ix2 * xstep
                  << "(" << ix2 - ix_start <<  " " << xstep << " m layers)"
                  << " m out of  " << nx * xstep << " m"
                  << " indices " << ind_start << " -- " << ind_stop
                  << " ( " << ind_stop - ind_start << " )"
                  << " out of " << nelem << std::endl;
    }

    return;
}
