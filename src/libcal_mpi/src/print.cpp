/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <cal_mpi_internal.hpp>

void cal::mpi_atm_sim::print(std::ostream & out) const
{
    for (int i = 0; i < ntask; ++i) {
        MPI_Barrier(comm);
        if (rank != i) continue;
        out << rank << " : cachedir " << cachedir << std::endl;
        out << rank << " : ntask = " << ntask
            << ", nthread = " << nthread << std::endl;
        out << rank << " : verbosity = " << verbosity
            << ", key1 = " << key1
            << ", key2 = " << key2
            << ", counter1 = " << counter1
            << ", counter2 = " << counter2
            << ", counter1start = " << counter1start
            << ", counter2start = " << counter2start << std::endl;
        out << rank << " : azmin = " << azmin
            << ", axmax = " << azmax
            << ", elmin = " << elmin
            << ", elmax = " << elmax
            << ", tmin = " << tmin
            << ", tmax = " << tmax
            << ", sinel0 = " << sinel0
            << ", cosel0 = " << cosel0
            << ", tanmin = " << tanmin
            << ", tanmax = " << tanmax << std::endl;
        out << rank << " : lmin_center = " << lmin_center
            << ", lmax_center = " << lmax_center
            << ", w_center = " << w_center
            << ", w_sigma = " << w_sigma
            << ", wdir_center = " << wdir_center
            << ", wdir_sigma = " << wdir_sigma
            << ", z0_center = " << z0_center
            << ", z0_sigma = " << z0_sigma
            << ", T0_center = " << T0_center
            << ", T0_sigma = " << T0_sigma
            << ", z0inv = " << z0inv << std::endl;
    }
    MPI_Barrier(comm);
}
