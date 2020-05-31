
// Copyright (c) 2015-2020 by the parties listed in the AUTHORS file.
// All rights reserved.  Use of this source code is governed by
// a BSD-style license that can be found in the LICENSE file.

#include <mpi_init.hpp>

/**
* If MPI is not yet initialized (by mpi4py or some
* other place), then initialize it here.
*/
void cal::mpi_init(int argc, char * argv[])
{
    std::cout << "Initialize MPI" << std::endl;
    int ret;
    int initialized;
    int threadprovided;
    int rank;

    ret = MPI_Initialized(&initialized);

    if (!initialized) {
        ret = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED,
                              &threadprovided);
    }

    return;
}

void cal::mpi_finalize() {
    std::cout << "Finalize MPI" << std::endl;
    int ret = MPI_Finalize();
    return;
}
