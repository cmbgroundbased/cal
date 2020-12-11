// Copyright (c) 2015-2020 by the parties listed in the AUTHORS file.
// All rights reserved.  Use of this source code is governed by
// a BSD-style license that can be found in the LICENSE file.

#include <cmath>

const size_t MPICALShmemTest::n = 10;

TEST_F(MPICALShmemTest, instantiate) {
    cal::mpi_shmem <double> shmem;
    shmem.allocate(n);

    cal::mpi_shmem <double> shmem2(n);

    size_t sz = shmem.size();

    EXPECT_EQ(shmem2.size(), n);

    EXPECT_EQ(shmem2.size(), shmem.size());
}

TEST_F(MPICALShmemTest, access) {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);

    cal::mpi_shmem <double> shmem(n);

    shmem.set(20);

    int ret = MPI_Barrier(comm);

    shmem.resize(2 * n);

    ret = MPI_Barrier(comm);

    double * p = shmem.data();

    if (rank == 0) {
        shmem[n - 1] = 10;
    }

    ret = MPI_Barrier(comm);

    EXPECT_FLOAT_EQ(shmem[n - 2], 20);
    EXPECT_FLOAT_EQ(shmem[n - 1], 10);
    EXPECT_EQ(shmem[n - 1], p[n - 1]);
}
