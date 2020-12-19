// Copyright (c) 2015-2020 by the parties listed in the AUTHORS file.
// All rights reserved.  Use of this source code is governed by
// a BSD-style license that can be found in the LICENSE file.

#ifndef CAL_MPI_TEST_TEST_HPP
#define CAL_MPI_TEST_TEST_HPP

#include <cal_mpi.hpp>
#include <gtest/gtest.h>


class MPICALShmemTest : public testing::Test {
    public:

        MPICALShmemTest() {}

        ~MPICALShmemTest() {}

        virtual void SetUp() {}

        virtual void TearDown() {}

        static const size_t n;
};


#endif // ifndef CAL_MPI_TEST_TEST_HPP
