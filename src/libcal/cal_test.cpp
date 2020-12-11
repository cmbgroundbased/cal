
// Copyright (c) 2015-2020 by the parties listed in the AUTHORS file.
// All rights reserved.  Use of this source code is governed by
// a BSD-style license that can be found in the LICENSE file.

#include <cal.hpp>
#include <gtest/gtest.h>
#include <tests/cal_test.hpp>

#include <tests/cal_env_test.hpp>
#include <tests/cal_healpix_test.hpp>
#include <tests/cal_qarray_test.hpp>
#include <tests/cal_rng_test.hpp>
#include <tests/cal_sf_test.hpp>
#include <tests/cal_utils_test.hpp>

int main(int argc, char * argv[]) {

    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
      
}
