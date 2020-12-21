
// Copyright (c) 2015-2020 by the parties listed in the AUTHORS file.
// All rights reserved.  Use of this source code is governed by
// a BSD-style license that can be found in the LICENSE file.

#ifndef CAL_TOD_POINTING_HPP
#define CAL_TOD_POINTING_HPP

#include <cal/math_healpix.hpp>
#include <string>

namespace cal {
void pointing_matrix_healpix(cal::HealpixPixels const & hpix,
                             bool nest, double eps, double cal,
                             std::string const & mode, size_t n,
                             double const * pdata, double const * hwpang,
                             uint8_t const * flags,
                             int64_t * pixels, double * weights);
}

#endif // ifndef CAL_TOD_POINTING_HPP