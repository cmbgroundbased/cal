
// Copyright (c) 2015-2020 by the parties listed in the AUTHORS file.
// All rights reserved.  Use of this source code is governed by
// a BSD-style license that can be found in the LICENSE file.

#include <_libcal.hpp>

using size_container = py::detail::any_container <ssize_t>;


PYBIND11_MODULE(_libcal, m) {
    m.doc() = R"(
    Interface to C++ CAL library.

    )";

    // Register aligned array types
    register_aligned <cal::AlignedI8> (m, "AlignedI8");
    register_aligned <cal::AlignedU8> (m, "AlignedU8");
    register_aligned <cal::AlignedI16> (m, "AlignedI16");
    register_aligned <cal::AlignedU16> (m, "AlignedU16");
    register_aligned <cal::AlignedI32> (m, "AlignedI32");
    register_aligned <cal::AlignedU32> (m, "AlignedU32");
    register_aligned <cal::AlignedI64> (m, "AlignedI64");
    register_aligned <cal::AlignedU64> (m, "AlignedU64");
    register_aligned <cal::AlignedF32> (m, "AlignedF32");
    register_aligned <cal::AlignedF64> (m, "AlignedF64");

    init_sys(m);
    init_math_sf(m);
    init_math_rng(m);
    init_math_healpix(m);
    init_atm(m);
    init_math_qarray(m);

}
