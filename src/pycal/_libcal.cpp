
// Copyright (c) 2015-2020 by the parties listed in the AUTHORS file.
// All rights reserved.  Use of this source code is governed by
// a BSD-style license that can be found in the LICENSE file.

#include <_libcal.hpp>

using size_container = py::detail::any_container <ssize_t>;


PYBIND11_MODULE(_libcal_mpi, m) {
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
    init_mpi_atm(m);

    // Internal unit test runner
    // m.def(
    //     "libcal_tests", [](py::list argv) {
    //         size_t narg = argv.size();
    //         std::vector <std::string> argbuffer;
    //         for (auto const & a : argv) {
    //             argbuffer.push_back(py::cast <std::string> (a));
    //         }
    //         char ** carg = (char **)std::malloc(narg * sizeof(char *));
    //         for (size_t i = 0; i < narg; ++i) {
    //             carg[i] = &(argbuffer[i][0]);
    //         }
    //         cal::test_runner(narg, carg);
    //         free(carg);
    //         return;
    //     }, py::arg(
    //         "argv"), R"(
    //     Run serial compiled tests from the internal libcal.
    //
    //     Args:
    //         argv (list):  The sys.argv or compatible list.
    //
    //     Returns:
    //         None
    //
    // )");
}
