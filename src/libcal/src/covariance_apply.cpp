/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <CALAtmSim.hpp>
#include <math_rng.hpp>
#include <fstream>

void cal::atm_sim::apply_sparse_covariance(cholmod_sparse * sqrt_cov,
                             long ind_start, long ind_stop)
{
    // Apply the Cholesky-decomposed (square-root) sparse covariance
    // matrix to a vector of Gaussian random numbers to impose the
    // desired correlation properties.

    cal::Timer tm;
    tm.start();

    size_t nelem = ind_stop - ind_start; // Number of elements in the slice

    // Draw the Gaussian variates in a single call

    cholmod_dense * noise_in = cholmod_allocate_dense(nelem, 1, nelem,
                                                      CHOLMOD_REAL, chcommon);
    cal::rng_dist_normal(nelem, key1, key2, counter1, counter2,
                           (double *)noise_in->x);

    cholmod_dense * noise_out = cholmod_allocate_dense(nelem, 1, nelem,
                                                       CHOLMOD_REAL, chcommon);

    // Apply the sqrt covariance to impose correlations

    int notranspose = 0;
    double one[2] = {1, 0};  // Complex one
    double zero[2] = {0, 0}; // Complex zero

    cholmod_sdmult(sqrt_cov, notranspose, one, zero, noise_in, noise_out, chcommon);
    if (chcommon->status != CHOLMOD_OK) throw std::runtime_error(
                  "cholmod_sdmult failed.");
    cholmod_free_dense(&noise_in, chcommon);

    // Subtract the mean of the slice to reduce step between the slices

    double * p = (double *)noise_out->x;
    double mean = 0, var = 0;
    for (long i = 0; i < nelem; ++i) {
        mean += p[i];
        var += p[i] * p[i];
    }
    mean /= nelem;
    var = var / nelem - mean * mean;
    for (long i = 0; i < nelem; ++i) p[i] -= mean;

    tm.stop();
    if (verbosity > 0) {
        std::ostringstream o;
        o << "Realization slice (" << ind_start << " -- "
          << ind_stop << ") var = " << var << ", constructed in";
        tm.report(o.str().c_str());
    }

    if (verbosity > 10) {
        std::ofstream f;
        std::ostringstream fname;
        fname << "realization_"
              << ind_start << "_" << ind_stop << ".txt";
        f.open(fname.str(), std::ios::out);
        for (long ielem = 0; ielem < nelem; ielem++) {
            double coord[3];
            ind2coord(ielem, coord);
            f << coord[0] << " " << coord[1] << " " << coord[2] << " "
              << p[ielem]  << std::endl;
        }
        f.close();
    }

    // Copy the slice realization over appropriate indices in
    // the full realization
    // FIXME: This is where we would blend slices

    for (long i = ind_start; i < ind_stop; ++i) {
        (*realization)[i] = p[i - ind_start];
    }

    cholmod_free_dense(&noise_out, chcommon);

    return;
}
