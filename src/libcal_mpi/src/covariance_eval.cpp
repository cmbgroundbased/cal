/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <cal_mpi_internal.hpp>

/**
* Evaluate the atmospheric absorption covariance between
* two coordinates Church (1995) Eq.(6) & (9).
* Coordinates are in the horizontal frame.
*/
double cal::mpi_atm_sim::cov_eval(double * coord1, double * coord2)
{
    const long nn = 1;
    const double ninv = 1.;

    double val = 0;

    for (int ii1 = 0; ii1 < nn; ++ii1) {
        double xx1 = coord1[0];
        double yy1 = coord1[1];
        double zz1 = coord1[2];

        for (int ii2 = 0; ii2 < nn; ++ii2) {
            double xx2 = coord2[0];
            double yy2 = coord2[1];
            double zz2 = coord2[2];

            double dx = xx1 - xx2;
            double dy = yy1 - yy2;
            double dz = zz1 - zz2;
            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < rcorrsq) {
                double r = sqrt(r2);

                // Water vapor altitude factor
                double chi1 = std::exp(-(zz1 + zz2) * z0inv);

                // Kolmogorov factor
                double chi2 = kolmogorov(r);

                val += chi1 * chi2;
            }
        }
    }

    return val * ninv;
}
