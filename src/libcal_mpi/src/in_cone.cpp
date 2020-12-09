/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <cal_mpi_internal.hpp>

/**
* Input coordinates are in the scan frame, rotate to
* horizontal frame
*/
bool cal::mpi_atm_sim::in_cone(double x, double y, double z, double t_in)
{

    double tstep = 1;

    for (double t = 0; t < delta_t; t += tstep) {
        if (t_in >= 0) {
            if (t != 0) break;
            t = t_in;
        }

        if ((t_in < 0) && (delta_t - t < tstep)) t = delta_t;

        double xtel_now = wx * t;
        double dx = x - xtel_now;

        // Is the point behind the telescope at this time?
        if (dx + xstep < 0) {
            if (t_in >= 0) std::cerr << "dx + xstep < 0: " << dx << std::endl;
            continue;
        }

        // Check the rest of the spherical coordinates
        double ytel_now = wy * t;
        double dy = y - ytel_now;

        double ztel_now = wz * t;
        double dz = z - ztel_now;

        double r = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (r > maxdist * 1.01) {
            if (t_in >= 0) std::cerr << "r > maxdist " << r << std::endl;
            continue;
        }

        if (dz > 0) {
            dz -= zstep;
        } else {
            dz += zstep;
        }

        if ((std::abs(dy) < 2 * ystep) && (std::abs(dz) < 2 * zstep)) return true;

        double dxx = dx * cosel0 - dz * sinel0;
        double dyy = dy;
        double dzz = dx * sinel0 + dz * cosel0;

        double el = std::asin(dzz / r);
        if ((el < elmin) || (el > elmax)) {
            if (t_in >= 0)
                std::cerr << "el outside cone: "
                          << el * 180 / M_PI << " not in "
                          << elmin * 180 / M_PI << " - "
                          << elmax * 180 / M_PI << std::endl;
            continue;
        }

        dxx = (dx + xstep) * cosel0 - dz * sinel0;
        double az = std::atan2(dyy, dxx);
        if (std::abs(az) > 0.5 * delta_az) {
            if (t_in >= 0)
                std::cerr << "abs(az) > delta_az/2 "
                          << az * 180 / M_PI << " > "
                          << 0.5 * delta_az * 180 / M_PI << std::endl;
            continue;
        }

        // Passed all the checks
        return true;
    }

    return false;
}
