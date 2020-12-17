/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <cal_mpi_internal.hpp>

/**
* For each sample, integrate alogn the line of sight by
* summing the atmosphere values. See Church (1995) Section
* 2.2 first equation. We omit the optical depth fatction
* which is close to unity.
*
* Coordinates at distance r. The scan is centered on the
* X-axis. Check if the top of the focal plane hits Z-MAX at
* this distance.  This way all lines-of-sight get integrated
* to the same distance
*
* The observe frame is co-moving with the wind.
*
* Combine atmospheric emission (via interpolation) with the
* ambient temperature.
* Note that the r^2 (beam area) and 1/r^2 (source distance)
* factors cancel in the integral.
*/
int cal::mpi_atm_sim::observe(double * t, double * az, double * el, double * tod,
            long nsamp, double fixed_r)
{
    if(!cached){
        throw std::runtime_error("There is no cached observation to observe.");
    }

    double t1 = MPI_Wtime();

    double zatm_inv = 1. / zatm;

    std::ostringstream o;
    o.precision(16);
    int error = 0;

    # pragma omp parallel for schedule(static, 100)
    for (long i = 0; i < nsamp; i++) {
        # pragma omp flush(error)
        if(error) continue;

        if ((!((azmin <= az[i]) && (az[i] <= azmax)) && ! ((azmin <= az[i] - 2 * M_PI) && (az[i] - 2 * M_PI <= azmax))) || !((elmin <= el[i]) && (el[i] <= elmax))) {
            o.precision(16);
            o << "atmsim::observe : observation out of bounds (az, el, t)"
              << " = (" << az[i] << ",  " << el[i] << ", " << t[i]
              << ") allowed: (" << azmin << " - " << azmax << ", "
              << elmin << " - " << elmax << ", "
              << tmin << " - " << tmax << ")"
              << std::endl;
            error = 1;
            # pragma omp flush(error)
            continue;
        }

        // the _now variable are refered
        // to this specific sample
        // The variable t[i] moves the telescope position.
        double t_now = t[i] - tmin;

        /** Relative to the center of field */
        double az_now = az[i] - az0;
        double el_now = el[i];

        double xtel_now = wx * t_now;
        double ytel_now = wy * t_now;
        double ztel_now = wz * t_now;

        double sin_el = sin(el_now);
        double sin_el_max = sin(elmax);
        double cos_el = cos(el_now);
        double sin_az = sin(az_now);
        double cos_az = cos(az_now);

        double r = 1.5 * xstep;
        double rstep = xstep;

        while (r < rmin) r += rstep;

        std::vector <long> last_ind(3);
        std::vector <double> last_nodes(8);

        double val = 0;
        if (fixed_r > 0) r = fixed_r;

        while (true){
            if (r > rmax) break;

            double zz = r * sin_el_max;
            if (zz >= zmax) break;

            // Horizontal coordinates
            zz = r * sin_el;
            double rproj = r * cos_el;
            double xx = rproj * cos_az;
            double yy = rproj * sin_az;

            // Rotate to scan frame
            double x = xx * cosel0 + zz * sinel0;
            double y = yy;
            double z = -xx * sinel0 + zz * cosel0;

            // Translate by the wind
            x += xtel_now;
            y += ytel_now;
            z += ztel_now;

# ifdef DEBUG
            if ((x < xstart) || (x > xstart + delta_x) ||
                (y < ystart) || (y > ystart + delta_y) ||
                (z < zstart) || (z > zstart + delta_z)) {
                std::cerr << "atmsim::observe : (x,y,z) out of bounds: "
                          << std::endl
                          << "x = " << x << std::endl
                          << "y = " << y << std::endl
                          << "z = " << z << std::endl;
                error = 1;
                # pragma omp flush (error)
                break;
            }
# endif // ifdef DEBUG

            // Atmospheric emission with ambient temperature
            double step_val;
            try {

                step_val = interp(x, y, z, last_ind, last_nodes) * (1. - z * zatm_inv);

            } catch (const std::runtime_error & e) {
                std::ostringstream o;
                o << "atmsim::observe : interp failed at " << std::endl
                  << "xxyyzz = (" << xx << ", " << yy << ", " << zz << ")"
                  << std::endl
                  << "xyz = (" << x << ", " << y << ", " << z << ")"
                  << std::endl
                  << "r = " << r << std::endl
                  << "tele at (" << xtel_now << ", " << ytel_now << ", "
                  << ztel_now << ")" << std::endl
                  << "( t, az, el ) = " << "( " << t[i] - tmin << ", "
                  << az_now * 180 / M_PI
                  << " deg , " << el_now * 180 / M_PI << " deg) "
                  << " in_cone(t) = " << in_cone(x, y, z, t_now)
                  << " with "
                  << std::endl << e.what() << std::endl;
                error = 1;
                # pragma omp flush(error)
                break;
            }
            val += step_val;
            // Prepare for the next step
            r += rstep;

            if (fixed_r > 0) break;
        }
        tod[i] = val * rstep * T0;
    }

    double t2 = MPI_Wtime();

    if ((rank == 0) && (verbosity > 0)) {
        if (fixed_r > 0){
            std::cerr << nsamp
                      << " samples observed at r =  "
                      << fixed_r << " in " << t2 - t1
                      << " sec." << std::endl;
        }
        else {
            std::cerr << nsamp << " samples observed in "
                      << t2 - t1 << " sec."
                      << std::endl;
        }
    }

    if (error) {
        std::cerr << "WARNING: atm::observe failed with: \""
                  << o.str() << "\"" << std::endl;
        return -1;
    }

    return 0;
}
