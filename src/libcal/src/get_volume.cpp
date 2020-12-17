/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <cal/CALAtmSim.hpp>

/**
* Trim zmax if rmax sets a more stringent limit
* Horizontal volume
* Cone witdth
*/

void cal::atm_sim::get_volume()
{
    // Trim zmax if rmax sets a more stringent limit
    double zmax_from_rmax = rmax * sin(elmax);
    if (zmax > zmax_from_rmax) zmax = zmax_from_rmax;

    // Horizontal volume
    double delta_z_h = zmax;

    // Maximum distance observed through the simulated volume
    maxdist = delta_z_h / sinel0;

    // Volume length
    double delta_x_h = maxdist * cos(elmin);

    double x, y, z, xx, zz, r, rproj, z_min, z_max;
    r = maxdist;

    z = r * sin(elmin);
    rproj = r * cos(elmin);
    x = rproj * cos(0);
    z_min = -x * sinel0 + z * cosel0;

    z = r * sin(elmax);
    rproj = r * cos(elmax);
    x = rproj * cos(delta_az / 2);
    z_max = -x * sinel0 + z * cosel0;

    // Cone witdth
    rproj = r * cos(elmin);
    if (delta_az > M_PI) delta_y_cone = 2 * rproj;
    else delta_y_cone = 2 * rproj * cos(0.5 * (M_PI - delta_az));

    // Cone height
    delta_z_cone = z_max - z_min;

    // Rotate to obsvation plane
    delta_z = delta_z_cone;
    delta_x = maxdist;

    // The wind effect
    double wdx = std::abs(wx) * delta_t;
    double wdy = std::abs(wy) * delta_t;
    double wdz = std::abs(wz) * delta_t;

    // We can move the frame
    delta_x += wdx;
    delta_y = delta_y_cone + wdy;
    delta_z += wdz;

    // Margin for interpolation
    delta_x += xstep;
    delta_y += 2 * ystep;
    delta_z += 2 * zstep;

    // Translate the volume to allow for wind. Telescope
    // sits at (0, 0, 0) at t=0
    if (wx < 0) xstart = -wdx;
    else xstart = 0;

    if (wy < 0) ystart = -0.5 * delta_y_cone - wdy - ystep;
    else ystart = -0.5 * delta_y_cone - ystep;

    if (wz < 0) zstart = z_min - wdz - zstep;
    else zstart = z_min - zstep;

    // Grid points
    nx = delta_x / xstep + 1;
    ny = delta_y / ystep + 1;
    nz = delta_z / zstep + 1;

    nn = nx * ny * nz;

    // 1D storage of the volume
    zstride = 1;
    ystride = zstride * nz;
    xstride = ystride * ny;

    xstrideinv = 1. / xstride;
    ystrideinv = 1. / ystride;
    zstrideinv = 1. / zstride;

    if ((rank == 0) && (verbosity > 0)) {
        std::cerr << std::endl;
        std::cerr << "Simulation volume:" << std::endl;
        std::cerr << "   delta_x = " << delta_x << " m" << std::endl;
        std::cerr << "   delta_y = " << delta_y << " m" << std::endl;
        std::cerr << "   delta_z = " << delta_z << " m" << std::endl;
        std::cerr << "Observation cone along the X-axis:" << std::endl;
        std::cerr << "   delta_y_cone = " << delta_y_cone << " m" << std::endl;
        std::cerr << "   delta_z_cone = " << delta_z_cone << " m" << std::endl;
        std::cerr << "    xstart = " << xstart << " m" << std::endl;
        std::cerr << "    ystart = " << ystart << " m" << std::endl;
        std::cerr << "    zstart = " << zstart << " m" << std::endl;
        std::cerr << "   maxdist = " << maxdist << " m" << std::endl;
        std::cerr << "        nx = " << nx << std::endl;
        std::cerr << "        ny = " << ny << std::endl;
        std::cerr << "        nz = " << nz << std::endl;
        std::cerr << "        nn = " << nn << std::endl;
    }

    initialize_kolmogorov();
}
