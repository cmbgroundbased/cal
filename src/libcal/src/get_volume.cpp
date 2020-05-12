#include <CALAtmSim.hpp>
#include <sys_utils.hpp>
#include <sys_env.hpp>
#include <math_rng.hpp>
// #inluce <qualcosa per PRNG>

#include <sstream>
#include <iostream>
#include <fstream>
#include <cstring>
#include <random>    // Ha un sacco di generatori
#include <functional>
#include <cmath>
#include <algorithm> // per fare il std::sort

void cal::atm_sim::get_volume()
{
    // Trim zmax if rmax sets a more stringent limit
    double zmax_from_rmax = rmax * sin(elmax);
    if (zmax > zmax_from_rmax) zmax = zmax_from_rmax;

    // Horizontal volume
    double delta_z_h = zmax;
    maxdist = delta_z_h / sinel0;

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

    delta_z_cone = z_max - z_min;
    delta_z = delta_z_cone;

    delta_x = maxdist;

    // The wind effect
    double wdx = std::abs(wx) * delta_t;
    double wdy = std::abs(wy) * delta_t;
    double wdz = std::abs(wz) * delta_t;

    delta_x += wdx;
    delta_y += wdy;
    delta_z += wdz;

    // Margin for interpolation

    delta_x += xstep;
    delta_y += 2 * ystep;
    delta_z += 2 * zstep;

    if (wx < 0) xstart = -wdx;
    else xstart = 0;

    if (wy < 0) ystart = -0.5 * delta_y_cone - wdy - ystep;
    else ystart = -0.5 * delta_y_cone - ystep;

    if (wy < 0) zstart = -0.5 * z_min - wdz - zstep;
    else zstart = -0.5 * z_min;

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

    initialize_kolmogorov();
}
