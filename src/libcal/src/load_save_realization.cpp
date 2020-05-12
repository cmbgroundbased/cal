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

void cal::atm_sim::load_realization()
{
    cached = false;

    std::ostringstream name;
    name << key1 << "_" << key2 << "_"
         << counter1start << "_" << counter2start;

    char success;

    if (rank == 0) {
        // Load metadata

        success = 1;

        std::ostringstream fname;
        fname << cachedir << "/" << name.str() << "_metadata.txt";

        std::ifstream f(fname.str());
        if (f.good()) {
            f >> nn;
            f >> nelem;
            f >> nx;
            f >> ny;
            f >> nz;
            f >> delta_x;
            f >> delta_y;
            f >> delta_z;
            f >> xstart;
            f >> ystart;
            f >> zstart;
            f >> maxdist;
            f >> wx;
            f >> wy;
            f >> wz;
            f >> lmin;
            f >> lmax;
            f >> w;
            f >> wdir;
            f >> z0;
            f >> T0;
            f.close();

            if (rank == 0 and verbosity > 0) {
                std::cerr << "Loaded metada from "
                          << fname.str() << std::endl;
            }
        } else success = 0;

        if ((verbosity > 0) && success) {
            std::cerr << std::endl;
            std::cerr << "Simulation volume:" << std::endl;
            std::cerr << "   delta_x = " << delta_x << " m" << std::endl;
            std::cerr << "   delta_y = " << delta_y << " m" << std::endl;
            std::cerr << "   delta_z = " << delta_z << " m" << std::endl;
            std::cerr << "    xstart = " << xstart << " m" << std::endl;
            std::cerr << "    ystart = " << ystart << " m" << std::endl;
            std::cerr << "    zstart = " << zstart << " m" << std::endl;
            std::cerr << "   maxdist = " << maxdist << " m" << std::endl;
            std::cerr << "        nx = " << nx << std::endl;
            std::cerr << "        ny = " << ny << std::endl;
            std::cerr << "        nz = " << nz << std::endl;
            std::cerr << "        nn = " << nn << std::endl;
            std::cerr << "Atmospheric realization parameters:" << std::endl;
            std::cerr << " lmin = " << lmin << " m" << std::endl;
            std::cerr << " lmax = " << lmax << " m" << std::endl;
            std::cerr << "    w = " << w << " m/s" << std::endl;
            std::cerr << "   wx = " << wx << " m/s" << std::endl;
            std::cerr << "   wy = " << wy << " m/s" << std::endl;
            std::cerr << "   wz = " << wz << " m/s" << std::endl;
            std::cerr << " wdir = " << wdir * 180. / M_PI << " degrees" << std::endl;
            std::cerr << "   z0 = " << z0 << " m" << std::endl;
            std::cerr << "   T0 = " << T0 << " K" << std::endl;
            std::cerr << "rcorr = " << rcorr << " m (corrlim = "
                      << corrlim << ")" << std::endl;
        }
    }

    if (!success) return;

    zstride = 1;
    ystride = zstride * nz;
    xstride = ystride * ny;

    xstrideinv = 1. / xstride;
    ystrideinv = 1. / ystride;
    zstrideinv = 1. / zstride;

    // Load realization

    try {
        compressed_index.reset(new AlignedVector <long> (nn));
        std::fill(compressed_index->begin(), compressed_index->end(), -1);

        full_index.reset(new AlignedVector <long> (nelem));
        std::fill(full_index->begin(), full_index->end(), -1);
    } catch (...) {
        std::cerr << rank
                  << " : Failed to allocate element indices. nn = "
                  << nn << std::endl;
        throw;
    }
    try {
        realization.reset(new AlignedVector <double> (nelem));
        std::fill(realization->begin(), realization->end(), 0.0);
    } catch (...) {
        std::cerr << rank
                  << " : Failed to allocate realization. nelem = "
                  << nelem << std::endl;
        throw;
    }
    // if (full_index->rank() == 0) {
    std::ostringstream fname_real;
    fname_real << cachedir << "/" << name.str() << "_realization.dat";
    std::ifstream freal(fname_real.str(),
                        std::ios::in | std::ios::binary);

    freal.read((char *)&(*full_index)[0],
               full_index->size() * sizeof(long));
    for (int i = 0; i < nelem; ++i) {
        long ifull = (*full_index)[i];
        (*compressed_index)[ifull] = i;
    }

    freal.read((char *)&(*realization)[0],
               realization->size() * sizeof(double));

    freal.close();

    if (verbosity > 0) std::cerr << "Loaded realization from "
                                 << fname_real.str() << std::endl;

    // }

    cached = true;

    return;
}
void cal::atm_sim::save_realization()
{
    if (rank == 0) {
        std::ostringstream name;
        name << key1 << "_" << key2 << "_"
             << counter1start << "_" << counter2start;

        // Save metadata

        std::ostringstream fname;
        fname << cachedir << "/" << name.str() << "_metadata.txt";

        std::ofstream f;
        f.precision(16);
        f.open(fname.str());
        f << nn << std::endl;
        f << nelem << std::endl;
        f << nx << std::endl;
        f << ny << std::endl;
        f << nz << std::endl;
        f << delta_x << std::endl;
        f << delta_y << std::endl;
        f << delta_z << std::endl;
        f << xstart << std::endl;
        f << ystart << std::endl;
        f << zstart << std::endl;
        f << maxdist << std::endl;
        f << wx << std::endl;
        f << wy << std::endl;
        f << wz << std::endl;
        f << lmin << std::endl;
        f << lmax << std::endl;
        f << w << std::endl;
        f << wdir << std::endl;
        f << z0 << std::endl;
        f << T0 << std::endl;
        f.close();

        if (verbosity > 0) std::cerr << "Saved metadata to "
                                     << fname.str() << std::endl;

        // Save realization

        std::ostringstream fname_real;
        fname_real << cachedir << "/" << name.str() << "_realization.dat";
        std::ofstream freal(fname_real.str(),
                            std::ios::out | std::ios::binary);

        freal.write((char *)&(*full_index)[0],
                    full_index->size() * sizeof(long));

        freal.write((char *)&(*realization)[0],
                    realization->size() * sizeof(double));

        freal.close();

        if (verbosity > 0) std::cerr << "Saved realization to "
                                     << fname_real.str() << std::endl;
    }

    return;
}
