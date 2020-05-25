/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <CAL_MPI_AtmSim.hpp>
#include <fstream>

void cal::mpi_atm_sim::load_realization() {
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
        bool there = f.good();

        success = 0;
        if (there) {
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
            success = f.good();
            f.close();
        }

        if (success) {
            if (verbosity > 0) {
                std::cerr << "Loaded metada from "
                          << fname.str() << std::endl;
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
                std::cerr << " wdir = " << wdir * 180. / M_PI << " degrees" <<
                    std::endl;
                std::cerr << "   z0 = " << z0 << " m" << std::endl;
                std::cerr << "   T0 = " << T0 << " K" << std::endl;
                std::cerr << "rcorr = " << rcorr << " m (corrlim = "
                          << corrlim << ")" << std::endl;
            }
        } else {
            if (there) {
                std::cerr << "FAILED to load metadata from "
                          << fname.str() << std::endl;
            }
        }
    }

    // Send out the atmosphere configuration parameters
    // to all processes
    if (MPI_Bcast(&success, 1, MPI_CHAR, 0, comm))
        throw std::runtime_error("Failed to bcast success");

    if (!success) return;

    if (MPI_Bcast(&nn, 1, MPI_LONG, 0, comm))
        throw std::runtime_error("Failed to bcast nn");

    if (MPI_Bcast(&nelem, 1, MPI_LONG, 0, comm))
        throw std::runtime_error("Failed to bcast nn");

    if (MPI_Bcast(&nx, 1, MPI_LONG, 0, comm))
        throw std::runtime_error("Failed to bcast nx");

    if (MPI_Bcast(&ny, 1, MPI_LONG, 0, comm))
        throw std::runtime_error("Failed to bcast ny");

    if (MPI_Bcast(&nz, 1, MPI_LONG, 0, comm))
        throw std::runtime_error("Failed to bcast nz");

    if (MPI_Bcast(&delta_x, 1, MPI_DOUBLE, 0, comm))
        throw std::runtime_error("Failed to bcast delta_x");

    if (MPI_Bcast(&delta_y, 1, MPI_DOUBLE, 0, comm))
        throw std::runtime_error("Failed to bcast delta_y");

    if (MPI_Bcast(&delta_z, 1, MPI_DOUBLE, 0, comm))
        throw std::runtime_error("Failed to bcast delta_z");

    if (MPI_Bcast(&xstart, 1, MPI_DOUBLE, 0, comm))
        throw std::runtime_error("Failed to bcast xstart");

    if (MPI_Bcast(&ystart, 1, MPI_DOUBLE, 0, comm))
        throw std::runtime_error("Failed to bcast ystart");

    if (MPI_Bcast(&zstart, 1, MPI_DOUBLE, 0, comm))
        throw std::runtime_error("Failed to bcast zstart");

    if (MPI_Bcast(&maxdist, 1, MPI_DOUBLE, 0, comm))
        throw std::runtime_error("Failed to bcast maxdist");

    if (MPI_Bcast(&wx, 1, MPI_DOUBLE, 0, comm))
        throw std::runtime_error("Failed to bcast wx");

    if (MPI_Bcast(&wy, 1, MPI_DOUBLE, 0, comm))
        throw std::runtime_error("Failed to bcast wy");

    if (MPI_Bcast(&wz, 1, MPI_DOUBLE, 0, comm))
        throw std::runtime_error("Failed to bcast wz");

    if (MPI_Bcast(&lmin, 1, MPI_DOUBLE, 0, comm))
        throw std::runtime_error("Failed to bcast lmin");

    if (MPI_Bcast(&lmax, 1, MPI_DOUBLE, 0, comm))
        throw std::runtime_error("Failed to bcast lmax");

    if (MPI_Bcast(&w, 1, MPI_DOUBLE, 0, comm))
        throw std::runtime_error("Failed to bcast w");

    if (MPI_Bcast(&wdir, 1, MPI_DOUBLE, 0, comm))
        throw std::runtime_error("Failed to bcast wdir");

    if (MPI_Bcast(&z0, 1, MPI_DOUBLE, 0, comm))
        throw std::runtime_error("Failed to bcast z0");

    if (MPI_Bcast(&T0, 1, MPI_DOUBLE, 0, comm))
        throw std::runtime_error("Failed to bcast T0");

    zstride = 1;
    ystride = zstride * nz;
    xstride = ystride * ny;

    xstrideinv = 1. / xstride;
    ystrideinv = 1. / ystride;
    zstrideinv = 1. / zstride;

    // Load realization

    try {
        compressed_index = new mpi_shmem_long(nn, comm);
        compressed_index->set(-1);

        full_index = new mpi_shmem_long(nelem, comm);
        full_index->set(-1);
    } catch (...) {
        std::cerr << rank
                  << " : Failed to allocate element indices. nn = "
                  << nn << std::endl;
        throw;
    }
    try {
        realization = new mpi_shmem_double(nelem, comm);
        realization->set(0);
    } catch (...) {
        std::cerr << rank
                  << " : Failed to allocate realization. nelem = "
                  << nelem << std::endl;
        throw;
    }
    if (full_index->rank() == 0) {
        std::ostringstream fname_real;
        fname_real << cachedir << "/" << name.str() << "_realization.dat";
        std::ifstream freal(fname_real.str(),
                            std::ios::in | std::ios::binary);

        freal.read((char *)&(*full_index)[0],
                   full_index->size() * sizeof(long));
        success = freal.good();

        if (success) {
            for (int i = 0; i < nelem; ++i) {
                long ifull = (*full_index)[i];
                if ((ifull < 0) || (ifull > compressed_index->size() - 1)) {
                    // Cached file must be corrupt
                    success = false;
                    break;
                }
                (*compressed_index)[ifull] = i;
            }
            if (success) {
                freal.read((char *)&(*realization)[0],
                           realization->size() * sizeof(double));
                success = freal.good();
            }
        }
        freal.close();

        if (success) {
            if (verbosity > 0) std::cerr << "Loaded realization from "
                                         << fname_real.str() << std::endl;
        } else {
            std::cerr << "FAILED to load realization from "
                      << fname_real.str() << std::endl;
        }
    }

    if (MPI_Allreduce(MPI_IN_PLACE, &success, 1, MPI_CHAR, MPI_MIN, comm))
        throw std::runtime_error("Failed to allreduce success");

    if (!success) {
        delete compressed_index;
        delete full_index;
        delete realization;
    } else {
        cached = true;
    }

    return;
}


void cal::mpi_atm_sim::save_realization() {
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
