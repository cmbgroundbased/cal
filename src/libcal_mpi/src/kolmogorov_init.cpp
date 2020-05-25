/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <CAL_MPI_AtmSim.hpp>
#include <fstream>
/**
* Numerically integrate the modified Kolmogorov correlation
* function at grid points. We integrate down from 10*kappamax
* to 0 for numerical precision.
*/
void cal::mpi_atm_sim::initialize_kolmogorov()
{
    MPI_Barrier(comm);
    double t1 = MPI_Wtime();

    rmin_kolmo = 0;
    double diag = sqrt(delta_x * delta_x + delta_y * delta_y);
    rmax_kolmo = sqrt(diag * diag + delta_z * delta_z) * 1.01;

    // Size of the interpolation grid;
    nr = 1000;

#ifdef DEBUG
    nr /= 10;
#endif // ifdef DEBUG

    rstep = (rmax_kolmo - rmin_kolmo) / (nr - 1);
    rstep_inv = 1. / rstep;

    kolmo_x.clear();
    kolmo_x.resize(nr, 0);
    kolmo_y.clear();
    kolmo_y.resize(nr, 0);

    double kappamin = 1. / lmax;
    double kappamax = 1. / lmin;
    double kappal = 0.9 * kappamax;
    double invkappal = 1. / kappal;
    double kappa0 = 0.75 * kappamin;
    double kappa0sq = kappa0 * kappa0;

    // Integration steps (has to be optimize!)
    long nkappa = 1000000;

    double upper_limit = 10 * kappamax;
    double kappastep = upper_limit / (nkappa - 1);
    double slope1 = 7. / 6.;
    double slope2 = -11. / 6.;

    if((rank==0) && (verbosity>0)){
        std::cerr << std::endl;
        std::cerr << "Evaluating Kolmogorov correlation at " << nr
                  << " different separations in range " << rmin_kolmo
                  << " - " << rmax_kolmo << " m" << std::endl;
        std::cerr << "kappamin = " << kappamin
                  << " 1/m, kappamax =  " << kappamax
                  << " 1/m. nkappa = " << nkappa << std::endl;
    }

    // Use the Newton's method to integrate the correlation function
    long nkappa_task = nkappa / ntask + 1;
    long first_kappa = nkappa_task * rank;
    long last_kappa = first_kappa + nkappa_task;
    if (last_kappa > nkappa) last_kappa = nkappa;

    // Precalculate the power spectrum function
    std::vector <double> phi(last_kappa - first_kappa);
    # pragma omp parallel for schedule(static, 10)
    for (long ikappa = first_kappa; ikappa < last_kappa; ++ikappa) {
        double kappa = ikappa * kappastep;
        double kkl = kappa * invkappal;
        phi[ikappa - first_kappa] =
            (1. + 1.802 * kkl - 0.254 * pow(kkl, slope1))
            * exp(-kkl * kkl) * pow(kappa * kappa + kappa0sq, slope2);
    }

    // Newton's method factors, not part of the power spectrum
    if (first_kappa == 0) phi[0] /= 2;
    if (last_kappa == nkappa) phi[last_kappa - first_kappa - 1] /= 2;

    // Integrate the power spectrum for a spherically symmetric
    // correlation function
    double nri = 1. / (nr - 1);
    double tau = 10.;
    double enorm = 1. / (exp(tau) - 1.);
    double ifac3 = 1. / (2. * 3.);

    # pragma omp parallel for schedule(static, 10)
    for (long ir = 0; ir < nr; ++ir) {
        double r = rmin_kolmo
                   + (exp(ir * nri * tau) - 1) * enorm * (rmax_kolmo - rmin_kolmo);
        double val = 0;
        if (r * kappamax < 1e-2) {
            // special limit r -> 0,
            // sin(kappa.r)/r -> kappa - kappa^3*r^2/3!
            for (long ikappa = first_kappa; ikappa < last_kappa; ++ikappa) {
                double kappa = ikappa * kappastep;
                double kappa2 = kappa * kappa;
                double kappa4 = kappa2 * kappa2;
                double r2 = r * r;
                val += phi[ikappa - first_kappa] * (kappa2 - r2 * kappa4 * ifac3);
            }
        } else {
            for (long ikappa = first_kappa; ikappa < last_kappa; ++ikappa) {
                double kappa = ikappa * kappastep;
                val += phi[ikappa - first_kappa] * sin(kappa * r) * kappa;
            }
            val /= r;
        }
        val *= kappastep;
        kolmo_x[ir] = r;
        kolmo_y[ir] = val;
    }

    if (MPI_Allreduce(MPI_IN_PLACE, kolmo_y.data(), (int)nr,
                      MPI_DOUBLE, MPI_SUM, comm))
        throw  std::runtime_error("Failed to allreduce kolmo_y");

    // Normalize
    double norm = 1. / kolmo_y[0];
    for (int i = 0; i < nr; ++i) kolmo_y[i] *= norm;

    if ((rank == 0) && (verbosity > 0)) {
        std::ofstream f;
        std::ostringstream fname;
        fname << "kolmogorov.txt";
        f.open(fname.str(), std::ios::out);
        for (int ir = 0; ir < nr;
             ir++) f << kolmo_x[ir] << " " << kolmo_y[ir] << std::endl;
        f.close();
    }

    // Measure the correlation length
    long icorr = nr - 1;
    while (fabs(kolmo_y[icorr]) < corrlim) --icorr;
    rcorr = kolmo_x[icorr];
    rcorrsq = rcorr * rcorr;

    double t2 = MPI_Wtime();

    if ((rank == 0) && (verbosity > 0)) {
        std::cerr << "rcorr = " << rcorr << " m (corrlim = "
                  << corrlim << ")" << std::endl;
        std::cerr << "Kolmogorov initialized in " << t2 - t1 << " s." << std::endl;
    }

    return;
}
