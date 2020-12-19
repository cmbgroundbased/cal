/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <cal_mpi_internal.hpp>

/**
* Cholesky-factorize the provided sparse matrix and return
* the sparse matrix representation of the factorization
*
* Extract band diagonal of the matrix and try
* factorizing again int ndiag = ntry - itry - 1;
*/
cholmod_sparse * cal::mpi_atm_sim::sqrt_sparse_covariance(cholmod_sparse * cov,
                                        long ind_start, long ind_stop)
{
    // Number of elements
    size_t nelem = ind_stop - ind_start;

    double t1 = MPI_Wtime();

    if (verbosity > 0) {
        std::cerr << rank
                  << " : Analyzing sparse covariance ... " << std::endl;
    }

    cholmod_factor * factorization;
    const int ntry = 4;
    for (int itry = 0; itry < ntry; ++itry) {
        factorization = cholmod_analyze(cov, chcommon);
        if (chcommon->status != CHOLMOD_OK) throw std::runtime_error(
                      "cholmod_analyze failed.");
        if (verbosity > 0) {
            std::cerr << rank
                      << " : Factorizing sparse covariance ... " << std::endl;
        }
        cholmod_factorize(cov, factorization, chcommon);
        if (chcommon->status != CHOLMOD_OK) {
            cholmod_free_factor(&factorization, chcommon);
            if (itry < ntry - 1) {
                // Extract band diagonal of the matrix and try
                // factorizing again
                // int ndiag = ntry - itry - 1;
                int ndiag = nelem - nelem * (itry + 1) / ntry;
                if (ndiag < 3) ndiag = 3;
                int iupper = ndiag - 1;
                int ilower = -iupper;
                if (verbosity > 0) {
                    cholmod_print_sparse(cov, "Covariance matrix", chcommon);

                    // DEBUG begin
                    if (itry > 2) {
                        FILE * covfile = fopen("failed_covmat.mtx", "w");
                        cholmod_write_sparse(covfile, cov, NULL, NULL, chcommon);
                        fclose(covfile);
                        exit(-1);
                    }

                    // DEBUG end
                    std::cerr << rank
                              << " : Factorization failed, trying a band "
                              << "diagonal matrix. ndiag = " << ndiag
                              << std::endl;
                }

                // Numerical (not pattern) matrix
                int mode = 1;
                cholmod_band_inplace(ilower, iupper, mode, cov, chcommon);
                if (chcommon->status != CHOLMOD_OK) throw std::runtime_error(
                              "cholmod_band_inplace failed.");
            } else throw std::runtime_error("cholmod_factorize failed.");
        } else {
            break;
        }
    }

    double t2 = MPI_Wtime();
    if (verbosity > 0) {
        std::cerr << std::endl;
        std::cerr << rank
                  << " : Cholesky decomposition done in " << t2 - t1
                  << " sec. N = " << nelem << std::endl;
    }

    // Report memory usage (only counting the non-zero elements, no
    // supernode information)
    size_t nnz = factorization->nzmax;
    double tot_mem = (nelem * sizeof(int) + nnz * (sizeof(int) + sizeof(double)))
                     / pow(2.0, 20.0);
    if (verbosity > 0) {
        std::cerr << std::endl;
        std::cerr << rank << " : Allocated " << tot_mem
                  << " MB for the sparse factorization." << std::endl;
    }

    cholmod_sparse * sqrt_cov = cholmod_factor_to_sparse(factorization, chcommon);
    if (chcommon->status != CHOLMOD_OK) throw std::runtime_error(
                  "cholmod_factor_to_sparse failed.");
    cholmod_free_factor(&factorization, chcommon);

    // Report memory usage
    nnz = sqrt_cov->nzmax;
    tot_mem = (nelem * sizeof(int) + nnz * (sizeof(int) + sizeof(double)))
              / pow(2.0, 20.0);
    double max_mem = (nelem * nelem * sizeof(double)) / pow(2.0, 20.0);
    if (verbosity > 0) {
        std::cerr << std::endl;
        std::cerr << rank << " : Allocated " << tot_mem
                  << " MB for the sparse sqrt covariance matrix. "
                  << "Compression: " << tot_mem / max_mem << std::endl;
    }

    return sqrt_cov;
}
