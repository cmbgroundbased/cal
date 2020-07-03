/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <CAL_MPI_AtmSim.hpp>
#include <cstring>

/**
* Build a sparse covariance matrix first in the triplet form
* then cast it to the column-packed format.
*/
cholmod_sparse * cal::mpi_atm_sim::build_sparse_covariance(long ind_start, long ind_stop)
{

    double t1 = MPI_Wtime();

    std::vector <int> rows;
    std::vector <int> cols;
    std::vector <double> vals;

    // Number of elements in the slice
    size_t nelem = ind_stop - ind_start;
    std::vector <double> diagonal(nelem);

    // Fill the elements of the covariance matrix.
    # pragma omp parallel
    {
        std::vector <int> myrows, mycols;
        std::vector <double> myvals;

        # pragma omp for schedule(static, 10)
        for (int64_t i = 0; i < nelem; ++i) {
            double coord[3];
            ind2coord(i + ind_start, coord);
            diagonal[i] = cov_eval(coord, coord);
        }

        # pragma omp for schedule(static, 10)
        for (int64_t icol = 0; icol < nelem; ++icol) {
            // Translate indices into coordinates
            double colcoord[3];
            ind2coord(icol + ind_start, colcoord);
            for (int64_t irow = icol; irow < nelem; ++irow) {
                // Evaluate the covariance between the two coordinates
                double rowcoord[3];
                ind2coord(irow + ind_start, rowcoord);
                if (fabs(colcoord[0] - rowcoord[0]) > rcorr) continue;
                if (fabs(colcoord[1] - rowcoord[1]) > rcorr) continue;
                if (fabs(colcoord[2] - rowcoord[2]) > rcorr) continue;

                double val = cov_eval(colcoord, rowcoord);

                if(icol == irow){
                    // Regularize the matrix promoting the diagonal
                    val *= 1.01;
                }

                // If the covariance exceeds the threshold, add it to the
                // sparse matrix
                if (val * val > 1e-6 * diagonal[icol] * diagonal[irow]) {
                    myrows.push_back(irow);
                    mycols.push_back(icol);
                    myvals.push_back(val);
                }
            }
        }
        # pragma omp critical
        {
            rows.insert(rows.end(), myrows.begin(), myrows.end());
            cols.insert(cols.end(), mycols.begin(), mycols.end());
            vals.insert(vals.end(), myvals.begin(), myvals.end());
        }
    }

    double t2 = MPI_Wtime();

    if (verbosity > 0) {
        std::cerr << rank
                  << " : Sparse covariance evaluated in "
                  << t2 - t1 << " s." << std::endl;
    }

    // stype > 0 means that only the lower diagonal
    // elements of the symmetric matrix are needed.
    int stype = 1;
    size_t nnz = vals.size();

    cholmod_triplet * cov_triplet = cholmod_allocate_triplet(nelem,
                                                             nelem,
                                                             nnz,
                                                             stype,
                                                             CHOLMOD_REAL,
                                                             chcommon);
    memcpy(cov_triplet->i, rows.data(), nnz * sizeof(int));
    memcpy(cov_triplet->j, cols.data(), nnz * sizeof(int));
    memcpy(cov_triplet->x, vals.data(), nnz * sizeof(double));

    // Ensure vectors are freed
    std::vector <int>().swap(rows);
    std::vector <int>().swap(cols);
    std::vector <double>().swap(vals);
    cov_triplet->nnz = nnz;

    cholmod_sparse * cov_sparse = cholmod_triplet_to_sparse(cov_triplet,
                                                            nnz,
                                                            chcommon);
    if (chcommon->status != CHOLMOD_OK) {
         throw std::runtime_error("cholmod_triplet_to_sparse failed.");
     }

    cholmod_free_triplet(&cov_triplet, chcommon);

    t2 = MPI_Wtime();

    if (verbosity > 0) {
        std::cerr << rank << " : Sparse covariance constructed in "
                  << t2 - t1 << " s." << std::endl;
    }

    // Report memory usage

    double tot_mem = (nelem * sizeof(int) + nnz * (sizeof(int) + sizeof(double)))
                     / pow(2.0, 20.0);
    double max_mem = (nelem * nelem * sizeof(double)) / pow(2.0, 20.0);
    if (verbosity > 0) {
        std::cerr << std::endl;
        std::cerr << rank << " : Allocated " << tot_mem
                  << " MB for the sparse covariance matrix. "
                  << "Compression: " << tot_mem / max_mem << std::endl;
    }

    return cov_sparse;
}
