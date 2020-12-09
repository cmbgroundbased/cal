
// Copyright (c) 2015-2020 by the parties listed in the AUTHORS file.
// All rights reserved.  Use of this source code is governed by
// a BSD-style license that can be found in the LICENSE file.

#ifndef CAL_MPI_SHM_HPP
#define CAL_MPI_SHM_HPP

#include <mpi.h>
#include <sstream>
#include <cal/sys_utils.hpp>


namespace cal {
/**
* Split the provided communicator into groups that share
* memory (are on the same node).
*/
template <typename T>
class mpi_shmem {
    public:

        mpi_shmem(MPI_Comm comm = MPI_COMM_WORLD) : comm_(comm) {
            int ret = MPI_Comm_rank(comm, &world_rank_);
            if (ret != MPI_SUCCESS) {
                auto here = cal_HERE();
                auto log = cal::Logger::get();
                std::string msg("Failed to get rank from input comm");
                log.error(msg.c_str(), here);
                throw std::runtime_error(msg.c_str());
            }

            ret = MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0,
                                      MPI_INFO_NULL, &shmcomm_);
            if (ret != MPI_SUCCESS) {
                if (world_rank_ == 0) {
                    auto here = cal_HERE();
                    auto log = cal::Logger::get();
                    std::string msg("Failed to split communicator by node.");
                    log.error(msg.c_str(), here);
                    throw std::runtime_error(msg.c_str());
                }
            }

            ret = MPI_Comm_size(shmcomm_, &ntasks_);
            if (ret != MPI_SUCCESS) {
                if (world_rank_ == 0) {
                    auto here = cal_HERE();
                    auto log = cal::Logger::get();
                    std::string msg("Failed to get node communicator size");
                    log.error(msg.c_str(), here);
                    throw std::runtime_error(msg.c_str());
                }
            }

            ret = MPI_Comm_rank(shmcomm_, &rank_);
            if (ret != MPI_SUCCESS) {
                if (world_rank_ == 0) {
                    auto here = cal_HERE();
                    auto log = cal::Logger::get();
                    std::string msg("Failed to get node communicator rank");
                    log.error(msg.c_str(), here);
                    throw std::runtime_error(msg.c_str());
                }
            }
        }

        mpi_shmem(int n, MPI_Comm comm = MPI_COMM_WORLD) : mpi_shmem(comm) {
            allocate(n);
        }

        T operator[](int i) const {
            return global_[i];
        }

        T & operator[](int i) {
            return global_[i];
        }

        /**
        * Determine the number of elements each process should offer
        * for the shared allocation
        */
        T * allocate(int n) {


            nlocal_ = n / ntasks_;

            if (nlocal_ * ntasks_ < n) nlocal_ += 1;
            if (nlocal_ * (rank_ + 1) > n) {
                nlocal_ = n - nlocal_ * rank_;
            }
            if (nlocal_ < 0) nlocal_ = 0;

            /** Allocate the shared memory */
            int ret = MPI_Win_allocate_shared(
                nlocal_ * sizeof(T), sizeof(T), MPI_INFO_NULL, shmcomm_,
                &local_, &win_);
            if (ret != MPI_SUCCESS) {
                if (rank_ == 0) {
                    auto here = cal_HERE();
                    auto log = cal::Logger::get();
                    std::ostringstream o;
                    o << " Failed to allocate " << n / 1024. / 1024.
                      << " MB of shared memory with " << ntasks_
                      << " tasks.";
                    log.error(o.str().c_str(), here);
                    throw std::runtime_error(o.str().c_str());
                }
            }
            n_ = n;

            MPI_Aint nn;
            int disp;

            /**
            * Get a pointer to the beginning of the shared memory
            * on rank # 0
            */
            ret = MPI_Win_shared_query(win_, 0, &nn, &disp, &global_);
            if (ret != MPI_SUCCESS) {
                if (rank_ == 0) {
                    auto here = cal_HERE();
                    auto log = cal::Logger::get();
                    std::string msg("Failed to query shared memory address.");
                    log.error(msg.c_str(), here);
                    throw std::runtime_error(msg.c_str());
                }
            }

            return global_;
        }

        void free() {
            if (global_) {
                int ret = MPI_Win_free(&win_);
                if (ret != MPI_SUCCESS) {
                    if (rank_ == 0) {
                        auto here = cal_HERE();
                        auto log = cal::Logger::get();
                        std::string msg("Failed to free shared memory.");
                        log.error(msg.c_str(), here);
                        throw std::runtime_error(msg.c_str());
                    }
                }
                local_ = NULL;
                global_ = NULL;
                n_ = 0;
                nlocal_ = 0;
            }
            return;
        }

        void set(T val) {
            for (int64_t i = 0; i < nlocal_; ++i) {
                local_[i] = val;
            }
            return;
        }

        /**
        * If there is memory already allocated, preserve its
        * contents.
        */
        T * resize(int64_t n) {

            int64_t n_copy = 0;
            cal::AlignedVector <T> temp;

            if (n < n_) {
                // We are shrinking the memory
                n_copy = n;
            } else if (n_ > 0) {
                // We are expanding the memory AND
                // we currently have a non-zero memory size
                n_copy = n_;
            }

            if ((n_copy > 0) && (rank_ == 0)) {
                temp.resize(n_copy);
                std::copy(global_, global_ + n_copy, temp.begin());
            }

            free();
            allocate(n);

            if ((n_copy > 0) && (rank_ == 0)) {
                std::copy(temp.begin(), temp.end(), global_);
            }

            return global_;
        }

        T * data() {
            return global_;
        }

        size_t size() {
            return n_;
        }

        int rank() {
            return rank_;
        }

        int ntasks() {
            return ntasks_;
        }

        ~mpi_shmem() {
            free();
            if (shmcomm_ != MPI_COMM_NULL) {
                MPI_Comm_free(&shmcomm_);
                shmcomm_ = MPI_COMM_NULL;
            }
        }

    private:

        T * local_ = NULL;
        T * global_ = NULL;
        int n_ = 0;
        int nlocal_ = 0;
        MPI_Comm comm_ = MPI_COMM_NULL;
        MPI_Comm shmcomm_ = MPI_COMM_NULL;
        MPI_Win win_ = MPI_WIN_NULL;
        int ntasks_;
        int rank_;
        int world_rank_;
};
}

#endif // ifndef CAL_MPI_SHM_HPP
