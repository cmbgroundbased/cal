/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#ifndef CAL_ATM_SIM_HPP
#define CAL_ATM_SIM_HPP

extern "C" {
#include <cholmod.h>
}
#include <sys_env.hpp>
#include <sys_utils.hpp>

/**
*@namespace cal
*@brief The API to create, simulate and observe an atmosphere are contained in the cal namespace.
*/
namespace cal{


using vec_long = std::unique_ptr <AlignedVector <long> >;
using vec_double = std::unique_ptr <AlignedVector <double> >;

/**
* \class atm_sim
* \brief Atmosphere creation, evolution and observation
*/
class atm_sim{
    public:

        /*!Helper name for shared pointer (Objects of shared_ptr types have the ability of taking ownership of a pointer and share that ownership: once they take ownership, the group of owners of a pointer become responsible for its deletion when the last one of them releases that ownership.) */
        typedef std::shared_ptr <atm_sim> pshr;

        /*!Helper name for unique pointer (These objects have the ability of taking ownership of a pointer: once they take ownership they manage the pointed object by becoming responsible for its deletion at some point. The unique_ptr objects automatically delete the object they manage (using a deleter) as soon as they themselves are destroyed, or as soon as their value changes either by an assignment operation or by an explicit call to unique_ptr::reset. )*/
        typedef std::unique_ptr <atm_sim> puniq;

        atm_sim(
            /**CES azimutal range*/
            double azmin,
            /**CES azimutal range*/
            double azmax,
            /**CES elevation range*/
            double elmin,
            /**CES elevation range*/
            double elmax,
            /**CES time range*/
            double tmin,
            /**CES time range*/
            double tmax,

            /** Dissipation scale of the Kolmogorov turbulence [m] with its sigma*/
            double lmin_center = .01, double lmin_sigma = .001,

            /** Injection scale of the Kolmogorov turbulence [m] with its sigma*/
            double lmax_center = 10, double lmax_sigma = 10,

            /** Wind speed [m/s] with its sigma*/
            double w_center = 25, double w_sigma = 10,

            /** Wind direction [radians] with its sigma */
            double wdir_center = 0, double wdir_sigma = 100,

            /** Water vapor distribution [m] with its sigma*/
            double z0_center = 2000, double z0_sigma = 0,

            /** Ground temperature [K] with its sigma*/
            double T0_center = 280, double T0_sigma = 10,

            /**Atmosphere extent for temperature profile [m]*/
            double zatm = 40000,

            /** Water vapor extent for integration [m]*/
            double zmax = 2000,

            /**Size of the volume elements [elements]*/
            double xstep = 100, double ystep = 100, double zstep = 100,

            /** Size of the simulation slices [elements]*/
            long nelem_sim_max = 1000,

            /**Verbosity [bool]*/
            int verbosity = 0,

            /**The RNG keys. In this constructur the atm_sim RNG has two different keys, key1 and key2*/
            uint64_t key1 = 0, uint64_t key2 = 0,

            /**The RNG counters: counterval1 and 2. Define the quantity of numbers that every thrads have to extract.*/
            uint64_t counterval1 = 0, uint64_t counterval2 = 0,

            /**The chacedir is used by the class in order to save atmospheric realizations in order to save time for multi-color simulation. The same atmosphere is just scaling for the extintion (provided by lib AATM)*/
            std::string cachedir = std::string(),

            /**Line-of-sight observing limit [m] (rmin and rmax by default 0, 10000)*/
            double rmin = 0, double rmax = 10000
        );

        ~atm_sim();

        /**Simulate the atmosphere time evolution*/
        int simulate(bool use_cache);

        /**Observe the atmosphere at that specific time along a specific line of sight*/
        int observe(double * t, double * az, double * el, double * tod,
                    long nsamp, double fixed_r = -1);

        /**Helper function for print*/
        void print(std::ostream & out = std::cout) const;

    private:

        std::string cachedir;
        int rank, ntask, nthread;
        int verbosity;
        uint64_t key1, key2, counter1, counter2, counter1start, counter2start;

        double azmin, azmax, elmin, elmax, tmin, tmax, sinel0, cosel0;

        /**Helper coordinate for the in-cone calculation*/
        double tanmin, tanmax;

        double lmin_center, lmin_sigma, lmax_center, lmax_sigma,
               w_center, w_sigma, wdir_center, wdir_sigma,
               z0_center, z0_sigma, T0_center, T0_sigma, z0inv;

        double az0, el0, delta_az, delta_el, delta_t;

        double zatm, zmax;

        double xstep, ystep, zstep, delta_x, delta_y, delta_z;

        double xstart, ystart, zstart, xxstep, yystep, zzstep;

        double delta_y_cone, delta_z_cone, maxdist;

        double xstepinv, ystepinv, zstepinv;

        long nx, ny, nz, nn, xstride, ystride, zstride;

        double xstrideinv, ystrideinv, zstrideinv;

        size_t nelem;

        bool cached = false;

        double lmin, lmax, w, wdir, z0, T0, wx, wy, wz;

        /**Number of steps in the Kolmogorov grid*/
        long nr;

        /**Size of the independent X-direction slices*/
        long nelem_sim_max;

        /**Kolmogorov correlation grid*/
        double rmin_kolmo, rmax_kolmo, rstep, rstep_inv;
        /**Correlation length*/
        double rcorr, rcorrsq, corrlim;
        /**Line-of-sight integration limits*/
        double rmin, rmax;

        /**Mapping between full volume and observation cone*/
        vec_long compressed_index;

        /**Inverse mapping between full volume and observation cone*/
        vec_long full_index;

        cholmod_common cholcommon;
        cholmod_common * chcommon;
        /**Draw values of lmin, lmax, w, wdir T0 (and optionally z0).*/
        void draw();
        /**Determine the rectangular volume needed*/
        void get_volume();

        /**determine of the given coordinates are within observed volume*/
        bool in_cone(double x, double y, double z, double t_in = -1);
        /** Find the volume elements really needed*/
        void compress_volume();

        vec_double realization;

        /**Find the next range of compressed indices to simulate*/
        void get_slice(long & ind_start, long & ind_stop);

        /**
        * Use the atmospheric parameters for volume element covariance
        * Cholesky decompose (square root) the covariance matrix
        */
        cholmod_sparse * sqrt_sparse_covariance(cholmod_sparse * cov,
                                                long ind_start, long ind_stop);
        cholmod_sparse * build_sparse_covariance(long ind_start, long ind_stop);

        /** Create a realization out of the square root covariance matrix */
        void apply_sparse_covariance(cholmod_sparse * sqrt_cov,
                                     long ind_start, long ind_stop);

        /** Compressed index to X Y Z - coordinates*/
        void ind2coord(long i, double * coord);

        /** X Y Z - coordinates to Compressed index*/
        long coord2ind(double x, double y, double z);

        /** Interpolate the realization value to given coordinates */
        double interp(double x, double y, double z, std::vector <long> & last_ind,
                      std::vector <double> & last_nodes);

        /** Evaluate the covariance matrix */
        double cov_eval(double * coord1, double * coord2);

        /** Integrate the correlation from the Kolmorov spectrum */
        void initialize_kolmogorov();

        /** Interpolate the correlation from precomputed grid*/
        double kolmogorov(double r);

        /** Smooths the realization*/
        void smooth();

        std::vector <double> kolmo_x;
        std::vector <double> kolmo_y;
        void load_realization();
        void save_realization();
};
}

#endif // End of  CAL_ATM_SIM Class
