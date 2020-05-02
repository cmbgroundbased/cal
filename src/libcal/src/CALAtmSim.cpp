#include <CALAtmSim.hpp>
#include <sys_utils.hpp>
#include <sys_env.hpp>
#include <iostream>


// public methods

cal::atm_sim::atm_sim(double azmin, double azmax, double elmin, double elmax,
                        double tmin, double tmax,
                        double lmin_center, double lmin_sigma,
                        double lmax_center, double lmax_sigma,
                        double w_center, double w_sigma,
                        double wdir_center, double wdir_sigma,
                        double z0_center, double z0_sigma,
                        double T0_center, double T0_sigma,
                        double zatm, double zmax,
                        double xstep, double ystep, double zstep,
                        long nelem_sim_max,
                        int verbosity,
                        uint64_t key1, uint64_t key2,
                        uint64_t counterval1, uint64_t counterval2,
                        std::string cachedir,
                        double rmin, double rmax)
                        :
                        cachedir(cachedir),
                        verbosity(verbosity),
                        key1(key1), key2(key2),
                        counter1start(counterval1), counter2start(counterval2),
                        azmin(azmin), azmax(azmax),
                        elmin(elmin), elmax(elmax), tmin(tmin), tmax(tmax),
                        lmin_center(lmin_center), lmin_sigma(lmin_sigma),
                        lmax_center(lmax_center), lmax_sigma(lmax_sigma),
                        w_center(w_center), w_sigma(w_sigma),
                        wdir_center(wdir_center), wdir_sigma(wdir_sigma),
                        z0_center(z0_center), z0_sigma(z0_sigma),
                        T0_center(T0_center), T0_sigma(T0_sigma),
                        zatm(zatm), zmax(zmax),
                        xstep(xstep), ystep(ystep), zstep(zstep),
                        nelem_sim_max(nelem_sim_max),
                        rmin(rmin), rmax(rmax)
{
    std::cout << "CTOR atm_sim class" << std::endl;
}

cal::atm_sim::~atm_sim()
{
    std::cout << "DTOR atm_sim class" << std::endl;
}

int cal::atm_sim::simulate(bool use_cache)
{
    if (use_cache) load_realization();
    return 0;
}

int cal::atm_sim::observe(double * t, double * az, double * el, double * tod,
            long nsamp, double fixed_r)
{
    return 0;
}

void cal::atm_sim::print(std::ostream & out) const
{
    std::cout << "print" << std::endl;
}

// private methods


void cal::atm_sim::draw()
{
    std::cout << "draw function" << std::endl;
}

void cal::atm_sim::get_volume()
{
    std::cout << "get the atmopheric volume cone" << std::endl;
}

bool cal::atm_sim::in_cone(double x, double y, double z, double t_in){
    return true;
}
void cal::atm_sim::compress_volume()
{
    std::cout << "Reduct the volume" << std::endl;
}

void cal::atm_sim::get_slice(long & ind_start, long & ind_stop)
{
    std::cout << "Get the atmospheric slab" << std::endl;
}


cholmod_sparse * cal::atm_sim::sqrt_sparse_covariance(cholmod_sparse * cov,
                                        long ind_start, long ind_stop)
{
    return ?;
}
cholmod_sparse * cal::atm_sim::build_sparse_covariance(long ind_start, long ind_stop)
{
    return ?;
}

void apply_sparse_covariance(cholmod_sparse * cov,
                             long ind_start, long ind_stop)
{
    std::cout << "apply sparse covariance" << std::endl;
}

void ind2coord(long i, double * coord)
{
    std::cout << "ind2coord" << std::endl;
}

long coord2ind(double x, double y, double z)
{
    return 0;
}

double interp(double x, double y, double z, std::vector <long> & last_ind,
              std::vector <double> & last_nodes)
{
    return 0.0;
}

double cov_eval(double * coord1, double * coord2)
{
    return 1.0;
}

void initialize_kolmogorov()
{
    std::cout << "Start Kolmogorov" << std::endl;
}

double kolmogorov(double r)
{
    return 0;
}

void smooth()
{
    std::cout << "Smooth Kernel" << std::endl;
}

void load_realization()
{
    std::cout << "load realization" << std::endl;
}
void save_realization()
{
    std::cout << "save realization" << std::endl;
}
