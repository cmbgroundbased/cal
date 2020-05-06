#include <CALAtmSim.hpp>
#include <sys_utils.hpp>
#include <sys_env.hpp>
// #inluce <qualcosa per PRNG>

#include <sstream>
#include <iostream>
#include <fstream>
#include <cstring>
#include <random>    // Ha un sacco di generatori
#include <functional>
#include <cmath>
#include <algorithm> // per fare il std::sort


double median(std::vector <double> vec) {
    if (vec.size() == 0) return 0;

    std::sort(vec.begin(), vec.end());
    int half1 = (vec.size() - 1) * .5;
    int half2 = vec.size() * .5;

    return .5 * (vec[half1] + vec[half2]);
}

double mean(std::vector <double> vec) {
    if (vec.size() == 0) return 0;

    double sum = 0;
    for (auto & val : vec) sum += val;

    return sum / vec.size();
}

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

    /* Costruiamo l'oggetto simulazione iniziando col definire il numero di processi e il numero di threads per processo, dedicati.*/
    std::cout << "CTOR atm_sim class" << std::endl;
    counter1 = counter1start;
    counter2 = counter2start;

    corrlim = 1e-3;

    ntask = 1;
    rank = 0;

    auto &  env = cal::Environment::get()
    nthread = env.max_threads();
    if ((rank == 0) && (verbosity > 0))
    {
        std:cerr << "atmsim constructed with " << ntask << " process, " << nthread << " threads per process." << std::endl;
    }

    /*Controlliamo i limiti entro i quali fare la scansione dell'atmosfera.*/

    if (azmin >= azmax) throw std::runtime_error("CTOR_atmsim ERROR: azmin >= azmax."); //Here we have to modify!

    if (elmin < 0) throw std::runtime_error("CTOR_atmsim ERROR: elmin < 0, you are going to observe the ground ... ");

    if (elmax > M_PI_2) throw std::runtime_error("CTOR_atmsim ERROR: elmax > pi/2.");

    if (elmin > elmax) throw std::runtime_error("CTOR_atmsim ERROR: elmin > elmax.");

    if (tmin > tmax) throw std::runtime_error("CTOR_atmsim ERROR: tmin > tmax.");

    if (lmin_center > lmax_center) throw std::runtime_error("CTOR_atmsim ERROR: lmin_center > lmax_center");

    // Griglia dell'atmosfera, inizializzazione di alcune variabili di supporto.

    xstepinv = 1 / xstep;
    ystepinv = 1 / ystep;
    zstepinv = 1 / zstep;

    // Span angolare che viene osservato

    delta_az = (azmax - azmin);
    delta_el = (elmax - elmin);
    delta_t = (tmax - tmin);

    // Starting point
    az0 = azmin + delta_az / 2;
    el0 = elmin + delta_el / 2;
    sinel0 = sin(el0);
    cosel0 = cos(el0);

    // Rotate the coordinate system. Align \hat{z} axis with \hat{r}. I don't undestand the utility. (I have to go deeper.)

    xxstep = xstep * cosel0 - zstep * sinel0;
    yystep = ystep;
    zzstep = xstep * sinel0 + zstep * cosel0;

    // Some prints

    /*bla bla bla */

    // Cholesky decomposition

    chcommon = &cholcommon;
    cholmod_start(chcommon);
    if (verbosity > 1){
        chcommon->print = 3;
    }
    else{
        chcommon->print = 1;
    }
    chcommon->itype = CHOLMOD_INT;
    chcommon->dtype = CHOLMOD_DOUBLE;
    // The factorization is LL' no LDL'
    chcommon->final_ll = 1;
}

cal::atm_sim::~atm_sim()
{
    // We don't neet a DTOR carefully definition. The unique_ptr(s) are authomatically free once they are derefenced.
    std::cout << "DTOR atm_sim class" << std::endl;
    compressed_index.reset();
    full_index.reset();
    realization.reset();
    cholmod_finish(chcommon);
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
