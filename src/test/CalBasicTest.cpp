#include <iostream>
#include <CALAtmSim.hpp>
using namespace std;

int main(){

    double azmin = 0.0;
    double azmax = M_PI/2;
    double elmin = M_PI/3;
    double elmax = M_PI/2.5;

    double tmin = 0;
    double tmax = 60*60*2;

    double lmin_center = 0.01;
    double lmin_sigma = 0.001;

    double lmax_center = 10;
    double lmax_sigma = 10;

    double w_center = 0;
    double w_sigma = 10;

    double wdir_center = 0;
    double wdir_sigma = 100;

    double z0_center = 2000;
    double z0_sigma = 0;

    double T0_center = 280;
    double T0_sigma = 10;

    double zatm = 40000;
    double zmax = 2000;

    double xstep = 100;
    double ystep = 100;
    double zstep = 100;

    long nelem_sim_max = 1000;

    int verbosity = 0;

    uint64_t key1 = 0;
    uint64_t key2 = 0;
    uint64_t counterval1 = 0;
    uint64_t counterval2 = 0;
    std::string cachedir = std::string("/tmp/");

    double rmin = 0;
    double rmax = 10000;


    cal::atm_sim *atm_strip = new cal::atm_sim(azmin, azmax, elmin, elmax, tmin, tmax, lmin_center, lmin_sigma, lmax_center, lmax_sigma, w_center, w_sigma, wdir_center, wdir_sigma, z0_center, z0_sigma, T0_center, T0_sigma, zatm, zmax, xstep, ystep, zstep, nelem_sim_max, verbosity, key1, key2, counterval1, counterval2, cachedir, rmin, rmax);




    atm_strip->simulate(true);
    // atm_strop->observe() prova


    return 0;
}
