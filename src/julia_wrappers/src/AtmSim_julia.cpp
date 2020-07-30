#include "jlcxx/jlcxx.hpp"
#include <CALAtmSim.hpp>

using namespace cal;

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.add_type<atm_sim>("atm_sim")
    .constructor<double, // azmin
                 double, // azmax
                 double, // elmin
                 double, // elmax
                 double, // tmin
                 double, // tmax
                 double, double, // lmin_c, lmin_sigma
                 double, double, // lmax_c, lmax_sigma
                 double, double, // w_c, w_sigma
                 double, double, // wdir_c, wdir_sigma
                 double, double, // z0, z0_sigma
                 double, double, // T0_c, T0_sigma
                 double, // zatm
                 double, // zmax
                 double, double, double, // xstep, ystep, zstep
                 uint64_t, // nelem_sim_max
                 int, // verbosity
                 uint64_t, uint64_t, // key1, key1
                 uint64_t, uint64_t, // counterval1, counterval2
                 std::string, // chache dir
                 double, double // rmin, rmax
                 >()
    .method("observe", &atm_sim::observe)
    .method("simulate", &atm_sim::observe)
    .method("print", &atm_sim::print);
}
