#include "jlcxx/jlcxx.hpp"
#include <AATM_fun.hpp>

using namespace cal;

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.method("atm_get_absorption_coefficient", &cal::atm_get_absorption_coefficient);
  mod.method("atm_get_absorption_coefficient_vec", &cal::atm_get_absorption_coefficient_vec);
  mod.method("atm_get_atmospheric_loading", &cal::atm_get_atmospheric_loading);
  mod.method("atm_get_atmospheric_loading_vec", &cal::atm_get_atmospheric_loading_vec);
}


