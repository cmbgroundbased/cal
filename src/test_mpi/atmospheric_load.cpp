#include <iostream>
#include <fstream>
#include <CAL_MPI_AtmSim.hpp>
#include <AATM_fun.hpp>
#include <mpi_init.hpp>
#include <mpi.h>

int main(){
    int nsamp = 1001;
    double *abs_coef;
    double *bright_temp;

    abs_coef=(double *)malloc(nsamp*sizeof(double));
    bright_temp=(double *)malloc(nsamp*sizeof(double));

    double altitude = 2039.0; // m
    double pressure = 100000; // Pa?
    double temperature = 280.0; // K
    double pwv = 5.0; // mm
    double freqmin = 20.0; // GHz
    double freqmax = 70.0; // GHz


    cal::atm_get_absorption_coefficient_vec(altitude, temperature,
                                            pressure, pwv, freqmin,
                                            freqmax, nsamp, abs_coef);

    cal::atm_get_atmospheric_loading_vec(altitude, temperature,
                                         pressure, pwv, freqmin,
                                         freqmax, nsamp, bright_temp);
    std::ofstream fabs;
    std::ofstream fTbr;

    fabs.open("absorption.txt", std::ios::out);
    fTbr.open("bright_tem.txt", std::ios::out);

    for (int i = 0; i < nsamp; ++i) {
        fabs << abs_coef[i] << std::endl;
        fTbr << bright_temp[i] << std::endl;
    }
    fabs.close();
    fTbr.close();

    return 0;
}
