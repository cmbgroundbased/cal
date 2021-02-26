/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

/**
* \namespace atm
* \brief The atm namespace contain all the classes, methods and other application
* that are provided by libAATM
*/

#include <cal/sys_utils.hpp>
#include <cal/AATM_fun.hpp>

#ifdef HAVE_AATM

# include "ATMRefractiveIndexProfile.h"
# include "ATMPercent.h"
# include "ATMPressure.h"
# include "ATMNumberDensity.h"
# include "ATMMassDensity.h"
# include "ATMTemperature.h"
# include "ATMLength.h"
# include "ATMInverseLength.h"
# include "ATMOpacity.h"
# include "ATMHumidity.h"
# include "ATMFrequency.h"
# include "ATMWaterVaporRadiometer.h"
# include "ATMWVRMeasurement.h"
# include "ATMProfile.h"
# include "ATMSpectralGrid.h"
# include "ATMRefractiveIndex.h"
# include "ATMSkyStatus.h"
# include "ATMAngle.h"

/**
* Return the atmospheric profile
*/
atm::AtmProfile get_atmprofile(double altitude, double temperature,
                               double pressure) {
    /**Atmospheric type (to reproduce behavior above the tropopause)*/
    unsigned int atmType = 1;

    /**Ground level temperature [K]*/
    atm::Temperature T(temperature, "K");
    /**Altitude of the site [m]*/
    atm::Length Alt(altitude, "m");
    /**Water vapor scale height [km]*/
    atm::Length WVL(2, "km");
    /**The Troposheric lapse rate [K/km]*/
    double TLR = -5.6;
    /**Upper atmospherci boundary for calculations*/
    atm::Length topAtm(48.0, "km");
    /**Primary pressure step*/
    atm::Pressure Pstep(1.0, "mb");
    /**Pressure step ratio between two consecutive layer*/
    double PstepFact = 1.2;
    /**Pressure [Pa]*/
    atm::Pressure P(pressure, "Pa");
    /**Placeholder humidity (overridden by the PWV)*/
    atm::Humidity H(10, "%");


    return atm::AtmProfile(Alt, P, T, TLR, H, WVL, Pstep, PstepFact, topAtm, atmType);
}

/**
* Return the atmospheric load for a given location and observe frequency (monochromatic)
*
*/
atm::SkyStatus get_sky_status(double altitude, double temperature,
                              double pressure, double freq) {

    atm::AtmProfile atmo = get_atmprofile(altitude, temperature, pressure);
    atm::Frequency Freq(freq, "GHz");
    atm::RefractiveIndexProfile rip(Freq, atmo);

    return atm::SkyStatus(rip);
}

/**
* See get_sky_status. Return the frequency-band response.
*/
atm::SkyStatus get_sky_status_vec(double altitude, double temperature,
                                  double pressure,
                                  double freqmin, double freqmax,
                                  size_t nfreq) {
    /*
       Create an ATM SkyStatus object for the observing altitude and frequency.
     */
    atm::AtmProfile atmo = get_atmprofile(altitude, temperature, pressure);
    double freqstep = 0;
    if (nfreq > 1) freqstep = (freqmax - freqmin) / (nfreq - 1);

    // aatm SpectralGrid seems to have a bug.  The first grid point is
    // a whole grid step after the reference frequency.
    atm::SpectralGrid grid(nfreq, 0,
                           atm::Frequency(freqmin - freqstep, "GHz"),
                           atm::Frequency(freqstep, "GHz"));
    atm::RefractiveIndexProfile rip(grid, atmo);
    atm::SkyStatus ss(rip);
    return ss;
}

double cal::atm_get_absorption_coefficient(double altitude,
                                             double temperature,
                                             double pressure,
                                             double pwv,
                                             double freq) {

    atm::SkyStatus ss = get_sky_status(altitude, temperature, pressure, freq);
    ss.setUserWH2O(pwv, "mm");
    // double opacity = ss.getWetOpacity().get();
    double opacity = ss.getTotalOpacity().get(); // ???

    return 1 - exp(-opacity);
}

int cal::atm_get_absorption_coefficient_vec(double altitude,
                                              double temperature,
                                              double pressure,
                                              double pwv,
                                              double freqmin, double freqmax,
                                              size_t nfreq,
                                              double * absorption) {

    atm::SkyStatus ss = get_sky_status_vec(altitude, temperature, pressure,
                                           freqmin, freqmax, nfreq);
    ss.setUserWH2O(pwv, "mm");
    for (size_t i = 0; i < nfreq; ++i) {
        // double opacity = ss.getWetOpacity(i).get();
        double opacity = ss.getTotalOpacity(i).get(); // ???
      
        absorption[i] = 1 - exp(-opacity);
    }

    return 0;
}

double cal::atm_get_atmospheric_loading(double altitude,
                                          double temperature,
                                          double pressure,
                                          double pwv,
                                          double freq) {

    atm::SkyStatus ss = get_sky_status(altitude, temperature, pressure, freq);
    ss.setUserWH2O(pwv, "mm");

    return ss.getTebbSky().get();
}

int cal::atm_get_atmospheric_loading_vec(double altitude,
                                           double temperature,
                                           double pressure,
                                           double pwv,
                                           double freqmin, double freqmax,
                                           size_t nfreq, double * loading) {

    atm::SkyStatus ss = get_sky_status_vec(altitude, temperature, pressure,
                                           freqmin, freqmax, nfreq);
    ss.setUserWH2O(pwv, "mm");
    for (size_t i = 0; i < nfreq; ++i) {
        loading[i] = ss.getTebbSky(i).get();
    }

    return 0;
}

#else // ifdef HAVE_AATM

double cal::atm_get_absorption_coefficient(double altitude,
                                             double temperature,
                                             double pressure,
                                             double pwv,
                                             double freq) {
    auto here = cal_HERE();
    auto log = cal::Logger::get();
    std::string msg = "Atmosphere utilities require libaatm";
    log.error(msg.c_str(), here);
    throw std::runtime_error(msg.c_str());
    return 0.0;
}

int cal::atm_get_absorption_coefficient_vec(double altitude,
                                              double temperature,
                                              double pressure,
                                              double pwv,
                                              double freqmin, double freqmax,
                                              size_t nfreq,
                                              double * absorption) {
    auto here = cal_HERE();
    auto log = cal::Logger::get();
    std::string msg = "Atmosphere utilities require libaatm";
    log.error(msg.c_str(), here);
    throw std::runtime_error(msg.c_str());
    return 0;
}

double cal::atm_get_atmospheric_loading(double altitude,
                                          double temperature,
                                          double pressure,
                                          double pwv,
                                          double freq) {
    auto here = cal_HERE();
    auto log = cal::Logger::get();
    std::string msg = "Atmosphere utilities require libaatm";
    log.error(msg.c_str(), here);
    throw std::runtime_error(msg.c_str());
    return 0.0;
}

int cal::atm_get_atmospheric_loading_vec(double altitude,
                                           double temperature,
                                           double pressure,
                                           double pwv,
                                           double freqmin, double freqmax,
                                           size_t nfreq, double * loading) {
    auto here = cal_HERE();
    auto log = cal::Logger::get();
    std::string msg = "Atmosphere utilities require libaatm";
    log.error(msg.c_str(), here);
    throw std::runtime_error(msg.c_str());
    return 0;
}

#endif // ifdef HAVE_AATM
