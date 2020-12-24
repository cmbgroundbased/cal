"""
This script runs an atmosphere simulation
"""

import os
from datetime import datetime

if "STARTUP_DELAY_ATM_SIM" in os.environ:
    import numpy as np
    import time

    delay = np.float(os.environ["STARTUP_DELAY_ATM_SIM"])
    wait = np.random.rand() * delay
    time.sleep(wait)

import sys
import argparse
import traceback
import pickle

import numpy as np

from pycal.mpi import get_world, Comm

from pycal.dist import distribute_uniform, Data

from pycal.utils import Logger, Environment, memreport

from pycal.timing import function_timer, GlobalTimers, Timer, gather_timers
from pycal.timing import dump as dump_timing

from pycal.todmap import TODGround, OpSimAtmosphere, OpPointingHpix
from pycal.weather import Weather
from pycal.todmap.atm import atm_absorption_coefficient, atm_absorption_coefficient_vec

from pycal.tests._helpers import boresight_focalplane

@function_timer
def setup_output(outdir, comm, mc, freq):
    outpath = "{}/{:08}/{:03}".format(outdir, mc, int(freq))
    if comm.world_rank == 0:
        os.makedirs(outpath, exist_ok=True)
    return outpath

def main():
    log = Logger.get()
    gt = GlobalTimers.get()
    env = Environment.get()

    gt.start("Atmospheric simulation (globbal timer)")
    timer0 = Timer()
    timer0.start()
    
    # Get the communicator
    mpiworld, procs, rank = get_world()
    if rank == 0:
        print(env)
    if mpiworld is None:
        log.info("Running serially with one process at {}".format(str(datetime.now())))
    else:
        if rank == 0:
            log.info(
                "Running with {} processes at {}".format(procs, str(datetime.now()))
            )
    comm = Comm(world=mpiworld)
        
    # Args
    ces_name = "Test1"
    scan="spin"
    subscan="spin_1day"
    ces_stop_time = 100
    ces_start_time = 0
    sample_rate = 20
    site_name= "Tenerife"
    site_lon = "-16:31:00"
    site_lat = "28:20:00"
    site_alt = 2390.0
    coord = "C"
    ces_azmin = 45
    ces_azmax = 55
    ces_el = 70
    scanrate = 1.0
    scan_accel = 3.0
    CES_start = None
    NSIDE=256
    outdir="out_directory"
        
    
    # Load a focalplane, single test detector
    dnames, dquat, depsilon, drate, _, _, _, _ = boresight_focalplane(1, 20)


    # Create the TOD structure
    data = Data(comm)
    weather = "weather_strip.fits"
    totsamples = int((ces_stop_time - ces_start_time) * sample_rate)

    # create the TOD for this observation
    if comm.comm_group is not None:
        ndetrank = comm.comm_group.size
    else:
        ndetrank = 1

    try:
        tod = TODGround(
            comm.comm_group,
            dquat,
            totsamples,
            detranks=ndetrank,
            firsttime=ces_start_time,
            rate=sample_rate,
            site_lon=site_lon,
            site_lat=site_lat,
            site_alt=site_alt,
            azmin=ces_azmin,
            azmax=ces_azmax,
            el=ces_el,
            scanrate=scanrate,
            scan_accel=scan_accel,
            sinc_modulation=None,
            CES_start=None,
            CES_stop=None,
            sun_angle_min=None,
            coord=coord,
            sampsizes=None,
            report_timing=None,
            hwprpm=None,
            hwpstep=None,
            hwpsteptime=None,
        )
    except RuntimeError as e:
        raise RuntimeError(
            'Failed to create TOD for {}-{}-{}: "{}"'
            "".format(ces_name, scan, subscan, e)
        )

    # Create the observation, and append the tod
    obs = {}
    obs["name"] = "CES-{}-{}-{}-{}".format(
        site_name, ces_name, scan, subscan
    )
    obs["tod"] = tod
    obs["id"] = data.comm.group
    obs["telescope_id"] = 1
    obs["site"] = "Tenerife"
    obs["site_name"] = site_name
    obs["site_id"] = 123
    obs["altitude"] = site_alt
    obs["weather"] = Weather(weather, site=123)
    obs["fpradius"] = 10.0
    obs["start_time"] = ces_start_time
    
    data.obs.append(obs)
    
    if comm.comm_world is not None:
        comm.comm_world.barrier()
    timer0.stop()
    if comm.world_rank == 0:
        timer0.report("Simulated scans")
   
    # Expand the pointing, interpolating the quaternions
    if comm.world_rank == 0:
        log.info("Expanding pointing")

    pointing = OpPointingHpix(
        nside=NSIDE,
        nest=True,
        mode="IQU",
        single_precision=1e-7,
        nside_submap=128,
    )

    pointing.exec(data)

    if comm.comm_world is not None:
        comm.comm_world.barrier()
    if comm.world_rank == 0:
        timer0.report_clear("Pointing generation")

    # Atmospheric MC simulation 
    start_mc = 0
    nsimu = 1
    cache_name = ""
    atm_cache="atm_cache"
    verbose = 15
    freq = 43 # GHz
    
    for mc in range(start_mc, start_mc + nsimu):
    
        log = Logger.get()
        tmr = Timer()
        tmr.start()
        if comm.world_rank == 0 and verbose:
            log.info("Simulating atmosphere")
            if atm_cache and not os.path.isdir(atm_cache):
                try:
                    os.makedirs(atm_cache)
                except FileExistsError:
                    pass
                
        common_atm_params = {
        "realization": mc,
        "component": 123456,
        "lmin_center": 0.001,
        "lmin_sigma": 0.0001,
        "lmax_center": 1,
        "lmax_sigma": 0.1,
        "zatm": 40000.0,
        "zmax": 2000.0,
        "xstep": 100.0,
        "ystep": 100.0,
        "zstep": 100.0,
        "nelem_sim_max": 10000,
        "verbosity": 20,
        "gain": 1,
        "z0_center": 2000,
        "z0_sigma": 0,
        "apply_flags": True,
        "common_flag_name": None,
        "common_flag_mask": tod.TURNAROUND,
        "flag_name": None,
        "flag_mask": 255,
        "report_timing": True,
        "wind_dist": 10000,
        "flush": False,
        }
        
        # Simulate the atmosphere signal
        atm = OpSimAtmosphere(out="atm", cachedir=atm_cache, freq=freq, **common_atm_params)
        atm.exec(data)

        if comm.comm_world is not None:
            comm.comm_world.barrier()
        tmr.stop()
        if comm.world_rank == 0:
            tmr.report("Atmosphere simulation")
        
        
        if comm.world_rank == 0:
                log.info(
                    "Processing frequency {}GHz, MC = {}".format(freq, mc))
        
        """Scale atmospheric fluctuations by frequency.

        Assume that cached signal under totalname_freq is pure atmosphere
        and scale the absorption coefficient according to the frequency.

        If the focalplane is included in the observation and defines
        bandpasses for the detectors, the scaling is computed for each
        detector separately and `freq` is ignored.

        """
        log = Logger.get()
        if comm.world_rank == 0 and verbose:
            log.info("Scaling atmosphere by frequency")
        timer = Timer()
        timer.start()
        for obs in data.obs:
            tod = obs["tod"]
            todcomm = tod.mpicomm
            
            weather = obs["weather"]
            # focalplane = obs["focalplane"] # ????????????????
            
            start_time = obs["start_time"]
            weather.set(123, mc, start_time)
            altitude = obs["altitude"]
            air_temperature = weather.air_temperature
            surface_pressure = weather.surface_pressure
            pwv = weather.pwv
            # Use the entire processing group to sample the absorption
            # coefficient as a function of frequency
            freqmin = 0
            freqmax = 2 * freq
            nfreq = 1001
            freqstep = (freqmax - freqmin) / (nfreq - 1)
            if todcomm is None:
                nfreq_task = nfreq
                my_ifreq_min = 0
                my_ifreq_max = nfreq
            else:
                nfreq_task = int(nfreq // todcomm.size) + 1
                my_ifreq_min = nfreq_task * todcomm.rank
                my_ifreq_max = min(nfreq, nfreq_task * (todcomm.rank + 1))
            my_nfreq = my_ifreq_max - my_ifreq_min
            my_freqs = freqmin + np.arange(my_ifreq_min, my_ifreq_max) * freqstep
            my_absorption = atm_absorption_coefficient_vec(
                        altitude,
                        air_temperature,
                        surface_pressure,
                        pwv,
                        my_freqs[0],
                        my_freqs[-1],
                        my_nfreq,
                    )
   
   
            if todcomm is None:
                freqs = my_freqs
                absorption = my_absorption
            else:
                freqs = np.hstack(todcomm.allgather(my_freqs))
                absorption = np.hstack(todcomm.allgather(my_absorption))
            
            for det in tod.local_dets:
                try:
                    # Use detector bandpass from the focalplane
                    center = focalplane[det]["bandcenter_ghz"]
                    width = focalplane[det]["bandwidth_ghz"]
                except Exception:
                    # Use default values for the entire focalplane
                    center = freq
                    width = 0.2 * freq
                nstep = 101
                # Interpolate the absorption coefficient to do a top hat
                # integral across the bandpass
                det_freqs = np.linspace(center - width / 2, center + width / 2, nstep)
                absorption_det = np.mean(np.interp(det_freqs, freqs, absorption))
                cachename = "{}_{}".format(cache_name, det)
                ref = tod.cache.reference(cachename)
                ref *= absorption_det
                del ref

        if comm.comm_world is not None:
            comm.comm_world.barrier()
        timer0.stop()
        if comm.world_rank == 0 and verbose:
            timer0.report("Atmosphere scaling")
            
        # Update weights
        """Update atmospheric noise weights.

        Estimate the atmospheric noise level from weather parameters and
        encode it as a noise_scale in the observation.  Madam will apply
        the noise_scale to the detector weights.  This approach assumes
        that the atmospheric noise dominates over detector noise.  To be
        more precise, we would have to add the squared noise weights but
        we do not have their relative calibration.

        """

        log = Logger.get()
        if comm.world_rank == 0 and verbose:
            log.info("Updating atmospheric noise weights")
        timer = Timer()
        timer.start()


        site_id = obs["site_id"]
        weather = obs["weather"]
        start_time = obs["start_time"]
        weather.set(site_id, mc, start_time)
        altitude = obs["altitude"]
        absorption = atm_absorption_coefficient(
            altitude,
            weather.air_temperature,
            weather.surface_pressure,
            weather.pwv,
            freq,
        )
        obs["noise_scale"] = absorption * weather.air_temperature
        
        if comm.comm_world is not None:
            comm.comm_world.barrier()
        timer.stop()
        if comm.world_rank == 0 and verbose:
            timer.report("Atmosphere weighting")
        
        # Set up the output directory
        mcoffset = freq * 1000000
        outpath = setup_output(outdir, comm, mc + mcoffset, freq)
        
    gt.stop_all()
    if mpiworld is not None:
        mpiworld.barrier()  
    timer = Timer()
    timer.start()
    alltimers = gather_timers(comm=mpiworld)
    if comm.world_rank == 0:
        out = os.path.join(outdir, "timing")
        dump_timing(alltimers, out)
        timer.stop()
        timer.report("Gather and dump timing info")
        timer0.report_clear("single_module_simulation.py")

    return
        
    
if __name__ == "__main__":
    try:
        main()
    except Exception:
        # A sort of stack-trace...
        mpiworld, procs, rank = get_world()
        if procs == 1:
            raise
        exc_type, exc_value, exc_traceback = sys.exc_info()
        lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
        lines = ["Proc {}: {}".format(rank, x) for x in lines]
        print("".join(lines), flush=True)
        if mpiworld is not None:
            mpiworld.Abort(6)

    
