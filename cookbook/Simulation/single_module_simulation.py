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

from pycal.todmap import TODGround, OpSimAtmosphere
from pycal.weather import Weather

from pycal.tests._helpers import boresight_focalplane

@function_timer
def create_observation(args, comm, quat, ces, site, noise=1.0, verbose=True):
    """ Create an observation.

    Create an observation for the CES or SPIN scan 

    Args:
        args :  argparse arguments
        comm :  communicator
        ces (CES) :  One constant elevation scan

    """
    weather = site.weather
    totsamples = int((ces.stop_time - ces.start_time) * args.sample_rate)

    # create the TOD for this observation

    if comm.comm_group is not None:
        ndetrank = comm.comm_group.size
    else:
        ndetrank = 1

    try:
        tod = TODGround(
            comm.comm_group,
            quat,
            totsamples,
            detranks=ndetrank,
            firsttime=ces.start_time,
            rate=args.sample_rate,
            site_lon=site.lon,
            site_lat=site.lat,
            site_alt=site.alt,
            azmin=ces.azmin,
            azmax=ces.azmax,
            el=ces.el,
            scanrate=args.scan_rate,
            scan_accel=args.scan_accel,
            sinc_modulation=args.scan_sinc_modulate,
            CES_start=None,
            CES_stop=None,
            sun_angle_min=args.sun_angle_min,
            coord=args.coord,
            sampsizes=None,
            report_timing=args.debug,
            hwprpm=args.hwp_rpm,
            hwpstep=args.hwp_step_deg,
            hwpsteptime=args.hwp_step_time_s,
        )
    except RuntimeError as e:
        raise RuntimeError(
            'Failed to create TOD for {}-{}-{}: "{}"'
            "".format(ces.name, ces.scan, ces.subscan, e)
        )

    # Create the observation

    obs = {}
    obs["name"] = "CES-{}-{}-{}-{}".format(
        site.name, ces.name, ces.scan, ces.subscan
    )
    obs["tod"] = tod
    obs["baselines"] = None
    obs["noise"] = noise
    obs["id"] = int(ces.mjdstart * 10000)
    obs["intervals"] = tod.subscans
    obs["site"] = site
    obs["site_name"] = site.name
    obs["site_id"] = site.id
    obs["altitude"] = site.alt
    obs["weather"] = site.weather
    # obs["telescope"] = telescope
    # obs["telescope_name"] = telescope.name
    # obs["telescope_id"] = telescope.id
    # obs["focalplane"] = telescope.focalplane.detector_data
    # obs["fpradius"] = telescope.focalplane.radius
    obs["start_time"] = ces.start_time
    obs["season"] = ces.season
    obs["date"] = ces.start_date
    obs["MJD"] = ces.mjdstart
    obs["rising"] = ces.rising
    obs["mindist_sun"] = ces.mindist_sun
    obs["mindist_moon"] = ces.mindist_moon
    obs["el_sun"] = ces.el_sun
    return obs


@function_timer
def create_data_structure(args, comm, quat, ces):
    log = Logger.get()
    timer = Timer()
    timer.start()
    
    data = Data(comm)
    
    obs = create_observation(args, comm, quat, ces)
    
    if comm.comm_world is None or comm.comm_group.rank == 0:
        log.info("Group # {:4} has {} observations.".format(comm.group, len(data.obs)))
    
    if comm.comm_world is not None:
        comm.comm_world.barrier()
    if comm.world_rank == 0:
        timer.report("Simulation scans")

    return data

@function_timer
def setup_output(args, comm, mc, freq):
    outpath = "{}/{:08}/{:03}".format(args.outdir, mc, int(freq))
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
    site_lon = "-67:47:10"
    site_lat = "-22:57:30"
    site_alt = 5200.0
    coord = "C"
    ces_azmin = 45
    ces_azmax = 55
    ces_el = 70
    scanrate = 1.0
    scan_accel = 3.0
    CES_start = None
        
    
    # Load a focalplane
    # Single test detector
    
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

    # Create the observation

    obs = {}
    obs["name"] = "CES-{}-{}-{}-{}".format(
        site_name, ces_name, scan, subscan
    )
    obs["tod"] = tod
    obs["baselines"] = None
    obs["noise"] = noise
    obs["id"] = int(ces.mjdstart * 10000)
    obs["intervals"] = tod.subscans
    obs["site"] = site
    obs["site_name"] = site_name
    obs["site_id"] = 123
    obs["altitude"] = site_alt
    obs["weather"] = Weather(weather, site=123)
    obs["start_time"] = ces_start_time
    
   
    # I have to expand the pointing, interpolating the quaternions
    
    # ...
    
    start_mc = 0
    nsimu = 1000
    cache_name = ""
    verbose = 15
    freq = 43 # GHz
    
    for mc in range(start_mc, start_mc + nsimu):
    
        log = Logger.get()
        tmr = Timer()
        tmr.start()
        if comm.world_rank == 0 and verbose:
            log.info("Simulating atmosphere")
            if args.atm_cache and not os.path.isdir(args.atm_cache):
                try:
                    os.makedirs(args.atm_cache)
                except FileExistsError:
                    pass

        # Simulate the atmosphere signal
        atm = OpSimAtmosphere(
            out=cache_name,
            realization=mc,
            lmin_center=args.atm_lmin_center,
            lmin_sigma=args.atm_lmin_sigma,
            lmax_center=args.atm_lmax_center,
            gain=args.atm_gain,
            lmax_sigma=args.atm_lmax_sigma,
            zatm=args.atm_zatm,
            zmax=args.atm_zmax,
            xstep=args.atm_xstep,
            ystep=args.atm_ystep,
            zstep=args.atm_zstep,
            nelem_sim_max=args.atm_nelem_sim_max,
            verbosity=args.atm_verbosity,
            z0_center=args.atm_z0_center,
            z0_sigma=args.atm_z0_sigma,
            apply_flags=False,
            common_flag_mask=args.common_flag_mask,
            cachedir=args.atm_cache,
            flush=args.flush,
            wind_dist=args.atm_wind_dist,
        )
        atm.exec(data)

        if comm.comm_world is not None:
            comm.comm_world.barrier()
        tmr.stop()
        if comm.world_rank == 0:
            tmr.report("Atmosphere simulation")
        
        
        if comm.world_rank == 0:
                log.info(
                    "Processing frequency {}GHz, MC = {}".format(freq, mc))
        
        # Copy signal
        for obs in data.obs:
            tod = obs["tod"]
            for det in tod.local_dets:
                inref = tod.local_signal(det, input_dir)
                outname = "{}_{}".format(output_dir, det)
                outref = tod.cache.put(outname, inref, replace=force_cache)
                del outref
                del inref
        
        """Scale atmospheric fluctuations by frequency.

        Assume that cached signal under totalname_freq is pure atmosphere
        and scale the absorption coefficient according to the frequency.

        If the focalplane is included in the observation and defines
        bandpasses for the detectors, the scaling is computed for each
        detector separately and `freq` is ignored.

        """
        if not args.simulate_atmosphere:
            return

        log = Logger.get()
        if comm.world_rank == 0 and verbose:
            log.info("Scaling atmosphere by frequency")
        timer = Timer()
        timer.start()
        for obs in data.obs:
            tod = obs["tod"]
            todcomm = tod.mpicomm
            site_id = obs["site_id"]
            weather = obs["weather"]
            if "focalplane" in obs:
                focalplane = obs["focalplane"]
            else:
                focalplane = None
            start_time = obs["start_time"]
            weather.set(site_id, mc, start_time)
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
            if my_nfreq > 0:
                if atm_available_utils:
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
                else:
                    raise RuntimeError(
                        "Atmosphere utilities from libaatm are not available"
                    )
            else:
                my_freqs = np.array([])
                my_absorption = np.array([])
            if todcomm is None:
                freqs = my_freqs
                absorption = my_absorption
            else:
                freqs = np.hstack(todcomm.allgather(my_freqs))
                absorption = np.hstack(todcomm.allgather(my_absorption))
            # loading = atm_atmospheric_loading(altitude, pwv, freq)
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
        timer.stop()
        if comm.world_rank == 0 and verbose:
            timer.report("Atmosphere scaling")
            
        # Update weights
        """Update atmospheric noise weights.

        Estimate the atmospheric noise level from weather parameters and
        encode it as a noise_scale in the observation.  Madam will apply
        the noise_scale to the detector weights.  This approach assumes
        that the atmospheric noise dominates over detector noise.  To be
        more precise, we would have to add the squared noise weights but
        we do not have their relative calibration.

        """
        if not args.simulate_atmosphere:
            for obs in data.obs:
                obs["noise_scale"] = 1.0
            return
        log = Logger.get()
        if comm.world_rank == 0 and verbose:
            log.info("Updating atmospheric noise weights")
        timer = Timer()
        timer.start()
        if atm_available_utils:
            for obs in data.obs:
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
        else:
            raise RuntimeError("Atmosphere utilities from libaatm are not available")

        if comm.comm_world is not None:
            comm.comm_world.barrier()
        timer.stop()
        if comm.world_rank == 0 and verbose:
            timer.report("Atmosphere weighting")
        
        # Set up the output directory
        
        outpath = setup_output(args, comm, mc + mcoffset, freq)
        
        
    gt.stop_all()
    if mpiworld is not None:
        mpiworld.barrier()
    timer = Timer()
    timer.start()
    alltimers = gather_timers(comm=mpiworld)
    if comm.world_rank == 0:
        out = os.path.join(args.outdir, "timing")
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

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    