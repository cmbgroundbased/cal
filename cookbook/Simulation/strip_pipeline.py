# import system packages
import os
from datetime import datetime
import sys
import argparse
import traceback
import pickle
import yaml
from argparse import ArgumentParser, RawTextHelpFormatter

# import numpy
import numpy as np
import healpy as hp

# the MPI comunicator class, customized for pycal data storage
from pycal.mpi import get_world, Comm

# the Data class 
from pycal.dist import distribute_uniform, Data

# some pycal utils to share informations with the environment
from pycal.utils import Logger, Environment, memreport

# some timers
from pycal.timing import function_timer, GlobalTimers, Timer, gather_timers
from pycal.timing import dump as dump_timing

# the simulate the pointing, the atmosphere and put all the information in the TODGround class
from pycal.todmap import TODGround, OpSimAtmosphere, OpPointingHpix
from pycal.weather import Weather

# Some wrapper to libaatm, they solve the radiative transfer equation in local thermodynamic equilibrium
from pycal.todmap.atm import atm_absorption_coefficient, atm_absorption_coefficient_vec

# helper functions
from pycal.tests._helpers import boresight_focalplane
import pycal.qarray as qa

# focal plane and telescope calsses
from pycal.todmap import Focalplane
from pycal.todmap import Telescope

# set up the output directory for each mc iterations

@function_timer
def setup_output(outdir, comm, mc, freq):
    outpath = "{}/{:08}/{:03}".format(outdir, mc, int(freq))
    if comm.world_rank == 0:
        print("Creating the outpath: {}".format(outpath))
        os.makedirs(outpath, exist_ok=True)
    return outpath

def load_focalplane(args, comm):
    focalplane = None

    # Load focalplane information

    if comm.comm_world is None or comm.comm_world.rank == 0:
        if focalplane is None:
            detector_data = {}
            with open(r'./strip_focal_plane.yaml') as file:
                focalplane=yaml.full_load(file)

            detecotrs=focalplane['horns'].keys()
            for i in detecotrs:
                directions=focalplane['horns'][i]['orientation']
                l=np.arctan(directions[0]/directions[2])
                u=np.arctan(directions[1]/directions[2])
                zaxis = np.array([0, 0, 1], dtype=np.float64)
                
                angrot = qa.rotation(zaxis, 0 * np.pi / 180.0)
                wx = np.rad2deg(l) * np.pi / 180.0
                wy = np.rad2deg(u) * np.pi / 180.0
                wz = np.sqrt(1.0 - (wx * wx + wy * wy))
                wdir = np.array([wx, wy, wz])
                strip_quat = qa.from_vectors(zaxis, wdir)
                
                strip = {}
                strip["quat"] = strip_quat
                strip["fwhm"] = 20.0
                strip["fknee"] = 0.0
                strip["fmin"] = 1e-9
                strip["alpha"] = 1.0
                strip["NET"] = 1.0
                strip["color"] = "r"
            
                detector_data[i] = strip
            
            focalplane = Focalplane(
                detector_data=detector_data, sample_rate=args.sample_rate
            )
        else:
            focalplane = Focalplane(
                fname_pickle=args.focalplane, sample_rate=args.sample_rate
            )
    if comm.comm_world is not None:
        focalplane = comm.comm_world.bcast(focalplane, root=0)

    if args.debug:
        if comm.comm_world is None or comm.comm_world.rank == 0:
            outfile = "{}/focalplane.png".format(args.outdir)
            focalplane._plot_fp(12, 12, outfile)

    #schedule.telescope.focalplane = focalplane
    #detweights = focalplane.detweights

    return focalplane

def main():
    parser = ArgumentParser(formatter_class=RawTextHelpFormatter, description="Simulation arguments", fromfile_prefix_chars="@")
    
    # Important and required args:
    parser.add_argument("-ces_name",  required=True, type=str, help="name of the scanning (useful for sky patch)")
    parser.add_argument("-ces_start_time", default="2022,9,1,0,0,0", required=True)
    parser.add_argument("-ces_stop_time", default="2022,9,1,1,0,0", required=True)
    parser.add_argument("-sample_rate", required=True, type=int, help="Frequency sample rate [Hz]")

    # Scan parameters
    parser.add_argument("-ces_azmin", default=0, type=int)
    parser.add_argument("-ces_azmax", default=0, type=int)
    parser.add_argument("-ces_el", default=70.0, type=float)    
    parser.add_argument("-scan", default="spin", type=str, help="Type of scanning strategy")
    parser.add_argument("-subscan", default="spin_1hour", type=str, help="Type of scanning strategy")
    parser.add_argument("-scanrate", default=1.0, type=float)
    parser.add_argument("-scan_accel", default=1.0, type=float)

    
    # The focalplane. The default is None. In this case, it uses the Strip focal plane.
    parser.add_argument("-focalplane",  default=None, type=str, help="PKL file with the focalplane quats. If none, it uses the STRIP focalplane.")
    
    # Site parameters
    parser.add_argument("-site_name", default="Tenerife", type=str)
    parser.add_argument("-site_lon", default="-16:31:00", type=str)
    parser.add_argument("-site_lat", default="28:20:00", type=str)
    parser.add_argument("-site_alt", default=2390.0, type=float)
    parser.add_argument("-coord", default="C", type=str)

    # Map parameters
    parser.add_argument("-CES_start", default=None)
    parser.add_argument("-NSIDE", default=128, type=int)
    parser.add_argument("-debug", action='store_true')
    parser.add_argument("-outdir", required=True, default="out_directory", type=str)
    
    #Atmosphere arguments
    parser.add_argument("-start_mc", default=0, type=int)
    parser.add_argument("-nsimu", default=1, type=int)
    parser.add_argument("-cache_name", default="atm_")
    parser.add_argument("-atm_cache", default="atm_cache_")
    parser.add_argument("-verbose", default=0, type=int)
    parser.add_argument("-freq", default=43.0, type=float)
    
    
    # Arguments of the simulation
    # class args:
    #     sample_rate=20
    #     focalplane=None
    #     ces_name = "Test1"
    #     scan="spin"
    #     subscan="spin_1hour"
    #     ces_stop_time = datetime(2022, 9, 2, 0, 0, 0).timestamp()
    #     ces_start_time = datetime(2022, 9, 1, 0, 0, 0).timestamp()
    #     site_name= "Tenerife"
    #     site_lon = "-16:31:00"
    #     site_lat = "28:20:00"
    #     site_alt = 2390.0
    #     coord = "C"
    #     ces_azmin = 0
    #     ces_azmax = 0
    #     ces_el = 70
    #     scanrate = 1.0
    #     scan_accel = 1000
    #     CES_start = None
    #     NSIDE=128
    #     debug=True
    #     outdir="out_directory"
    
    # start_mc = 0
    # nsimu = 1
    # cache_name = "atm_"
    # atm_cache="atm_cache_"
    # verbose = 0
    # freq = 43 # GHz
    
    # definition of the logger, the global timer and the environment
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
    args = parser.parse_args()
    
    args.outdir = args.outdir+args.ces_start_time

    if comm.world_rank == 0:
        print("Creating the outdir: {}".format(args.outdir))
        os.makedirs(args.outdir, exist_ok=True)
    
    fp = load_focalplane(args, comm)

    # Create the TOD structure
    data = Data(comm)
    weather = "weather_STRIP.fits"
    
    sta = str(args.ces_start_time).split(",")
    sto = str(args.ces_stop_time).split(",")
    start_time = datetime(int(sta[0]), int(sta[1]), int(sta[2]), int(sta[3]), int(sta[4]), int(sta[5])).timestamp()
    stop_time = datetime(int(sto[0]), int(sto[1]), int(sto[2]), int(sto[3]), int(sto[4]), int(sto[5])).timestamp()
    
    totsamples = int((stop_time - start_time) * args.sample_rate)
    # create the TOD for this observation
    if comm.comm_group is not None:
        ndetrank = comm.comm_group.size
    else:
        ndetrank = 1

    try:
        tod = TODGround(
            comm.comm_group,
            fp.detquats,
            totsamples,
            detranks=ndetrank,
            firsttime=start_time,
            rate=args.sample_rate,
            site_lon=args.site_lon,
            site_lat=args.site_lat,
            site_alt=args.site_alt,
            azmin=args.ces_azmin,
            azmax=args.ces_azmax,
            el=args.ces_el,
            scanrate=args.scanrate,
            scan_accel=args.scan_accel,
            sinc_modulation=None,
            CES_start=None,
            CES_stop=None,
            sun_angle_min=None,
            coord=args.coord,
            sampsizes=None,
            report_timing=None,
            hwprpm=None,
            hwpstep=None,
            hwpsteptime=None,
        )
    except RuntimeError as e:
        raise RuntimeError(
            'Failed to create TOD for {}-{}-{}: "{}"'
            "".format(args.ces_name, args.scan, args.subscan, e)
        )

    # Create the observation, and append the tod
    obs = {}
    obs["name"] = "CES-{}-{}-{}-{}".format(
        args.site_name, args.ces_name, args.scan, args.subscan
    )
    obs["tod"] = tod
    obs["id"] = data.comm.group
    obs["telescope_id"] = 1
    obs["site"] = "Tenerife"
    obs["site_name"] = args.site_name
    obs["site_id"] = 123
    obs["altitude"] = args.site_alt
    obs["weather"] = Weather(weather, site=123)
    obs["fpradius"] = 10.0
    obs["start_time"] = start_time
    obs["focalplane"] = fp

    data.obs.append(obs)

    # Expand the pointing, interpolating the quaternions
    if comm.comm_world is not None:
        comm.comm_world.barrier()
    timer0.stop()
    if comm.world_rank == 0:
        timer0.report("Simulated scans")

    if comm.world_rank == 0:
        log.info("Expanding pointing")

    pointing = OpPointingHpix(
        nside=128,
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
    
    poin={}
    for i in obs['tod'].local_dets:
        p = obs['tod'].cache.reference("pixels_{}".format(i))
        poin[i]=p
    np.save(args.outdir+'/pointings', poin)
    
        
    # Atmospheric MC simulation 

    for mc in range(args.start_mc, args.start_mc + args.nsimu):
        
        timer_MC_iter = Timer()
        timer_MC_iter.start()

        log = Logger.get()
        tmr = Timer()
        tmr.start()
        if comm.world_rank == 0 and args.verbose:
            log.info("Simulating atmosphere")
            if args.atm_cache and not os.path.isdir(args.atm_cache):
                try:
                    os.makedirs(args.atm_cache)
                except FileExistsError:
                    pass

        common_atm_params = {
        "realization": mc,
        "component": 123456,
        "lmin_center": 0.01, # in m?
        "lmin_sigma": 0.001, 
        "lmax_center": 100,   # in m?
        "lmax_sigma": 10,
        "zatm": 40000.0,
        "zmax": 2000.0,
        "xstep": 100.0,
        "ystep": 100.0,
        "zstep": 100.0,
        "nelem_sim_max": 10000,
        "verbosity": 0,
        "gain": 1,
        "z0_center": 2000,
        "z0_sigma": 0,
        "apply_flags": False,
        "common_flag_name": None,
        "common_flag_mask": 255,
        "flag_name": None,
        "flag_mask": 255,
        "report_timing": True,
        "wind_dist": 10000,
        "flush": False,
        }

        # Simulate the atmosphere signal
        atm = OpSimAtmosphere(out="atm", cachedir=args.atm_cache, freq=args.freq, **common_atm_params)
        atm.exec(data)

        if comm.comm_world is not None:
            comm.comm_world.barrier()
        tmr.stop()
        if comm.world_rank == 0:
            tmr.report("Atmosphere simulation")


        if comm.world_rank == 0:
                log.info(
                    "Processing frequency {}GHz, MC = {}".format(args.freq, mc))
    
        # Set up the output directory
        mcoffset = args.freq * 1000000
        outpath = setup_output(args.outdir, comm, mc + mcoffset, args.freq)

        cache_name = "atm"
        log = Logger.get()
        if comm.world_rank == 0 and args.verbose:
            log.info("Scaling atmosphere by frequency")
        timer = Timer()
        timer.start()
        for obs in data.obs: # Now we have only one observation
            tod = obs["tod"]
            todcomm = tod.mpicomm

            weather = obs["weather"]
            focalplane = obs["focalplane"]

            start_time = obs["start_time"]
            weather.set(123, mc, start_time)
            altitude = obs["altitude"]
            air_temperature = weather.air_temperature
            surface_pressure = weather.surface_pressure
            pwv = weather.pwv
            # Use the entire processing group to sample the absorption
            # coefficient as a function of frequency
            freqmin = 0
            freqmax = 2 * args.freq
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
                    center = args.freq
                    width = 0.2 * args.freq
                nstep = 101
                # Interpolate the absorption coefficient to do a top hat
                # integral across the bandpass
                det_freqs = np.linspace(center - width / 2, center + width / 2, nstep)
                absorption_det = np.mean(np.interp(det_freqs, freqs, absorption))
                cachename = "{}_{}".format(cache_name, det)
                # print("{}_{}".format(cache_name, det))
                ref = tod.cache.reference(cachename)
                ref *= absorption_det
                del ref

        if comm.comm_world is not None:
            comm.comm_world.barrier()
        timer0.stop()
        if comm.world_rank == 0 and args.verbose:
            timer0.report("Atmosphere scaling")
        log = Logger.get()
        if comm.world_rank == 0 and args.verbose:
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
            args.freq,
        )
        obs["noise_scale"] = absorption * weather.air_temperature

        if comm.comm_world is not None:
            comm.comm_world.barrier()
        timer.stop()
        if comm.world_rank == 0 and args.verbose:
            timer.report("Atmosphere weighting")
            
        # Questa iterazione montecarlo puo` essere salvata in outhpath, no?
        tods = {}
        for i in obs['tod'].local_dets:
            t = obs['tod'].cache.reference("atm_{}".format(i))
            tods[i]=t    
        np.save(outpath+'/tod_mc_'+str(mc), tods)
        
        timer_MC_iter.stop()
        timer_MC_iter.report("Monte Carlo iteration completed in ")

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

    
