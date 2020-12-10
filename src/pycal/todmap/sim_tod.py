# Copyright (c) 2015-2020 by the parties listed in the AUTHORS file.
# All rights reserved.  Use of this source code is governed by
# a BSD-style license that can be found in the LICENSE file.

import numpy as np

from scipy.constants import degree

import healpy as hp

try:
    import ephem
except:
    ephem = None

from .. import qarray as qa

from ..timing import function_timer, Timer

from ..tod import Interval, TOD

from ..healpix import ang2vec

from .pointing_math import quat_equ2ecl, quat_equ2gal, quat_ecl2gal


XAXIS, YAXIS, ZAXIS = np.eye(3)

# FIXME: Once the global package control interface is added,
# move this to that unified location.
tod_buffer_length = 1048576


@function_timer
def simulate_hwp(tod, hwprpm, hwpstep, hwpsteptime):
    """ Simulate and cache HWP angle in the TOD
    """
    if hwprpm is None and hwpsteptime is None and hwpstep is None:
        # No HWP
        return

    if hwprpm is not None and hwpsteptime is not None:
        raise RuntimeError("choose either continuously rotating or stepped HWP")

    if hwpstep is not None and hwpsteptime is None:
        raise RuntimeError("for a stepped HWP, you must specify the time between steps")

    # compute effective sample rate

    times = tod.local_times()
    dt = np.mean(times[1:-1] - times[0:-2])
    rate = 1.0 / dt
    del times

    if hwprpm is not None:
        # convert to radians / second
        tod._hwprate = hwprpm * 2.0 * np.pi / 60.0
    else:
        tod._hwprate = None

    if hwpstep is not None:
        # convert to radians and seconds
        tod._hwpstep = hwpstep * np.pi / 180.0
        tod._hwpsteptime = hwpsteptime * 60.0
    else:
        tod._hwpstep = None
        tod._hwpsteptime = None

    offset, nsamp = tod.local_samples
    hwp_angle = None
    if tod._hwprate is not None:
        # continuous HWP
        # HWP increment per sample is:
        # (hwprate / samplerate)
        hwpincr = tod._hwprate / rate
        startang = np.fmod(offset * hwpincr, 2 * np.pi)
        hwp_angle = hwpincr * np.arange(nsamp, dtype=np.float64)
        hwp_angle += startang
    elif tod._hwpstep is not None:
        # stepped HWP
        hwp_angle = np.ones(nsamp, dtype=np.float64)
        stepsamples = int(tod._hwpsteptime * rate)
        wholesteps = int(offset / stepsamples)
        remsamples = offset - wholesteps * stepsamples
        curang = np.fmod(wholesteps * tod._hwpstep, 2 * np.pi)
        curoff = 0
        fill = remsamples
        while curoff < nsamp:
            if curoff + fill > nsamp:
                fill = nsamp - curoff
            hwp_angle[curoff:fill] *= curang
            curang += tod._hwpstep
            curoff += fill
            fill = stepsamples

    if hwp_angle is not None:
        # Choose the HWP angle between [0, 2*pi)
        hwp_angle %= 2 * np.pi
        tod.cache.put(tod.HWP_ANGLE_NAME, hwp_angle)

    return

class TODGround(TOD):
    """Provide a simple generator of ground-based detector pointing.

    Detector focalplane offsets are specified as a dictionary of
    4-element ndarrays.  The boresight pointing is a generic
    2-angle model.

    Args:
        mpicomm (mpi4py.MPI.Comm): the MPI communicator over which the
            data is distributed (or None).
        detectors (dictionary): each key is the detector name, and each
            value is a quaternion tuple.
        samples (int):  The total number of samples.
        firsttime (float): starting time of data.
        rate (float): sample rate in Hz.
        site_lon (float/str): Observing site Earth longitude in radians
            or a pyEphem string.
        site_lat (float/str): Observing site Earth latitude in radians
            or a pyEphem string.
        site_alt (float/str): Observing site Earth altitude in meters.
        scanrate (float): Sky scanning rate in degrees / second.
        scan_accel (float): Sky scanning rate acceleration in
            degrees / second^2 for the turnarounds.
        sinc_modulation (bool): Modulate the scan rate according to
             1/sin(az) to achieve uniform integration depth.
        CES_start (float): Start time of the constant elevation scan
        CES_stop (float): Stop time of the constant elevation scan
        sun_angle_min (float): Minimum angular distance for the scan and
            the Sun [degrees].
        sampsizes (list):  Optional list of sample chunk sizes which
            cannot be split.
        sampbreaks (list):  Optional list of hard breaks in the sample
            distribution.
        coord (str):  Sky coordinate system.  One of
            C (Equatorial), E (Ecliptic) or G (Galactic)
        report_timing (bool):  Report the time spent simulating the scan
            and translating the pointing.
        hwprpm (float): If None, a constantly rotating HWP is not
            included. Otherwise it is the rate (in RPM) of constant
            rotation.
        hwpstep (float): If None, then a stepped HWP is not included.
            Otherwise, this is the step in degrees.
        hwpsteptime (float): The time in minutes between HWP steps.
        All other keyword arguments are passed to the parent constructor.

    """

    TURNAROUND = 1
    LEFTRIGHT_SCAN = 2
    RIGHTLEFT_SCAN = 4
    LEFTRIGHT_TURNAROUND = LEFTRIGHT_SCAN + TURNAROUND
    RIGHTLEFT_TURNAROUND = RIGHTLEFT_SCAN + TURNAROUND
    SUN_UP = 8
    SUN_CLOSE = 16

    @function_timer
    def __init__(
        self,
        mpicomm,
        detectors,
        samples,
        firsttime=0.0,
        rate=100.0,
        site_lon=0,
        site_lat=0,
        site_alt=0,
        azmin=0,
        azmax=0,
        el=0,
        scanrate=1,
        scan_accel=0.1,
        CES_start=None,
        CES_stop=None,
        el_min=0,
        sun_angle_min=90,
        sampsizes=None,
        sampbreaks=None,
        coord="C",
        report_timing=True,
        hwprpm=None,
        hwpstep=None,
        hwpsteptime=None,
        sinc_modulation=False,
        **kwargs
    ):
        if samples < 1:
            raise RuntimeError(
                "TODGround must be instantiated with a positive number of "
                "samples, not samples == {}".format(samples)
            )

        if ephem is None:
            raise RuntimeError("Cannot instantiate a TODGround object without pyephem.")

        if sampsizes is not None or sampbreaks is not None:
            raise RuntimeError(
                "TODGround will synthesize the sizes to match the subscans."
            )

        if CES_start is None:
            CES_start = firsttime
        elif firsttime < CES_start:
            raise RuntimeError(
                "TODGround: firsttime < CES_start: {} < {}"
                "".format(firsttime, CES_start)
            )
        lasttime = firsttime + samples / rate
        if CES_stop is None:
            CES_stop = lasttime
        elif lasttime > CES_stop:
            raise RuntimeError(
                "TODGround: lasttime > CES_stop: {} > {}" "".format(lasttime, CES_stop)
            )

        self._firsttime = firsttime
        self._lasttime = lasttime
        self._rate = rate
        self._site_lon = site_lon
        self._site_lat = site_lat
        self._site_alt = site_alt
        #if azmin == azmax:
        #    raise RuntimeError("TODGround requires non-empty azimuth range")
        self._azmin = azmin * degree
        self._azmax = azmax * degree
        if el < 1 or el > 89:
            raise RuntimeError("Impossible CES at {:.2f} degrees".format(el))
        self._el = el * degree
        self._scanrate = scanrate * degree
        self._scan_accel = scan_accel * degree
        self._CES_start = CES_start
        self._CES_stop = CES_stop
        self._el_min = el_min
        self._sun_angle_min = sun_angle_min
        if coord not in "CEG":
            raise RuntimeError("Unknown coordinate system: {}".format(coord))
        self._coord = coord
        self._report_timing = report_timing
        self._sinc_modulation = sinc_modulation

        self._observer = ephem.Observer()
        self._observer.lon = self._site_lon
        self._observer.lat = self._site_lat
        self._observer.elevation = self._site_alt  # In meters
        self._observer.epoch = ephem.J2000  # "2000"
        # self._observer.epoch = -9786 # EOD
        self._observer.compute_pressure()

        self._min_az = None
        self._max_az = None
        self._min_el = None
        self._min_el = None

        self._az = None
        self._commonflags = None
        self._boresight_azel = None
        self._boresight = None

        # Set the boresight pointing based on the given scan parameters

        tm = Timer()
        if self._report_timing:
            if mpicomm is not None:
                mpicomm.Barrier()
            tm.start()

        # sizes, starts = self.simulate_scan(samples)
        if azmin == azmax:
            self._azmin = 0
            self._azmax = 2 * np.pi
            sizes, starts = self.simulate_scan(samples)
        else:
            sizes, starts = self.simulate_scan(samples)

        if self._report_timing:
            if mpicomm is not None:
                mpicomm.Barrier()
            tm.stop()
            if (mpicomm is None) or (mpicomm.rank == 0):
                tm.report("TODGround: simulate scan")
            tm.clear()
            tm.start()

        # Create a list of subscans that excludes the turnarounds.
        # All processes in the group still have all samples.

        self.subscans = []
        self._subscan_min_length = 10  # in samples
        for istart, istop in zip(self._stable_starts, self._stable_stops):
            if istop - istart < self._subscan_min_length:
                self._commonflags[istart:istop] |= self.TURNAROUND
                continue
            start = self._firsttime + istart / self._rate
            stop = self._firsttime + istop / self._rate
            self.subscans.append(
                Interval(start=start, stop=stop, first=istart, last=istop - 1)
            )

        self._commonflags[istop:] |= self.TURNAROUND

        if np.sum((self._commonflags & self.TURNAROUND) == 0) == 0:
            raise RuntimeError(
                "The entire TOD is flagged as turnaround. Samplerate too low "
                "({} Hz) or scanrate too high ({} deg/s)?"
                "".format(rate, scanrate)
            )

        if self._report_timing:
            if mpicomm is not None:
                mpicomm.Barrier()
            tm.stop()
            if (mpicomm is None) or (mpicomm.rank == 0):
                tm.report("TODGround: list valid intervals")
            tm.clear()
            tm.start()

        self._fp = detectors
        self._detlist = sorted(list(self._fp.keys()))

        # call base class constructor to distribute data

        props = {
            "site_lon": site_lon,
            "site_lat": site_lat,
            "site_alt": site_alt,
            "azmin": azmin,
            "azmax": azmax,
            "el": el,
            "scanrate": scanrate,
            "scan_accel": scan_accel,
            "el_min": el_min,
            "sun_angle_min": sun_angle_min,
        }
        super().__init__(
            mpicomm,
            self._detlist,
            samples,
            sampsizes=[samples],
            sampbreaks=None,
            meta=props,
            **kwargs
        )

        self._AU = 149597870.7
        self._radperday = 0.01720209895
        self._radpersec = self._radperday / 86400.0
        self._radinc = self._radpersec / self._rate
        self._earthspeed = self._radpersec * self._AU

        if self._report_timing:
            if mpicomm is not None:
                mpicomm.Barrier()
            tm.stop()
            if (mpicomm is None) or (mpicomm.rank == 0):
                tm.report("TODGround: call base class constructor")
            tm.clear()
            tm.start()

        self.translate_pointing()

        self.crop_vectors()

        if self._report_timing:
            if mpicomm is not None:
                mpicomm.Barrier()
            tm.stop()
            if (mpicomm is None) or (mpicomm.rank == 0):
                tm.report("TODGround: translate scan pointing")

        # If HWP parameters are specified, simulate and cache HWP angle

        simulate_hwp(self, hwprpm, hwpstep, hwpsteptime)

        return

    @function_timer
    def __del__(self):
        try:
            del self._boresight_azel
        except:
            pass
        try:
            del self._boresight
        except:
            pass
        try:
            del self._az
        except:
            pass
        try:
            del self._commonflags
        except:
            pass
        try:
            self.cache.clear()
        except:
            pass

    def to_JD(self, t):
        """
        Convert CAL UTC time stamp to Julian date
        """
        return t / 86400.0 + 2440587.5

    def to_DJD(self, t):
        """
        Convert CAL UTC time stamp to Dublin Julian date used
        by pyEphem.
        """
        return self.to_JD(t) - 2415020

    @property
    def scan_range(self):
        """
        (tuple):  The extent of the boresight pointing as (min_az, max_az,
            min_el, max_el) in radians.  Includes turnarounds.
        """
        return self._min_az, self._max_az, self._min_el, self._max_el

    @function_timer
    def simulate_scan(self, samples):
        """ Simulate a constant elevation scan, either constant rate or
        1/sin(az)-modulated.

        """

        # Begin by simulating one full scan with turnarounds.
        # It will be used to interpolate the full CES.
        # `nstep` is the number of steps used for one sweep or one
        # turnaround.

        nstep = 10000

        azmin, azmax = [self._azmin, self._azmax]
        if self._sinc_modulation:
            # We always simulate a rising sinc scan and then
            # mirror it if necessary
            azmin %= np.pi
            azmax %= np.pi
            if azmin > azmax:
                raise RuntimeError(
                    "Cannot scan across zero meridian with sinc-modulated scan"
                )
        elif azmax < azmin:
            azmax += 2 * np.pi
        t = self._CES_start
        all_t = []
        all_az = []
        all_flags = []
        # translate scan rate from sky to mount coordinates
        base_rate = self._scanrate / np.cos(self._el)
        # scan acceleration is already in the mount coordinates
        scan_accel = self._scan_accel

        # left-to-right

        tvec = []
        azvec = []
        t0 = t
        if self._sinc_modulation:
            t1 = t0 + (np.cos(azmin) - np.cos(azmax)) / base_rate
            tvec = np.linspace(t0, t1, nstep, endpoint=True)
            azvec = np.arccos(np.cos(azmin) + base_rate * t0 - base_rate * tvec)
        else:
            # Constant scanning rate, only requires two data points
            t1 = t0 + (azmax - azmin) / base_rate
            tvec = np.array([t0, t1])
            azvec = np.array([azmin, azmax])
        all_t.append(np.array(tvec))
        all_az.append(np.array(azvec))
        all_flags.append(np.zeros(tvec.size, dtype=np.uint8) | self.LEFTRIGHT_SCAN)
        t = t1

        # turnaround

        t0 = t
        if self._sinc_modulation:
            dazdt = base_rate / np.abs(np.sin(azmax))
        else:
            dazdt = base_rate
        t1 = t0 + 2 * dazdt / scan_accel
        tvec = np.linspace(t0, t1, nstep, endpoint=True)[1:]
        azvec = azmax + (tvec - t0) * dazdt - 0.5 * scan_accel * (tvec - t0) ** 2
        all_t.append(tvec[:-1])
        all_az.append(azvec[:-1])
        all_flags.append(
            np.zeros(tvec.size - 1, dtype=np.uint8) | self.LEFTRIGHT_TURNAROUND
        )
        t = t1

        # right-to-left

        tvec = []
        azvec = []
        t0 = t
        if self._sinc_modulation:
            t1 = t0 + (np.cos(azmin) - np.cos(azmax)) / base_rate
            tvec = np.linspace(t0, t1, nstep, endpoint=True)
            azvec = np.arccos(np.cos(azmax) - base_rate * t0 + base_rate * tvec)
        else:
            # Constant scanning rate, only requires two data points
            t1 = t0 + (azmax - azmin) / base_rate
            tvec = np.array([t0, t1])
            azvec = np.array([azmax, azmin])
        all_t.append(np.array(tvec))
        all_az.append(np.array(azvec))
        all_flags.append(np.zeros(tvec.size, dtype=np.uint8) | self.RIGHTLEFT_SCAN)
        t = t1

        # turnaround

        t0 = t
        if self._sinc_modulation:
            dazdt = base_rate / np.abs(np.sin(azmin))
        else:
            dazdt = base_rate
        t1 = t0 + 2 * dazdt / scan_accel
        tvec = np.linspace(t0, t1, nstep, endpoint=True)[1:]
        azvec = azmin - (tvec - t0) * dazdt + 0.5 * scan_accel * (tvec - t0) ** 2
        all_t.append(tvec)
        all_az.append(azvec)
        all_flags.append(
            np.zeros(tvec.size, dtype=np.uint8) | self.RIGHTLEFT_TURNAROUND
        )

        # Concatenate

        tvec = np.hstack(all_t)
        azvec = np.hstack(all_az)
        flags = np.hstack(all_flags)

        # Now interpolate the simulated scan to timestamps

        times = self._CES_start + np.arange(samples) / self._rate
        tmin, tmax = tvec[0], tvec[-1]
        tdelta = tmax - tmin
        self._az = np.interp((times - tmin) % tdelta, tvec - tmin, azvec)
        if self._sinc_modulation and self._azmin > np.pi:
            # We always simulate a rising sinc scan and then
            # mirror it if necessary
            self._az += np.pi
        ind = np.searchsorted(tvec - tmin, (times - tmin) % tdelta)
        ind[ind == tvec.size] = tvec.size - 1
        self._commonflags = flags[ind]

        # Subscan start indices

        daz = np.diff(self._az)
        starts = np.hstack(
            [
                [0],
                np.argwhere(np.logical_and(daz[:-1] < 0, daz[1:] > 0)).ravel() + 2,
                [samples],
            ]
        )
        turnflags = self._commonflags & self.TURNAROUND
        self._stable_starts = (
            np.argwhere(np.logical_and(turnflags[:-1] != 0, turnflags[1:] == 0)).ravel()
            + 1
        )
        if turnflags[0] == 0:
            self._stable_starts = np.hstack([[0], self._stable_starts])
        self._stable_stops = (
            np.argwhere(np.logical_and(turnflags[:-1] == 0, turnflags[1:] != 0)).ravel()
            + 2
        )
        if turnflags[-1] == 0:
            self._stable_stops = np.hstack([self._stable_stops, [samples]])

        sizes = np.diff(starts)
        if np.sum(sizes) != samples:
            raise RuntimeError("Subscans do not match samples")

        # Store the scan range

        self._az %= 2 * np.pi
        self._min_az = np.amin(self._az)
        self._max_az = np.amax(self._az)
        self._min_el = self._el
        self._max_el = self._el

        return sizes, starts[:-1]

    @function_timer
    def translate_pointing(self):
        """Translate Az/El into bore sight quaternions

        Translate the azimuth and elevation into bore sight quaternions.

        """
        # At this point, all processes still have all of the scan
        nsamp = len(self._az)
        rank = 0
        ntask = 1
        if self._mpicomm is not None:
            rank = self._mpicomm.rank
            ntask = self._mpicomm.size
        nsamp_task = nsamp // ntask + 1
        my_start = rank * nsamp_task
        my_stop = min(my_start + nsamp_task, nsamp)
        my_nsamp = max(0, my_stop - my_start)
        my_ind = slice(my_start, my_stop)

        # Remember that the azimuth is measured clockwise and the
        # longitude counter-clockwise
        my_azelquats = qa.from_angles(
            np.pi / 2 - np.ones(my_nsamp) * self._el,
            -(self._az[my_ind]),
            np.zeros(my_nsamp),
            IAU=False,
        )
        azelquats = None
        if self._mpicomm is None:
            azelquats = my_azelquats
        else:
            azelquats = np.vstack(self._mpicomm.allgather(my_azelquats))
        self._boresight_azel = azelquats

        my_times = self.local_times()[my_ind]
        azel2radec_times, azel2radec_quats = self._get_azel2radec_quats()
        my_azel2radec_quats = qa.slerp(my_times, azel2radec_times, azel2radec_quats)
        my_quats = qa.mult(my_azel2radec_quats, my_azelquats)
        del my_azelquats

        quats = None
        if self._mpicomm is None:
            quats = my_quats
        else:
            quats = np.vstack(self._mpicomm.allgather(my_quats))
        self._boresight = quats
        del my_quats
        return

    @function_timer
    def crop_vectors(self):
        """Crop the TOD vectors.

        Crop the TOD vectors to match the sample range assigned to this task.

        """
        offset, n = self.local_samples
        ind = slice(offset, offset + n)

        self._az = self.cache.put("az", self._az[ind])
        self._commonflags = self.cache.put(
            "common_flags", self._commonflags[ind], replace=True
        )
        self._boresight_azel = self.cache.put(
            "boresight_azel", self._boresight_azel[ind]
        )
        self._boresight = self.cache.put("boresight_radec", self._boresight[ind])
        return

    @function_timer
    def _get_azel2radec_quats(self):
        """Construct a sparsely sampled vector of Az/El->Ra/Dec quaternions.

        The interpolation times must be tied to the total observation so
        that the results do not change when data is distributed in time
        domain.

        """
        # One control point at least every 10 minutes.  Overkill but
        # costs nothing.
        n = max(2, 1 + int((self._lasttime - self._firsttime) / 600))
        times = np.linspace(self._firsttime, self._lasttime, n)
        quats = np.zeros([n, 4])
        for i, t in enumerate(times):
            quats[i] = self._get_coord_quat(t)
            # Make sure we have a consistent branch in the quaternions.
            # Otherwise we'll get interpolation issues.
            if i > 0 and (
                np.sum(np.abs(quats[i - 1] + quats[i]))
                < np.sum(np.abs(quats[i - 1] - quats[i]))
            ):
                quats[i] *= -1
        quats = qa.norm(quats)
        return times, quats

    @function_timer
    def _get_coord_quat(self, t):
        """Get the Az/El -> Ra/Dec conversion quaternion for boresight.

        We will apply atmospheric refraction and stellar aberration in
        the detector frame.

        """
        self._observer.date = self.to_DJD(t)
        # Set pressure to zero to disable atmospheric refraction.
        pressure = self._observer.pressure
        self._observer.pressure = 0
        # Rotate the X, Y and Z axes from horizontal to equatorial frame.
        # Strictly speaking, two coordinate axes would suffice but the
        # math is cleaner with three axes.
        #
        # PyEphem measures the azimuth East (clockwise) from North.
        # The direction is standard but opposite to ISO spherical coordinates.
        try:
            xra, xdec = self._observer.radec_of(0, 0, fixed=False)
            yra, ydec = self._observer.radec_of(-np.pi / 2, 0, fixed=False)
            zra, zdec = self._observer.radec_of(0, np.pi / 2, fixed=False)
        except Exception as e:
            # Modified pyephem not available.
            # Translated pointing will include stellar aberration.
            xra, xdec = self._observer.radec_of(0, 0)
            yra, ydec = self._observer.radec_of(-np.pi / 2, 0)
            zra, zdec = self._observer.radec_of(0, np.pi / 2)
        self._observer.pressure = pressure
        xvec, yvec, zvec = ang2vec(
            np.pi / 2 - np.array([xdec, ydec, zdec]), np.array([xra, yra, zra])
        )
        # Orthonormalize for numerical stability
        xvec /= np.sqrt(np.dot(xvec, xvec))
        yvec -= np.dot(xvec, yvec) * xvec
        yvec /= np.sqrt(np.dot(yvec, yvec))
        zvec -= np.dot(xvec, zvec) * xvec + np.dot(yvec, zvec) * yvec
        zvec /= np.sqrt(np.dot(zvec, zvec))
        # Solve for the quaternions from the transformed axes.
        X = (xvec[1] + yvec[0]) / 4
        Y = (xvec[2] + zvec[0]) / 4
        Z = (yvec[2] + zvec[1]) / 4
        """
        if np.abs(X) < 1e-6 and np.abs(Y) < 1e-6:
            # Avoid dividing with small numbers
            c = .5 * np.sqrt(1 - xvec[0] + yvec[1] - zvec[2])
            d = np.sqrt(c**2 + .5 * (zvec[2] - yvec[1]))
            b = np.sqrt(.5 * (1 - zvec[2]) - c**2)
            a = np.sqrt(1 - b**2 - c**2 - d**2)
        else:
        """
        d = np.sqrt(np.abs(Y * Z / X))  # Choose positive root
        c = d * X / Y
        b = X / c
        a = (xvec[1] / 2 - b * c) / d
        # qarray has the scalar part as the last index
        quat = qa.norm(np.array([b, c, d, a]))
        """
        # DEBUG begin
        errors = np.array(
            [
                np.dot(qa.rotate(quat, [1, 0, 0]), xvec),
                np.dot(qa.rotate(quat, [0, 1, 0]), yvec),
                np.dot(qa.rotate(quat, [0, 0, 1]), zvec),
            ]
        )
        errors[errors > 1] = 1
        errors = np.degrees(np.arccos(errors))
        if np.any(errors > 1) or np.any(np.isnan(errors)):
            raise RuntimeError(
                "Quaternion is not right: ({}), ({} {} {})" "".format(errors, X, Y, Z)
            )
        # DEBUG end
        """
        return quat

    @function_timer
    def free_azel_quats(self):
        self._boresight_azel = None
        self.cache.destroy("boresight_azel")

    @function_timer
    def free_radec_quats(self):
        self._boresight = None
        self.cache.destroy("boresight_radec")

    @function_timer
    def radec2quat(self, ra, dec, pa):
        qR = qa.rotation(ZAXIS, ra + np.pi / 2)
        qD = qa.rotation(XAXIS, np.pi / 2 - dec)
        qP = qa.rotation(ZAXIS, pa)  # FIXME: double-check this
        q = qa.mult(qR, qa.mult(qD, qP))

        if self._coord != "C":
            # Add the coordinate system rotation
            if self._coord == "G":
                q = qa.mult(quat_equ2gal, q)
            elif self._coord == "E":
                q = qa.mult(quat_equ2ecl, q)
            else:
                raise RuntimeError("Unknown coordinate system: {}".format(self._coord))
        return q

    def detoffset(self):
        return {d: np.asarray(self._fp[d]) for d in self._detlist}

    def _get(self, detector, start, n):
        # This class just returns data streams of zeros
        return np.zeros(n, dtype=np.float64)

    def _put(self, detector, start, data):
        raise RuntimeError("cannot write data to simulated data streams")
        return

    def _get_flags(self, detector, start, n):
        return np.zeros(n, dtype=np.uint8)

    def _put_flags(self, detector, start, flags):
        raise RuntimeError("cannot write flags to simulated data streams")
        return

    def _get_common_flags(self, start, n):
        return self._commonflags[start : start + n]

    def _put_common_flags(self, start, flags):
        raise RuntimeError("cannot write flags to simulated data streams")
        return

    def _get_hwp_angle(self, start, n):
        if self.cache.exists(self.HWP_ANGLE_NAME):
            angle = self.cache.reference(self.HWP_ANGLE_NAME)[start : start + n]
        else:
            angle = None
        return angle

    def _put_hwp_angle(self, start, flags):
        raise RuntimeError("cannot write HWP angle to simulated data streams")
        return

    def _get_times(self, start, n):
        start_abs = self.local_samples[0] + start
        start_time = self._firsttime + float(start_abs) / self._rate
        return start_time + np.arange(n) / self._rate

    def _put_times(self, start, stamps):
        raise RuntimeError("cannot write timestamps to simulated data streams")
        return

    def _get_boresight(self, start, n, azel=False):
        if azel:
            if self._boresight_azel is None:
                raise RuntimeError("Boresight azel pointing was purged.")
            return self._boresight_azel[start : start + n]
        else:
            if self._boresight is None:
                raise RuntimeError("Boresight radec pointing was purged.")
            return self._boresight[start : start + n]

    def _get_boresight_azel(self, start, n):
        return self._get_boresight(start, n, azel=True)

    def _put_boresight(self, start, data):
        raise RuntimeError("cannot write boresight to simulated data streams")
        return

    def _put_boresight_azel(self, start, data):
        raise RuntimeError("cannot write boresight to simulated data streams")
        return

    @function_timer
    def read_boresight_az(self, local_start=0, n=0):
        """Read the boresight azimuth.

        Args:
            local_start (int): the sample offset relative to the first locally
                assigned sample.
            n (int): the number of samples to read.  If zero, read to end.

        Returns:
            (array): a numpy array containing the timestamps.
        """
        if n == 0:
            n = self.local_samples[1] - local_start
        if self.local_samples[1] <= 0:
            raise RuntimeError(
                "cannot read boresight azimuth - process "
                "has no assigned local samples"
            )
        if (local_start < 0) or (local_start + n > self.local_samples[1]):
            raise ValueError(
                "local sample range {} - {} is invalid".format(
                    local_start, local_start + n - 1
                )
            )
        return self._az[local_start : local_start + n]

    @function_timer
    def _get_pntg(self, detector, start, n, azel=False):
        # FIXME: this is where we will apply atmospheric refraction and
        # stellar aberration corrections in the detector frame.  For
        # simulations they will only matter if we want to simulate the
        # error coming from ignoring them.
        boresight = self._get_boresight(start, n, azel=azel)
        detquat = self._fp[detector]
        return qa.mult(boresight, detquat)

    def _put_pntg(self, detector, start, data):
        raise RuntimeError("cannot write data to simulated pointing")
        return

    @function_timer
    def _get_position(self, start, n):
        # For this simple class, assume that the Earth is located
        # along the X axis at time == 0.0s.  We also just use the
        # mean values for distance and angular speed.  Classes for
        # real experiments should obviously use ephemeris data.
        rad = np.fmod((start - self._firsttime) * self._radpersec, 2.0 * np.pi)
        ang = self._radinc * np.arange(n, dtype=np.float64) + rad
        x = self._AU * np.cos(ang)
        y = self._AU * np.sin(ang)
        z = np.zeros_like(x)
        return np.ravel(np.column_stack((x, y, z))).reshape((-1, 3))

    def _put_position(self, start, pos):
        raise RuntimeError("cannot write data to simulated position")
        return

    @function_timer
    def _get_velocity(self, start, n):
        # For this simple class, assume that the Earth is located
        # along the X axis at time == 0.0s.  We also just use the
        # mean values for distance and angular speed.  Classes for
        # real experiments should obviously use ephemeris data.
        rad = np.fmod((start - self._firsttime) * self._radpersec, 2.0 * np.pi)
        ang = self._radinc * np.arange(n, dtype=np.float64) + rad + (0.5 * np.pi)
        x = self._earthspeed * np.cos(ang)
        y = self._earthspeed * np.sin(ang)
        z = np.zeros_like(x)
        return np.ravel(np.column_stack((x, y, z))).reshape((-1, 3))

    def _put_velocity(self, start, vel):
        raise RuntimeError("cannot write data to simulated velocity")
        return
