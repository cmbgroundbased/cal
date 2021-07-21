# Copyright (c) 2019-2020 by the parties listed in the AUTHORS file.
# All rights reserved.  Use of this source code is governed by
# a BSD-style license that can be found in the LICENSE file.

import os
import pickle
import sys

import numpy as np

from ..timing import function_timer, Timer
from ..tod import AnalyticNoise
from ..utils import Logger, Environment
from .. import qarray


def name2id(name, maxval=2 ** 16):
    """ Map a name into an index.
    """
    value = 0
    for c in name:
        value += ord(c)
    return value % maxval


class Focalplane:
    _detweights = None
    _detquats = None
    _noise = None

    def __init__(
        self, detector_data=None, fname_pickle=None, sample_rate=None, radius_deg=None
    ):
        """ Instantiate a focalplane

        Args:
            detector_data (dict) :  Dictionary of detector attributes, such
                as detector quaternions and noise parameters.
            fname_pickle (str) :  Pickle file containing the focal
                 plane dictionary.  If both `detector_data` and
                 `fname_pickle` are set, the dictionaries are merged.
            sample_rate (float) :  Default sampling rate for all
                detectors.  Will be overridden by 'fsample' fields
                if they exist for the detectors in the dictionary.
            radius_deg (float) :  force the radius of the focal plane.
                otherwise it will be calculated from the detector
                offsets.
        """
        self.detector_data = {}
        if detector_data is not None:
            self.detector_data.update(detector_data)
        if fname_pickle is not None:
            with open(fname_pickle, "rb") as picklefile:
                self.detector_data.update(pickle.load(picklefile))
        self.sample_rate = sample_rate
        self._radius = radius_deg
        self._get_pol_angles()
        self._get_pol_efficiency()

    def _get_pol_angles(self):
        """ Get the detector polarization angles from the quaternions
        """
        for detname, detdata in self.detector_data.items():
            if "pol_angle_deg" not in detdata and "pol_angle_rad" not in detdata:
                quat = detdata["quat"]
                psi = qarray.to_angles(quat)[2]
                detdata["pol_angle_rad"] = psi
        return

    def _get_pol_efficiency(self):
        """ Get the polarization efficiency from polarization leakage
        or vice versa
        """
        for detname, detdata in self.detector_data.items():
            if "pol_leakage" in detdata and "pol_efficiency" not in detdata:
                # Derive efficiency from leakage
                epsilon = detdata["pol_leakage"]
                eta = (1 - epsilon) / (1 + epsilon)
                detdata["pol_efficiency"] = eta
            elif "pol_leakage" not in detdata and "pol_efficiency" in detdata:
                # Derive leakage from efficiency
                eta = detdata["pol_effiency"]
                epsilon = (1 - eta) / (1 + eta)
                detdata["pol_leakage"] = epsilon
            elif "pol_leakage" not in detdata and "pol_efficiency" not in detdata:
                # Assume a perfectly polarized detector
                detdata["pol_efficiency"] = 1
                detdata["pol_leakage"] = 0
            else:
                # Check that efficiency and leakage are consistent
                epsilon = detdata["pol_leakage"]
                eta = detdata["pol_efficiency"]
                np.testing.assert_almost_equal(
                    eta,
                    (1 + epsilon) / (1 - epsilon),
                    err_msg="inconsistent polarization leakage and efficiency",
                )
        return

    def __contains__(self, key):
        return key in self.detector_data

    def __getitem__(self, key):
        return self.detector_data[key]

    def __setitem__(self, key, value):
        self.detector_data[key] = value

    def reset_properties(self):
        """ Clear automatic properties so they will be re-generated
        """
        self._detweights = None
        self._radius = None
        self._detquats = None
        self._noise = None

    @property
    def detweights(self):
        """ Return the inverse noise variance weights [K_CMB^-2]
        """
        if self._detweights is None:
            self._detweights = {}
            for detname, detdata in self.detector_data.items():
                net = detdata["NET"]
                if "fsample" in detdata:
                    fsample = detdata["fsample"]
                else:
                    fsample = self.sample_rate
                detweight = 1.0 / (fsample * net ** 2)
                self._detweights[detname] = detweight
        return self._detweights

    @property
    def radius(self):
        """ The focal plane radius in degrees
        """
        if self._radius is None:
            # Find the largest distance from the bore sight
            ZAXIS = np.array([0, 0, 1])
            cosangs = []
            for detname, detdata in self.detector_data.items():
                quat = detdata["quat"]
                vec = qarray.rotate(quat, ZAXIS)
                cosangs.append(np.dot(ZAXIS, vec))
            mincos = np.amin(cosangs)
            self._radius = np.degrees(np.arccos(mincos))
            # Add a very small margin to avoid numeric issues
            # in the atmospheric simulation
            self._radius *= 1.001
        return self._radius

    @property
    def detquats(self):
        if self._detquats is None:
            self._detquats = {}
            for detname, detdata in self.detector_data.items():
                self._detquats[detname] = detdata["quat"]
        return self._detquats

    @property
    def noise(self):
        if self._noise is None:
            detectors = sorted(self.detector_data.keys())
            fmin = {}
            fknee = {}
            alpha = {}
            NET = {}
            rates = {}
            for detname in detectors:
                detdata = self.detector_data[detname]
                if "fsample" in detdata:
                    rates[detname] = detdata["fsample"]
                else:
                    rates[detname] = self.sample_rate
                fmin[detname] = detdata["fmin"]
                fknee[detname] = detdata["fknee"]
                alpha[detname] = detdata["alpha"]
                NET[detname] = detdata["NET"]
            self._noise = AnalyticNoise(
                rate=rates,
                fmin=fmin,
                detectors=detectors,
                fknee=fknee,
                alpha=alpha,
                NET=NET,
            )
        return self._noise

    def __repr__(self):
        value = (
            "(Focalplane : {} detectors, sample_rate = {} Hz, radius = {} deg, "
            "detectors = ("
            "".format(len(self.detector_data), self.sample_rate, self.radius)
        )
        for detector_name, detector_data in self.detector_data.items():
            value += "{}, ".format(detector_name)
        value += "))"
        return value
    
    def _plot_fp(self, width, height, outfile, fwhm=None, facecolor=None, polcolor=None, labels=None):
        """Visualize a dictionary of detectors.

        This makes a simple plot of the detector positions on the projected
        focalplane.

        To avoid python overhead in large MPI jobs, we place the matplotlib
        import inside this function, so that it is only imported when the
        function is actually called.

        If the detector dictionary contains a key "fwhm", that will be assumed
        to be in arcminutes.  Otherwise a nominal value is used.

        If the detector dictionary contains a key "viscolor", then that color
        will be used.

        Args:
            dets (dict): dictionary of detector quaternions.
            width (float): width of plot in degrees.
            height (float): height of plot in degrees.
            outfile (str): output PNG path.  If None, then matplotlib will be
                used for inline plotting.
            fwhm (dict): dictionary of detector beam FWHM in arcminutes, used
                to draw the circles to scale.
            facecolor (dict): dictionary of color values for the face of each
                detector circle.
            polcolor (dict): dictionary of color values for the polarization
                arrows.
            labels (dict): plot this text in the center of each pixel.

        Returns:
            None

        """
        if outfile is not None:
            import matplotlib
            import warnings

            # Try to force matplotlib to not use any Xwindows backend.
            warnings.filterwarnings("ignore")
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        xfigsize = int(5 * width)
        yfigsize = int(5 * height)
        figdpi = 75

        # Compute the font size to use for detector labels
        fontpix = 0.2 * figdpi
        fontpt = int(0.75 * fontpix)

        fig = plt.figure(figsize=(xfigsize, yfigsize), dpi=figdpi)
        ax = fig.add_subplot(1, 1, 1)

        half_width = 0.5 * width
        half_height = 0.5 * height
        ax.set_xlabel("Degrees", fontsize="large")
        ax.set_ylabel("Degrees", fontsize="large")
        ax.set_xlim([-half_width, half_width])
        ax.set_ylim([-half_height, half_height])

        xaxis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        yaxis = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        zaxis = np.array([0.0, 0.0, 1.0], dtype=np.float64)

        for d, quat in self.detquats.items():

            # radius in degrees
            detradius = 0.5 * 5.0 / 60.0
            if fwhm is not None:
                detradius = 0.5 * fwhm[d] / 60.0

            # rotation from boresight
            rdir = qarray.rotate(quat, zaxis).flatten()
            ang = np.arctan2(rdir[1], rdir[0])

            orient = qarray.rotate(quat, xaxis).flatten()
            polang = np.arctan2(orient[1], orient[0])

            mag = np.arccos(rdir[2]) * 180.0 / np.pi
            xpos = mag * np.cos(ang)
            ypos = mag * np.sin(ang)

            detface = "none"
            if facecolor is not None:
                detface = facecolor[d]

            circ = plt.Circle((xpos, ypos), radius=detradius, fc=detface, ec="k")
            ax.add_artist(circ)

            ascale = 2.0

            xtail = xpos - ascale * detradius * np.cos(polang)
            ytail = ypos - ascale * detradius * np.sin(polang)
            dx = ascale * 2.0 * detradius * np.cos(polang)
            dy = ascale * 2.0 * detradius * np.sin(polang)

            detcolor = "black"
            if polcolor is not None:
                detcolor = polcolor[d]

            ax.arrow(
                xtail,
                ytail,
                dx,
                dy,
                width=0.1 * detradius,
                head_width=0.3 * detradius,
                head_length=0.3 * detradius,
                fc=detcolor,
                ec=detcolor,
                length_includes_head=True,
            )

            if labels is not None:
                xsgn = 1.0
                if dx < 0.0:
                    xsgn = -1.0
                labeloff = 0.05 * xsgn * fontpix * len(labels[d]) / figdpi
                ax.text(
                    (xtail + 1.1 * dx + labeloff),
                    (ytail + 1.1 * dy),
                    labels[d],
                    color="k",
                    fontsize=fontpt,
                    horizontalalignment="center",
                    verticalalignment="center",
                    bbox=dict(fc="w", ec="none", pad=1, alpha=1.0),
                )

        if outfile is None:
            plt.show()
        else:
            plt.savefig(outfile)
            plt.close()
        return fig
        


class Telescope(object):
    def __init__(self, name, focalplane=None, site=None):
        self.name = name
        self.id = name2id(name)
        self.focalplane = focalplane
        self.site = site

    def __repr__(self):
        value = "(Telescope '{}' : ID = {}, Site = {}, Focalplane = {}" "".format(
            self.name, self.id, self.site, self.focalplane
        )
        return value