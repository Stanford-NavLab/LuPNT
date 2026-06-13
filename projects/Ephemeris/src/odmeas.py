import numpy as np
import pandas as pd
import pylupnt as pnt

import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import os


class ODMeas:
    def __init__(self, name):
        self.name = name

    def get_meas(self, t, x):
        # no implementation
        nx = x.shape[0]
        n_meas = 0
        if t.ndim == 0:
            lent = 1
        else:
            lent = t.shape[0]
        y = np.zeros((lent, n_meas))
        H = np.zeros((lent, n_meas, nx))
        R = np.zeros((lent, n_meas, n_meas))

        return y, H, R


def lla2cart(lat, lon, alt, R):
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    x = (R + alt) * np.cos(lat) * np.cos(lon)
    y = (R + alt) * np.cos(lat) * np.sin(lon)
    z = (R + alt) * np.sin(lat)

    return np.array([x, y, z])


class EarthStationMeas(ODMeas):
    def __init__(self, min_elev_deg=10.0):
        super().__init__("EarthStation")
        self.min_elev_deg = min_elev_deg

        self.stations = {
            "LEGS-1": {
                "location": "WhiteSands, USA",
                "lat": 32.544863,
                "lon": 253.387496,
                "alt_m": 1464.0,
            },
            "LEGS-2": {
                "location": "Matjiesfontein, SouthAfrica",
                "lat": -33.231224,
                "lon": 20.58163,
                "alt_m": 899,
            },
            "LEGS-3": {
                "location": "Dongara, Australia",
                "lat": -29.0457,
                "lon": 115.3487,
                "alt_m": 250.6,
            },
        }

        # compute the cartesian coordinates of the stations
        for label, station in self.stations.items():
            lat = station["lat"]
            lon = station["lon"]
            alt = station["alt_m"]
            self.stations["ecef_xyz"] = lla2cart(lat, lon, alt)

    def get_meas(self, t, sat_x, sat_frame="mci", use_range=True, use_doppler=True):
        n_meas = 0
        if use_range:
            n_meas += 3  # for each station
        if use_doppler:
            n_meas += 3  # for each station

        if sat_x.ndim == 1:
            n_sat = 1
        else:
            n_sat = sat_x.shape[0]

        nx = sat_x.shape[0]
        n_meas = 0
