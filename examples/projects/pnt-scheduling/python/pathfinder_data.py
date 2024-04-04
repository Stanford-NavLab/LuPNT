# Lunar Pathfinder
# Data relay satellite on orbit around the Moon
# Service Guide

import pylupnt as pnt

# Orbital elements (a, e, i, W, w, M) [km, -, deg, deg, deg, deg]
orbital_elements = [5740, 0.58, 54.856, 0, 86.322, 0]

users = [
    {
        "id": 1,
        "desc": "Low Lunar Orbit (LLO) mission",
        "type": "orbital",
        "altitude": 100,  # [km]
        "orbital_elements": [pnt.R_MOON + 100, 0, 28.5, 0, 0, 0],
        "frame": pnt.MI,
        "duration": 30,  # [days]
        "EIRP": 9,  # [dBW]
        "GT": -23.4,  # [dB/K]
        "data": 500,  # [Mb] per day
        "dara_rate": 123,  # [kbps]
        "contact": 150,  # [min] per day
    },
    {
        "id": 2,
        "desc": "Surface operations",
        "type": "surface",
        "location": [-90, 0, 0],  # [deg, deg, km]
        "duration": 13,  # [days]
        "EIRP": 21.5,  # [dBW]
        "GT": -19.2,  # [dB/K]
        "data": 60e3,  # [Mb] per day
        "data_rate": 1966,  # [kbps]
        "contact": 467,  # [min] per day
    },
    {
        "id": 3,
        "desc": "Autonomous Rover",
        "type": "surface",
        "location": [-74.5, 135, 0],  # [deg, deg, km]
        "duration": 13,  # [days]
        "EIRP": 13,  # [dBW]
        "GT": -23,  # [dB/K]
        "data": 20e3,  # [Mb] per day
        "data_rate": 492,  # [kbps]
        "contact": 529,  # [min] per day
    },
    {
        "id": 4,
        "desc": "Tele-Operated Rover",
        "type": "surface",
        "location": [-80, 45, 0],  # [deg, deg, km]
        "duration": 100,  # [days]
        "EIRP": 26.5,  # [dBW]
        "GT": -6,  # [dB/K]
        "data": 30e3,  # [Mb] per day
        "contact": 248,  # [min] per day
        "data_rate": 1966,  # [kbps]
    },
    {
        "id": 5,
        "desc": "CubeSat on PCO",
        "type": "orbital",
        "altitude": 100,  # [km]
        "orbital_elements": [pnt.R_MOON + 100, 0, 60, 90, 0, 90],
        # "orbital_elements": [pnt.R_MOON + 3000, 0, 75, 90, 0, 90],
        "frame": pnt.MI,
        "duration": 181,  # [days]
        "EIRP": 5.7,  # [dBW]
        "GT": -23.4,  # [dB/K]
        "data": 500,  # [Mb] per day
        "contact": 100,  # [min] per day
        "data_rate": 61,  # [kbps]
    },
    {
        "id": 6,
        "desc": "Lander on the North Pole",
        "type": "surface",
        "location": [90, 0, 0],  # [deg, deg, km]
        "duration": 13,  # [days]
        "EIRP": 21.5,  # [dBW]
        "GT": -19.2,  # [dB/K]
        "data": 3e3,  # [Mb] per day
        "contact": 23,  # [min] per day
        "data_rate": 1966,  # [kbps]
    },
]
