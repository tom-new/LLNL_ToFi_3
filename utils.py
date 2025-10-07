#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python 3.11.9

"""
Author:   Bernhard Schuberth, LMU Munich, Germany (bernhard.schuberth@lmu.de)
Date:     2019-02-15
Modified: Tom New, The University of Sydney, Australia (tom.new@sydney.edu.au)
Date:     2024-08-12

LLNL_ToFi

Definition of constants.

    Original work Copyright (C) 2019 Bernhard Schuberth (bernhard.schuberth@lmu.de)
    Modified work Copyright (C) 2024 Tom New (tom.new@sydney.edu.au)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""

# -----------------------------------------------------------------------------

from pathlib import Path
import math
import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import griddata

# --------------------------------------------------------------------------
# Define some constants

# Radius of Earth in km for normalization
R_EARTH_KM = 6371.0
R_CMB_KM = 3501.0

# Path to LLNL-GD3-JPS data (relative path) and filenames:
LLNL_PATH = Path("./DATA/")
LLNL_COORD_FILE = "LLNL_G3D_JPS.Tessellated.Coordinates.txt"
LLNL_DEPTH_FILE = "LLNL_G3D_JPS.Layer_Depths_min_avg_max.txt"
LLNL_R_FILE_PREFIX = "R_Matrix_TomoFilt_Layer"
FIREDRAKE_PATH = Path(
    "/Volumes/Grey/firedrake_simulations/HT_4e8/Z22/0Ma/output_0.pvtu"
)

# Number of radial layers in the upper mantle and transition zone
nl_UM_TZ = 18

# Number of grid points in upper mantle and transition zone layers
np_UM_TZ = 40962
# Number of grid points in lower mantle
np_LM = 10242

# total number of entries in model vector
n_m = 1003608


# I/O paths and filenames
OUTPUT_PATH = Path("./OUTPUT_FILES/")
OUTFILE_FILT_PREFIX = "LLNL_G3D_JPS_ToFi_layer"
OUTFILE_PARM_PREFIX = "LLNL_G3D_JPS_Parm_layer"

# --------------------------------------------------------------------------


def dimensional_constants():
    return {
        "T_0": 3700,  # Used for dimensionalising temperature T_nd * T_0 + T_1
        "T_1": 300,  # See above
        "r_max": 2.208,  # Radius of the Earth in nd
        "r_min": 1.208,  # Radius of the CMB in nd
    }


def compute_mean_profile(array):
    """Compute the mean profile of an array over depth layers"""
    return (
        array.reshape(-1, 129)
        .T.mean(axis=1, keepdims=True)
        .repeat(array.shape[0] // 129, axis=1)
        .T.reshape(-1)
    )


def read_coords():

    # Define filename for coordinates of the LLNL_G3D_JPS parametrization
    coord_file = LLNL_PATH / LLNL_COORD_FILE

    f = open(coord_file, "r")
    cl = f.readlines()
    f.close()

    # The coordinate file of LLNL_G3D_JPS provides both geodetic and geocentric latitude
    gd_lat = np.arange(len(cl), dtype="float64")
    gc_lat = np.arange(len(cl), dtype="float64")
    lon = np.arange(len(cl), dtype="float64")
    top_radius = np.arange(len(cl), dtype="float64")

    for i in range(len(cl)):
        cline = cl[i].split()

        gd_lat[i] = float(cline[0])
        lon[i] = float(cline[1])
        gc_lat[i] = float(cline[2])
        top_radius[i] = float(cline[3])

    return gd_lat, lon, gc_lat, top_radius


# --------------------------------------------------------------------------


def read_radii():

    # Define filename for coordinates of the LLNL_G3D_JPS parametrization
    depth_file = LLNL_PATH / LLNL_DEPTH_FILE

    f = open(depth_file, "r")
    cl = f.readlines()
    f.close()

    radii = {}

    for i in range(len(cl)):
        cline = cl[i].split()

        radii[i] = {}
        radii[i]["min"] = R_EARTH_KM - float(cline[2])
        radii[i]["avg"] = R_EARTH_KM - float(cline[1])
        radii[i]["max"] = R_EARTH_KM - float(cline[0])

    return radii


# --------------------------------------------------------------------------


def read_radius():

    # Define filename for coordinates of the LLNL_G3D_JPS parametrization
    depth_file = LLNL_PATH / LLNL_DEPTH_FILE

    f = open(depth_file, "r")
    cl = f.readlines()
    f.close()

    radius = np.arange(len(cl), dtype=float)

    for i in range(len(cl)):
        cline = cl[i].split()

        radius[i] = R_EARTH_KM - float(cline[0])

    return radius


# --------------------------------------------------------------------------
def read_layer_R(ilyr):

    # Open file that contains the entries of R in current layer
    R_file_cl = LLNL_PATH / f"{LLNL_R_FILE_PREFIX}_{ilyr}.txt"

    f = open(R_file_cl, "r")
    cl = f.readlines()
    f.close()

    row_i = np.arange(len(cl), dtype=int)
    column_j = np.arange(len(cl), dtype=int)
    R_ij = np.arange(len(cl), dtype=float)

    for i in range(len(cl)):
        cline = cl[i].split()

        row_i[i] = int(cline[0])
        column_j[i] = int(cline[1])
        R_ij[i] = float(cline[2])

    return row_i, column_j, R_ij


# --------------------------------------------------------------------------
def read_layer(ilyr, radius, PREFIX):

    depth = R_EARTH_KM - radius

    # Open file that contains the entries of R in current layer
    file_cl = "".join([OUTPUT_PATH, PREFIX, "_%02d_d%04dkm.txt" % (ilyr, depth)])

    f = open(file_cl, "r")
    cl = f.readlines()
    f.close()

    lon = np.arange(len(cl), dtype=float)
    lat = np.arange(len(cl), dtype=float)
    value = np.arange(len(cl), dtype=float)

    for i in range(len(cl)):

        cline = cl[i].split()

        if i == 0 and "#" in cline[0]:
            header = cline

        else:
            lon[i] = float(cline[0])
            lat[i] = float(cline[1])
            value[i] = float(cline[2])

    return lon, lat, value, header


# --------------------------------------------------------------------------
def write_layer(ilyr, m_prime, radius, lon, lat, PREFIX, OUT_PATH, string=""):

    depth = int(math.floor(R_EARTH_KM - radius))

    # Open output file
    filt_file_cl = OUT_PATH / f"{PREFIX}_{ilyr:02d}_d{depth:04d}km.txt"

    f = open(filt_file_cl, "w+")

    f.write("#    Radius: %8.3f, Depth: %8.3f %s\n" % (radius, depth, string))

    for i in range(len(m_prime)):

        # Output for gmt (i.e., lon in first column)
        f.write(" %8.3f %8.3f %12.7f\n" % (lon[i], lat[i], m_prime[i]))

    f.close()


# --------------------------------------------------------------------------


def row_index_offset(ilyr):

    if ilyr <= nl_UM_TZ:
        offset = (ilyr - 1) * np_UM_TZ
    else:
        offset = nl_UM_TZ * np_UM_TZ + (ilyr - nl_UM_TZ - 1) * np_LM

    return offset


# --------------------------------------------------------------------------


def calculate_coord_index(cj):

    ntot_UM_TZ = nl_UM_TZ * np_UM_TZ

    if cj <= ntot_UM_TZ:
        c_index = np.mod(cj, np_UM_TZ) - 1
        l_index = cj // np_UM_TZ
    else:
        cj_rem = cj - ntot_UM_TZ
        c_index = np.mod(cj_rem, np_LM) - 1
        l_index = nl_UM_TZ + cj_rem // np_LM

    if cj == n_m:  # Last entry in model vector
        l_index = 43

    return c_index, l_index


# --------------------------------------------------------------------------


def get_coordinates(column_j, radius, gc_lat, lon):

    # Compute indices for current point
    [c_index, l_index] = calculate_coord_index(column_j)

    # Get current coordinates
    clat = gc_lat[c_index]
    clon = lon[c_index]

    # Set current radius
    crad = radius[l_index]

    return crad, clat, clon, c_index, l_index


# --------------------------------------------------------------------------


def parallelize(myrank, num_procs, ntot):

    n_sub = ntot // num_procs

    if myrank == num_procs - 1:
        my_ib = myrank * n_sub
        my_ie = ntot - 1
    else:
        my_ib = myrank * n_sub
        my_ie = (myrank + 1) * n_sub - 1

    return n_sub, my_ib, my_ie


# --------------------------------------------------------------------------


def convert_to_netcdf(txt_path, nc_path):

    def grid_llnl_from_txt(files, grid_lon, grid_lat):
        """Interpolate LLNL text grids onto a target lon/lat grid and return (V, depths)."""
        depths = []
        V = []
        files = sorted(files, key=lambda filename: int(str(filename).split("_")[-2]))
        for file in files:
            # read whitespace-delimited text, skipping the first header row
            arr = np.genfromtxt(file, skip_header=1, usecols=(0, 1, 2), dtype=float)

            lon = arr[:, 0]
            lat = arr[:, 1]
            val = arr[:, 2]

            # ensure longitudes in [-180, 180)
            lon = ((lon + 180.0) % 360.0) - 180.0

            # stitch across the dateline to avoid edge artifacts near ±180°
            left_mask = lon > 180.0 - 8.0
            right_mask = lon < -180.0 + 8.0

            ext_lon = np.concatenate(
                [lon, lon[left_mask] - 360.0, lon[right_mask] + 360.0]
            )
            ext_lat = np.concatenate([lat, lat[left_mask], lat[right_mask]])
            ext_val = np.concatenate([val, val[left_mask], val[right_mask]])

            # interpolate onto the provided grid
            interp = griddata(
                (ext_lon, ext_lat),
                ext_val,
                (grid_lon, grid_lat),
                method="cubic",
            )
            V.append(interp)

            # extract depth from filename token
            depth = float(str(file).split("_")[-1][1:-6])
            depths.append(depth)

        return np.asarray(V), np.asarray(depths)

    lats = np.linspace(-90, 90, 181)
    lons = np.linspace(-180, 179, 360)  # includes both -180 and +179
    grid_lon, grid_lat = np.meshgrid(lons, lats)

    types = ["Parm_layer_p", "ToFi_layer_p", "Parm_layer_s", "ToFi_layer_s"]
    names = [
        "dlnVp_reparam_percent",
        "dlnVp_tofi_percent",
        "dlnVs_reparam_percent",
        "dlnVs_tofi_percent",
    ]

    Vs = []
    for ftype in types:
        files = [
            p
            for p in txt_path.iterdir()
            if p.is_file() and ftype in p.name and not p.name.startswith(("._"))
        ]
        V, depths = grid_llnl_from_txt(files, grid_lon, grid_lat)
        Vs.append(V)

    depths = np.array(depths)
    radii = (R_EARTH_KM - depths) * 1e3

    # # drop the duplicate seam at +180°: slice off the last longitude column in the data
    # for i in range(len(Vs)):
    #     Vs[i] = Vs[i][..., :-1]

    # # and drop the +180° coordinate so lon runs [-180, 179] in 1° steps
    # lons = lons[:-1]

    # set up DataArrays for primary coordinates
    r = xr.DataArray(
        radii,
        dims="r",
        attrs={"long_name": "radius", "units": r"\metre", "positive": "up"},
    )
    lat = xr.DataArray(
        lats, dims="lat", attrs={"long_name": "latitude", "units": r"\degree"}
    )
    lon = xr.DataArray(
        lons,
        dims="lon",
        attrs={"long_name": "longitude", "units": r"\degree", "convention": "bipolar"},
    )

    # create dataset
    ds = xr.Dataset(coords={"r": r, "lat": lat, "lon": lon, "depth": ("r", depths)})

    for i in range(len(Vs)):
        Vs[i] *= 100  # convert to percent
        # assign attributes to depth
        ds["depth"] = ds["depth"].assign_attrs(
            {"long_name": "depth", "units": r"\kilo\metre", "positive": "down"}
        )
        # explicitly specify dims; data shape is (r, lat, lon)
        ds[f"{names[i]}"] = (("r", "lat", "lon"), Vs[i])

        # assign attributes to data
        ds[f"{names[i]}"] = ds[f"{names[i]}"].assign_attrs(
            {"long_name": "Velocity perturbation", "units": r"\percent"}
        )

    # write to disk
    ds.to_netcdf(nc_path)
