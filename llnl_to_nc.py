import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import griddata
from pathlib import Path

from utils import (
    R_EARTH_KM,
    LLNL_PATH,
    LLNL_COORD_FILE,
    LLNL_DEPTH_FILE,
    LLNL_R_FILE_PREFIX,
    nl_UM_TZ,
    np_UM_TZ,
    np_LM,
    n_m,
    OUTPUT_PATH,
    OUTFILE_FILT_PREFIX,
    OUTFILE_PARM_PREFIX,
)


def grid_llnl_from_txt(files, grid_lon, grid_lat):
    depths = []
    V = []
    files = sorted(files, key=lambda filename: int(str(filename).split("_")[-2]))
    for file in files:
        data = pd.read_csv(
            file, sep=r"\s+", skiprows=1, header=None, names=["lon", "lat", "V"]
        )
        data["lon"] = (
            (data["lon"] + 180) % 360
        ) - 180  # ensure longitudes in [-180, 180)
        left = data[data["lon"] > 180.0 - 8].copy()
        left["lon"] -= 360.0
        right = data[data["lon"] < -180.0 + 8].copy()
        right["lon"] += 360.0
        data = pd.concat([data, left, right], ignore_index=True)
        data = griddata(
            data[["lon", "lat"]].to_numpy(),
            data["V"].to_numpy(),
            (grid_lon, grid_lat),
            method="cubic",
        )
        V.append(data)
        depth = float(str(file).split("_")[-1][1:-6])
        depths.append(depth)
    return np.array(V), np.array(depths)


lats = np.linspace(-90, 90, 181)
lons = np.linspace(-180, 179, 360)  # includes both -180 and +179
grid_lon, grid_lat = np.meshgrid(lons, lats)

model = "HTK"
reconstruction = "Z22"
root_path = Path(
    f"/Volumes/Grey/firedrake_simulations/{model}/{reconstruction}/LLNL_ToFi_3"
)
types = ["Parm_layer_p", "ToFi_layer_p", "Parm_layer_s", "ToFi_layer_s"]
names = [
    "dVp_reparam_percent",
    "dVp_tofi_percent",
    "dVs_reparam_percent",
    "dVs_tofi_percent",
]

Vs = []
for type in types:
    files = [p for p in root_path.iterdir() if p.is_file() and type in p.name]
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
    radii, dims="r", attrs={"long_name": "radius", "units": r"\metre", "positive": "up"}
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
ds = xr.Dataset(
    coords={"r": r, "lat": lat, "lon": lon, "depth": ("r", depths)},
    attrs={"id": f"{reconstruction} LLNL ToFi"},
)

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
write_path = Path.home() / Path(
    f"OneDrive/phd/firedrake-models/{model}_{reconstruction}_LLNL_ToFi.nc"
)
ds.to_netcdf(write_path)
