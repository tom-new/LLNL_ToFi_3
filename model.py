#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python 3.11.9

"""
Author:   Bernhard Schuberth, LMU Munich, Germany (bernhard.schuberth@lmu.de)
Date:     2019-02-15
Modified: Tom New, The University of Sydney, Australia (tom.new@sydney.edu.au)
Date:     2024-08-12

LLNL_ToFi

Example routines for determining the values of a seismic velocity model on the 
grid points of the LLNL-G3D-JPS model.

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
#-----------------------------------------------------------------------------

from mpi4py import MPI
import numpy as np
import pyvista as pv
import gdrift
import spherical
from scipy.interpolate import RBFInterpolator, CubicSpline
from pathlib import Path
import sys

import ctypes as C


from utils import (R_EARTH_KM, LLNL_PATH, LLNL_COORD_FILE, LLNL_DEPTH_FILE,
                   LLNL_R_FILE_PREFIX, nl_UM_TZ, np_UM_TZ, np_LM, n_m,
                   OUTPUT_PATH, OUTFILE_FILT_PREFIX, OUTFILE_PARM_PREFIX,
                   FIREDRAKE_PATH)

import utils


#--------------------------------------------------------------------------    
def init_model_parallel(comm=0):

    myrank = comm.Get_rank()
    num_procs = comm.Get_size()

    comm.barrier()

    snd_model = None
    if myrank == 0:
        snd_model = read_model()

    keys = None
    if myrank == 0:
        keys = list(snd_model.keys())
    keys = comm.bcast(keys, root=0)

    rcv_model = {}

    for key in keys:
        snd_array = None
        if myrank == 0:
            snd_array = snd_model[key]
        snd_array = comm.bcast(snd_array, root=0)
    
        rcv_model[key] = snd_array

    model = {
        "du": RBFInterpolator(rcv_model["coords"], rcv_model["du"], neighbors=32, kernel="linear"),
        "v_1D": CubicSpline(rcv_model["radii"], rcv_model["v_1D"])
    }

    comm.barrier()

    return model
  
#--------------------------------------------------------------------------    
def read_model():

    # USER MODIFICATION REQUIRED
    # Please provide the code to read in your model
    model_path = Path(FIREDRAKE_PATH) / Path("Hall2002/Stage_27_Gplates") / Path("output_4.pvtu")
    model = pv.read(model_path)
    model = model.clean()
    model.points *= R_EARTH_KM / 2.22
    model.point_data['T'] = model['FullTemperature'] * 3700 + 300
    model.point_data['dT'] = model['DeltaT'] * (np.max(model['T']) - np.min(model['T']))
    model.point_data['T_av'] = model['T'] - model['dT']
    model.point_data['radii'] = np.linalg.norm(model.points, axis=1)

    model.point_data['depth'] = (R_EARTH_KM - model["radii"]) * 1.0e3

    slb_pyrolite = gdrift.ThermodynamicModel('SLB_16', 'pyrolite')

    model.point_data['vp'] = slb_pyrolite.temperature_to_vp(temperature=np.array(model['T']), depth=np.array(model['depth']))
    model.point_data['vp_av'] = slb_pyrolite.temperature_to_vp(temperature=np.array(model['T_av']), depth=np.array(model['depth']))

    model.point_data['vs'] = slb_pyrolite.temperature_to_vs(temperature=np.array(model['T']), depth=np.array(model['depth']))
    model.point_data['vs_av'] = slb_pyrolite.temperature_to_vs(temperature=np.array(model['T_av']), depth=np.array(model['depth']))

    keys = ["vp", "vp_av"]
    
    coords = np.array(model.points)
    du = np.array(1/model[keys[0]] - 1/model[keys[1]])
    radii = np.linspace(model["radii"].min(), model["radii"].max(), num=250)
    v_1D_coords = np.array([radii, np.zeros_like(radii), np.zeros_like(radii)]).T
    v_1D = pv.PolyData(v_1D_coords)
    v_1D = np.array(v_1D.sample(model)["vp_av"])

    model = {
        "coords": coords,
        "du": du,
        "radii": radii,
        "v_1D": v_1D
    }
    # END USER MODIFICATION REQUIRED

    return model

#--------------------------------------------------------------------------    
def project_slowness_3D(model,radius_avg,lat,lon,radius_min,radius_max,grid_spacing):

    # This is a dummy routine that needs to be modified by the user.

    # Please modify the code to obtain 3-D slowness perturbations (i.e., the absolute difference between 
    # 3-D and 1-D slowness at the current point) for your model.
    # Depending on your model parametrization (i.e., coarser or finer than LLNL-G3D-JPS), you will need
    # to either perform an interpolation (e.g., from the nearest neighbor's on your grid to the current 
    # point in the LLNL-G3D-JPS grid if coarser), or you will have to compute an average value in the
    # volume given by "radius_min" and "radius_max" in vertical direction and grid_spacing as
    # search radius in lateral direction (e.g., by an inverse-distance weighting algorithm).

    # NOTE: Perturbation in slowness is du = (1/v_3D - 1/v_1D), and du is approximately -(v - v_1D)/(v_1D**2) 
    #       => du = -dv/v_1D^2 = -dln(v)/v_1D

    # USER MODIFICATION REQUIRED
    spherical_coord = spherical.geo2sph([radius_avg,lon,lat])
    cart_coord = spherical.sph2cart(spherical_coord)
    du = model["du"]([cart_coord])[0]
    # END USER MODIFICATION REQUIRED

    return du

#--------------------------------------------------------------------------    
def model_1D(model,radius):

    # This is a dummy routine that needs to be modified by the user

    # Please modify the code to obtain the 1-D seismic velocity value for the given radius.


    # USER MODIFICATION REQUIRED
    v_1D = model["v_1D"](radius)
    # END USER MODIFICATION REQUIRED

    return float(v_1D)

#--------------------------------------------------------------------------    
def get_slowness_layer(model,radius_in,lat,lon,grid_spacing):

    # This is a dummy routine that illustrates how to get values of a seismic velocity 
    # model in terms of slowness perturbation du = 1/v_3D - 1/v_1D onto the grid 
    # of the LLNL-G3D-JPS tomographic model. 
    # Note: dv = -du*v_1D^2 => dv/v_1D = dln(v) = -du*v_1D; du = -dln(v)/v_1D


    # USER MODIFICATION REQUIRED
    # This routine expects radius to be given in km.
    # Thus, normalize the radii if necessary (uncomment the line below if applicable).
    #r_norm = R_EARTH_KM
    r_norm = 1. # no radius normalization by default
    # END USER MODIFICATION REQUIRED

    # turn input radius into a vector if not already provided in this form
    if np.size(radius_in["avg"]) != np.size(lat):
        radius_avg = np.ones(len(lat)) * radius_in["avg"] / r_norm
        radius_min = np.ones(len(lat)) * radius_in["min"] / r_norm
        radius_max = np.ones(len(lat)) * radius_in["max"] / r_norm
    else:
        radius_avg = radius_in["avg"] / r_norm
        radius_min = radius_in["min"] / r_norm
        radius_max = radius_in["max"] / r_norm

    # Get 1-D seismic velocity for that layer
    v_1D = model_1D(model,radius_avg[0])

    slowness_perturbation = np.ones(len(lat))
    # Loop over all points
    for ip in np.arange(len(lat)):
            
        # Get 3-D (absolute) slowness perturbation du = 1/v_3D - 1/v_1D
        slowness_perturbation[ip]  = project_slowness_3D(model,radius_avg[ip],lat[ip],lon[ip],radius_min[ip],radius_max[ip],grid_spacing)


    return slowness_perturbation, v_1D


#--------------------------------------------------------------------------    
def reparam(comm,radii,gc_lat,lon,reparam):
    
    myrank = comm.Get_rank()
    num_procs = comm.Get_size()

    # Get number of layers
    nl = len(radii)

    slowness_perturbation = {}
    
    if reparam:
        # USER MODIFICATION REQUIRED
        # Initialize the seismic model (if necessary)
        model = init_model_parallel(comm)
        # END USER MODIFICATION REQUIRED

    v_1D = np.zeros(nl)

    for ilyr in range(1,nl+1):


        # Initialize model vectors for that layer
        if ilyr <= nl_UM_TZ:
            cnp = np_UM_TZ
            # nominal grid spacing is 1 degree in the upper mantle and transition zone
            grid_spacing = 111.
        else:
            cnp = np_LM
            # nominal grid spacing is 2 degree in the lower mantle
            grid_spacing = 222.

        
        slowness_perturbation[ilyr-1] = np.zeros(cnp,dtype='float64')

        if reparam:
    
            if myrank == 0:
                if ilyr == 1:
                    print('#')
                    print('# reparametrising the model...')
                    print('#       ... layer %2d ...' % ilyr)
                elif ilyr == nl:
                    print('#       ... layer %2d' % ilyr)
                else:
                    print('#       ... layer %2d ...' % ilyr)


            # Distribute work load on all processors
            [ cnp_sub, my_ib, my_ie ] = utils.parallelize(myrank,num_procs,cnp)

            m_true = np.zeros(cnp,dtype='float64')
            tmp = np.zeros(cnp,dtype='float64')
        

            # Get slowness and 1-D velocity at current location
            [  tmp[my_ib:my_ie], v_1D_tmp ] = get_slowness_layer(model,radii[ilyr-1], gc_lat[my_ib:my_ie],lon[my_ib:my_ie],
                                                                            grid_spacing)

            v_1D[ilyr-1] = v_1D_tmp
        
            comm.Allreduce([tmp, MPI.DOUBLE], [slowness_perturbation[ilyr-1], MPI.DOUBLE], op = MPI.SUM)

            if myrank == 0:
                # Note: dv = -du*v_1D^2 => dv/v_1D = dln(v) = -du*v_1D
                # reparametrised model (dln(v))
                m_true = -1. * slowness_perturbation[ilyr-1] * v_1D[ilyr-1]

                # Output reparametrised model
                header = '# v1D: %12.7f ' % v_1D[ilyr-1]
                utils.write_layer(ilyr,m_true,radii[ilyr-1]["avg"],lon,gc_lat,OUTFILE_PARM_PREFIX,string=header)
    
        else:
    
            if myrank == 0:
                if ilyr == 1:
                    print('#')
                    print('# Reading the reparametrised model...')
                    print('#       ... layer %2d ...' % ilyr)
                elif ilyr == nl:
                    print('#       ... layer %2d' % ilyr)
                else:
                    print('#       ... layer %2d ...' % ilyr)

            m_true = []
            header = ''
            if myrank == 0: 
                # reparametrised model
                [lon_in, gc_lat_in, m_true, header] = utils.read_layer(ilyr,radii[ilyr-1]["avg"],OUTFILE_PARM_PREFIX)

            m_true = comm.bcast(m_true, root=0)
            header = comm.bcast(header, root=0)

            v_1D[ilyr-1] = header[-1]
    
            # Convert velocity to slowness perturbation du
            # dv = -du*v_1D^2 => dv/v_1D = dln(v) = -du*v_1D, du = -dln(v)/v_1D
            slowness_perturbation[ilyr-1] = -1. * m_true / v_1D[ilyr-1] 


    return slowness_perturbation, v_1D


