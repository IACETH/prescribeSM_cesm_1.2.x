#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Mathias Hauser
# Date: 08.2014

# USAGE: CONVERT clm history files in xy-format (lat, lon) to col format


try:
	import netCDF4 as nc
except: # netCDF4 is not available on cscs
	import scipy.io.netcdf as nc
	nc.Dataset = nc.netcdf_file

import shutil


# DEFINE FILE NAMES

# example file in column/ 3D form (see README)
# this file needs to exist
fN_3D_in = 'file_3D_in.nc'

# input file in lat-lon/ 4D form
# this file needs to exist
# this is what we want to transform
fN_4D = 'file_4D.nc'

# output file: this is the transformed fN_xy
# does not need to exist (will be overwritten)
fN_3D_out = 'file_3D_out.nc'

# START PROGRAM
with nc.Dataset(fN_3D_in) as ncf:

    # get indices
    cols1d_ixy = ncf.variables['cols1d_ixy'][:]
    cols1d_jxy = ncf.variables['cols1d_jxy'][:]
    cols1d_itype_lunit = ncf.variables['cols1d_itype_lunit'][:]

    shape = ncf.variables['SOILLIQ'].shape

istsoil = 1

# not all the points in SOILLIQ and SOILICE are SM
sel_soil = cols1d_itype_lunit == istsoil
# python uses 0 based indexing
col = cols1d_ixy[sel_soil] - 1
row = cols1d_jxy[sel_soil] - 1

# copy col/3D 'example file' to 'output file'
shutil.copyfile(fN_3D_in, fN_3D_out)

# read from fN_xy (only the first 10 lvl are SM levels)
with nc.Dataset(fN_4D) as ncf:
    SOILLIQ_xy = ncf.variables['SOILLIQ'][:, 0:10, :, :]
    SOILICE_xy = ncf.variables['SOILICE'][:, 0:10, :, :]

# write out the transformed SM
with nc.Dataset(fN_3D_out, 'a') as ncf:
    ncf.variables['SOILLIQ'][:, 0:10, sel_soil] = SOILLIQ_xy[:, :, row, col]
    ncf.variables['SOILICE'][:, 0:10, sel_soil] = SOILICE_xy[:, :, row, col]


# LOOP VERSION

# for day in [0]: # replace with 0:365
#     for lev in [0]: # replace with 0:10
#         for i in xrange(n_cols1d):
#             col = cols1d_ixy[i]
#             row = cols1d_jxy[i]

#             if cols1d_itype_lunit[i] == 1:
#                 SL[row, col] = SOILLIQ[day, lev, sel_soil]


# SL_v = np.empty(shape=(13693,))
# SL_v.fill(np.nan)

# for day in [0]:
#     for lev in [0]:
#         for i in xrange(n_cols1d):
#             col = cols1d_ixy[i]
#             row = cols1d_jxy[i]

#             if cols1d_itype_lunit[i] == 1:
#                 SL_v[i] = SOILLIQ_xy[day, lev, row - 1, col - 1]



