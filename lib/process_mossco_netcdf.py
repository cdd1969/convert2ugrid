#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# This file is part of "convert2ugrid" tool
#
# Author: Nikolai Chernikov, nikolai.chernikov.ru@gmail.com
#
# version = 0.1

'''
This module contains functions for processing data stored in
MOSSCO output netcdf file
'''

from __future__ import division
from netCDF4 import Dataset
import numpy as np
import time
import sys
import os
import matplotlib as plt

from pylab import *



def read_mossco_nc_1d(filename, varname):

    fullname = filename


    root_grp = Dataset(fullname, mode='r')
    Vars = root_grp.variables
    Dims = root_grp.dimensions

    v = Vars[varname]
    dims_n = v.dimensions
    dims_v = v.shape

    a = np.zeros(dims_v)
    a[:] = v[:]  # copying values from variable to numpy array

    root_grp.close()
    print 'read_mossco_nc_1d: returning {0} of shape {1}'.format(type(a), a.shape)
    return a, dims_v


def read_mossco_nc_2d(filename, varname, mask=None):

    fullname = filename

    root_grp = Dataset(fullname, mode='r')
    Vars = root_grp.variables
    Dims = root_grp.dimensions

    v = Vars[varname]

    if mask is None:
        dims_n = v.dimensions
        dims_v = v.shape

        ######
        # BRUTAL HARDCODE!!!!
        ######
        a = np.zeros(dims_v)
        a[:] = v[:].T  # copying values from variable to numpy array
        a = a.flatten(order='F')

    elif mask is not None:
        n_valid_2d = np.sum(np.invert(mask))  #number of valid elements in 2d part. invert - because True is an invalid element
        dims_v = tuple([n_valid_2d])
        
        var_masked = np.ma.array(v, mask=mask.T).T
        var_masked = var_masked.flatten(order='F').compressed()
        a = var_masked

    root_grp.close()
    print 'read_mossco_nc_2d: returning {0} of shape {1}'.format(type(a), a.shape)
    return  a, dims_v


def read_mossco_nc_3d(filename, varname, mask=None):

    fullname = filename

    root_grp = Dataset(fullname, mode='r')
    Vars = root_grp.variables
    Dims = root_grp.dimensions
    print varname
    v = Vars[varname]

    if mask is None:  # IF NOT MASKED ARRAY
        if v.dimensions[0] in ['nMesh2_data_time', 'time']:
            dims_v = tuple([v.shape[0], v.shape[1]*v.shape[2]])
            a = np.zeros(dims_v)
            for i in xrange(v.shape[0]):
                a[i, :] = v[i, ...].T.flatten(1)

        else:
            msg = "Attempting to write 3D variable which has no ['nMesh2_data_time', 'time'] dimension.\nNot implemented yet"
            raise KeyError(msg)

    elif mask is not None:  # IF MASKED ARRAY
        n_valid_2d = np.sum(np.invert(mask))  #number of valid elements in 2d part. invert - because True is an invalid element
        #print 'read_mossco_nc_3d: working with masked array. Mask shape {0}. N_valid elements {1}'.format(mask.shape, n_valid_2d)
        #print 'read_mossco_nc_3d: data shape {0}'.format(v.shape)
        if v.dimensions[0] in ['nMesh2_data_time', 'time']:
            dims_v = tuple([v.shape[0], n_valid_2d])
            a = np.zeros(dims_v)
            for i in xrange(v.shape[0]):
                var_masked = np.ma.array(v[i, ...], mask=mask.T).T
                var_masked = var_masked.flatten(order='F').compressed()
                a[i, :] = var_masked
                
                #plt.imshow(layer_thicness)
                #plt.colorbar()
                #plt.show()

        #--------------------------------------------------------------------------------------
        #--------------------------------------------------------------------------------------
        #--------------------------HARDCODE BELOW WARNIING-----------------------------------
        #--------------------------------------------------------------------------------------
        #--------------------------------------------------------------------------------------
        elif v.dimensions[0] in ['getmGrid3D_getm_3']:  # if we have variable (z, y, x)
            if varname in ['getmGrid3D_getm_layer', 'getmGrid3D_getm_z']:  # HARDCODE
                # we know that in mossco output this variable is 3D (z,y,x) but we need it as
                # a timedependent variable in davit.... therefore adding time dimension
                layer_thicness0 = np.ma.array(Vars['water_depth_at_soil_surface'][0, ...] , mask=mask.T).T

                #layer_thicness = layer_thicness.flatten(order='F').compressed()
                #for i in xrange(len(layer_thicness)):
                #    layer_thicness[i] = layer_thicness[i]/float(v.shape[0])
                # now we have an approximate layer thickness...

                dims_v = tuple([Vars['time'].size, v.shape[0], n_valid_2d])
                a = np.zeros(dims_v)
                print a.shape
                print layer_thicness0.shape
                for t in xrange(a.shape[0]):
                    i = 0
                    for z in xrange(a.shape[1]):
                        var_masked = np.ma.array(v[z, ...], mask=mask.T).T
                        var_masked = var_masked.flatten(order='F').compressed()
                        #layer_thicness = layer_thicness0*(1-i/30.)
                        i += 1
                        a[t, z, :] = var_masked
                        
                        #a[t, z, :] = layer_thicness.flatten(order='F').compressed()
                        #plt.imshow(layer_thicness)
                        #plt.colorbar()
                        #plt.show()
            else:
                msg = "Attempting to write 3D variable which first dimension is 'getmGrid3D_getm_3'\nNot implemented yet. Check code"
                raise KeyError(msg)
        #--------------------------------------------------------------------------------------
        #--------------------------------------------------------------------------------------
        #--------------------------HARDCODE ABOVE WARNIING-----------------------------------
        #--------------------------------------------------------------------------------------
        #--------------------------------------------------------------------------------------
        else:
            msg = "Attempting to write 3D variable which first dimension is not any of ['nMesh2_data_time', 'time', 'getmGrid3D_getm_3']\nNot implemented yet"
            raise KeyError(msg)

    root_grp.close()
    print 'read_mossco_nc_3d: returning {0} of shape {1}'.format(type(a), a.shape)
    return a, dims_v


def read_mossco_nc_4d(filename, varname, mask=None):
    # this will work only for 4d variables which have dimensions (time, z, y, x)
    #
    fullname = filename

    root_grp = Dataset(fullname, mode='r')
    Vars = root_grp.variables
    Dims = root_grp.dimensions

    v = Vars[varname]

    if mask is None:  # IF NOT MASKED ARRAY
        if v.dimensions[0] in ['nMesh2_data_time', 'time']:
            if v.dimensions[1] in ['getmGrid3D_getm_3']:

                    dims_v = tuple([v.shape[0], v.shape[1], v.shape[2]*v.shape[3]])
                    a = np.zeros(dims_v)
                    for t in xrange(v.shape[0]):
                        for z in xrange(v.shape[1]):
                            a[t, z, :] = v[t, z, ...].T.flatten(1)
            else:
                msg = "Attempting to write 4D variable which has no ['getmGrid3D_getm_3'] as 2nd dimension.\nNot implemented yet"
                raise KeyError(msg)
        else:
            msg = "Attempting to write 4D variable which has no ['nMesh2_data_time', 'time'] as 1st dimension.\nNot implemented yet"
            raise KeyError(msg)

    elif mask is not None:  # IF MASKED ARRAY
        n_valid_2d = np.sum(np.invert(mask))  #number of valid elements in 2d part. invert - because True is an invalid element
        #print 'read_mossco_nc_3d: working with masked array. Mask shape {0}. N_valid elements {1}'.format(mask.shape, n_valid_2d)
        #print 'read_mossco_nc_3d: data shape {0}'.format(v.shape)
        if v.dimensions[0] in ['nMesh2_data_time', 'time']:
            if v.dimensions[1] in ['getmGrid3D_getm_3']:
                dims_v = tuple([v.shape[0], v.shape[1], n_valid_2d])
                a = np.zeros(dims_v)
                for t in xrange(v.shape[0]):
                    for z in xrange(v.shape[1]):
                        var_masked = np.ma.array(v[t, z, ...], mask=mask.T).T
                        var_masked = var_masked.flatten(order='F').compressed()
                        a[t, z, :] = var_masked

            else:
                msg = "Attempting to write 4D variable which has no ['getmGrid3D_getm_3'] as 2nd dimension.\nNot implemented yet"
                raise KeyError(msg)
        else:
            msg = "Attempting to write 4D variable which has no ['nMesh2_data_time', 'time'] as 1st dimension.\nNot implemented yet"
            raise KeyError(msg)

    root_grp.close()
    #return a.flatten(1), dims_v
    print 'read_mossco_nc_4d: returning {0} of shape {1}'.format(type(a), a.shape)
    return a, dims_v


def read_mossco_lon_vector(filename, varname='lonc'):

    #path = os.path.dirname(sys.argv[0])
    #fullname = os.path.join(path, filename)
    fullname = filename

    root_grp = Dataset(fullname, mode='r')
    lon_nc = root_grp.variables[varname]
    lon = np.squeeze(lon_nc[:])
    
    root_grp.close()
    return lon


def read_mossco_lat_vector(filename, varname='latc'):

    #path = os.path.dirname(sys.argv[0])
    #fullname = os.path.join(path, filename)
    fullname = filename

    root_grp = Dataset(fullname, mode='r')
    lat_nc = root_grp.variables[varname]
    lat = np.squeeze(lat_nc[:])
    
    root_grp.close()
    return lat


def make_mask_array_from_mossco_bathymetry(filename, varname='bathymetry', fillvalue=None, transpose=False):
    ''' Function reads an 2D array (in MOSSCO - 'bathymetry') and generates a boolean mask
    array (True means that value is masked), based on 'fillvalue'

    input:
        filename        - string to filename (relative path)
        varname         - string containing name of the variable
        fillvalue       - float indicating the values to be masked (mask=True). By default (fillvalue=None)
                        gets _FillValue automatically
        transpose       - a boolean flag to return transposed array or not
    '''
    #path = os.path.dirname(sys.argv[0])
    #fullname = os.path.join(path, filename)
    fullname = filename

    root_grp = Dataset(fullname, mode='r')
    v_nc = root_grp.variables[varname]
    if fillvalue is None:
        try:
            fillvalue = v_nc._FillValue
        except:
            print 'make_mask_array_from_mossco_bathymetry(): attribute "_FillValue" not found. Looking for "missing_value"'
            fillvalue = v_nc.missing_value
    v = np.squeeze(v_nc[:])
    # now get masked array and its mask ( if value is equal to FillValue we assume that is False)
    msk = np.ma.masked_values(v, fillvalue).mask
    if transpose:
        msk = msk.T
    root_grp.close()
    return msk


def get_number_of_depth_layer_from_mossco(list_with_filenames, dimname='getmGrid3D_getm_3'):
    nLayers = False
    for nc_file in list_with_filenames:
        try:
            root_grp = Dataset(nc_file, mode='r')
            d_nc = root_grp.dimensions[dimname]
            nLayers = d_nc.__len__()
            root_grp.close()
        except KeyError:  # if dimension is not found in cuurent file > skip to next file
            pass
    if nLayers:
        print 'found vertical-layers:', nLayers, '\t(from variable <{0}>)'.format(dimname)
        return nLayers
    else:
        raise ValueError('Vertical layers not found. Variable <{0}> not found in files: {1}'.format(dimname, list_with_synoptic_nc))


def get_davit_friendly_variables(filename, tdim=['time'], zdim=['getmGrid3D_getm_3'],
                                    ydim=['getmGrid3D_getm_2', 'getmGrid2D_getm_2', 'y', 'yc', 'lat', 'latc'],
                                    xdim=['getmGrid3D_getm_1', 'getmGrid2D_getm_1', 'x', 'xc', 'lon', 'lonc'],
                                    log=False):
    '''
    function searches for 2D, 3D, 4D variables in mossco netcdf file , that are davit-friendly.
    Let Davit-friendly variables be those, which have dimensions:
        (tdim)
        (ydim, xdim)
        (tdim, ydim, xdim)
        (tdim, zdim, ydim, xdim)

                # (nMesh2_data_time)
                # (nMesh2_data_time, nMesh2_face)
                # (nMesh2_data_time, nMesh2_layer_2d, nMesh2_face)
                # (nMesh2_data_time, nMesh2_layer_3d, nMesh2_face)
                # (nMesh2_data_time, nMesh2_suspension_classes, nMesh2_layer_2d, nMesh2_face)
                # (nMesh2_data_time, nMesh2_suspension_classes, nMesh2_layer_3d, nMesh2_face)
                
    Where "tdim, zdim, ydim, xdim" are dimensions which are contained in the list specified as input parameter
    Returns a dictionary with varibles names as strings.
    '''
    FRIENDLY_VAR_DICT = dict()
    _1D = list()
    _2D = list()
    _3D = list()
    _4D = list()

    #path = os.path.dirname(sys.argv[0])
    #fullname = os.path.join(path, filename)
    fullname = filename

    root_grp = Dataset(fullname, mode='r')
    for var_name, var in root_grp.variables.iteritems():
        if (len(var.shape) == 1) and (var.dimensions[0] in tdim):
            _1D.append(var_name)
        if (len(var.shape) == 2) and (var.dimensions[0] in ydim) and (var.dimensions[1] in xdim):
            _2D.append(var_name)
        elif (len(var.shape) == 3) and (var.dimensions[0] in tdim) and (var.dimensions[1] in ydim) and (var.dimensions[2] in xdim):
            _3D.append(var_name)
        elif (len(var.shape) == 4) and (var.dimensions[0] in tdim) and (var.dimensions[1] in zdim) and (var.dimensions[2] in ydim) and (var.dimensions[3] in xdim):
            _4D.append(var_name)
    root_grp.close()
    FRIENDLY_VAR_DICT['1D'] = _1D
    FRIENDLY_VAR_DICT['2D'] = _2D
    FRIENDLY_VAR_DICT['3D'] = _3D
    FRIENDLY_VAR_DICT['4D'] = _4D
    
    if log:
        print '\n'+'-'*50
        print 'Displaying Davit-friendly variables found in MOSSCO output netcdf file'
        print 'Let Davit-friendly variables be those, which have dimensions:\n\t(tdim)\n\t(ydim, xdim)\n\t(tdim, ydim, xdim)\n\t(tdim, zdim, ydim, xdim)'
        print '-'*50
        i = 1
        for typ in ['1D', '2D', '3D', '4D']:
            for v in FRIENDLY_VAR_DICT[typ]:
                print i, ')\t', typ, ':', v
                i += 1
    return FRIENDLY_VAR_DICT
