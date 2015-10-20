#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# This file is part of "convert2ugrid" tool
#
# Author: Nikolai Chernikov, nikolai.chernikov.ru@gmail.com


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



def find_coordinate_vars(filename, bathymetry_vname='bathymetry', x_vname=None, y_vname=None):
    x_cartesian_vnames = ['x', 'xx', 'x_x' , 'xX', 'x_X', 'xt', 'x_t', 'xT', 'x_T']
    x_geographi_vnames = ['lon', 'lonx', 'lon_x' , 'lonX', 'lon_X', 'lont', 'lon_t', 'lonT', 'lon_T']
    y_cartesian_vnames = ['y', 'yx', 'y_x' , 'yX', 'y_X', 'yt', 'y_t', 'yT', 'y_T']
    y_geographi_vnames = ['lat', 'latx', 'lat_x' , 'latX', 'lat_X', 'latt', 'lat_t', 'latT', 'lat_T']
    
    _n = 'find_coordinate_vars():'
    print '-'*50
    print _n, 'Now searching for coordinates to create grid (X and Y coordinate variables)'
    print _n, 'Surching in file <{0}>'.format(filename)
    root_grp = Dataset(filename, mode='r')
    Vars = root_grp.variables

    if x_vname is not None and y_vname is not None:
        print _n, 'User has explicitly passed varnames for X:<{0}> and for Y:<{1}>'.format(x_vname, y_vname)
        if x_vname not in Vars.keys():
            raise KeyError('find_coordinate_vars(): Variable <{0}> not found in current file'.format(x_vname))
        if y_vname not in Vars.keys():
            raise KeyError('find_coordinate_vars(): Variable <{0}> not found in current file'.format(y_vname))
        print _n, 'Picking them now...'
        
        x_v_n = x_vname
        y_v_n = y_vname
        #return x, y
    else:
        if bathymetry_vname not in Vars.keys():
            raise KeyError('find_coordinate_vars(): Variable <{0}> used for determining X-Y coordinates not found in file <{1}>.\
                            \nTry passing manually var_names'.format(bathymetry_vname, filename))
        b = Vars[bathymetry_vname]
        bath_x_dimName = b.dimensions[1]
        bath_y_dimName = b.dimensions[0]

        # now try to get these variables , same as dim_names
        if bath_x_dimName in Vars.keys() and bath_y_dimName in Vars.keys():
            print _n, 'Using variables that have same names as dimensions of bathymetry variable <{0}>'.format(bathymetry_vname)
            print _n, '\t for X:<{0}>, for Y:<{1}>'.format(bath_x_dimName, bath_y_dimName)
            x_v_n = bath_x_dimName
            y_v_n = bath_y_dimName
        else:
            # nothing has been found... so try to match names from lists
            print _n, 'Bathymetry variable <{2}> has dimensions (Y,X): <({0}, {1})>.'.format(bath_y_dimName, bath_x_dimName, bathymetry_vname)
            print _n, 'But no variables found with these Y,X names. Now trying to search for other coord-vars'

            x_valid_names = list()  # list with X variable found...
            for x_vname in x_cartesian_vnames:
                if x_vname in Vars.keys():
                    x_valid_names.append(x_vname)
            for x_vname in x_geographi_vnames:
                if x_vname in Vars.keys():
                    x_valid_names.append(x_vname)


            y_valid_names = list()  # list with Y variables found
            for y_vname in y_cartesian_vnames:
                if y_vname in Vars.keys():
                    y_valid_names.append(y_vname)
            for y_vname in y_geographi_vnames:
                if y_vname in Vars.keys():
                    y_valid_names.append(y_vname)

            # at this point two lists <x_valid_names>, <y_valid_names> should contain
            # names of the variables to fetch coordinate data
            # There could happen 3 cases:
            #  1) nothing has been found => raise Error
            #  2) one var has been found for x and for y => take them
            #  3) more then one var has been found for x and for y => ask user to choose one

            # case (1)
            if len(x_valid_names) == 0 or len(y_valid_names) == 0:
                raise KeyError('find_coordinate_vars(): No coordinate variables found in current file. Try passing manually var_names')
            
            # case (2)
            elif len(x_valid_names) == 1 and len(y_valid_names) == 1:
                x_v_n = x_valid_names[0]
                y_v_n = y_valid_names[0]
            # case (3)
            else:
                if len(x_valid_names) > 1:
                    print _n, 'More than 1 X-coords variable exist; variables found:'
                    for x_v_n in x_valid_names:
                        print _n, '\t', x_v_n
                    # now promt user to choose correct one
                    x_v_n = None
                    while x_v_n not in Vars.keys():
                        x_v_n = raw_input(_n+'Type desired X-coords variable name:')
                else:
                    x_v_n = x_valid_names[0]

                if len(y_valid_names) > 1:
                    print _n, 'More than 1 Y-coords variable exist; variables found:'
                    for y_v_n in y_valid_names:
                        print _n, '\t', y_v_n
                    
                    # now promt user to choose correct one
                    y_v_n = None
                    while y_v_n not in Vars.keys():
                        y_v_n = raw_input(_n+'Type desired Y-coords variable name:')
                else:
                    y_v_n = y_valid_names[0]

    # -------------------------------------------------------
    # now we know var-names for X and Y
    # -------------------------------------------------------

    # determine mode.... (Geographic or Cartesian)
    if x_v_n in x_cartesian_vnames and y_v_n in y_cartesian_vnames:
        coord_mode = 'cartesian'
    elif x_v_n in x_geographi_vnames and y_v_n in y_geographi_vnames:
        coord_mode = 'geographic'
    else:
        # now show user found vars
        print _n, 'Using variable for X-coords <{0}>, of shape <{1}>'.format(x_v_n, Vars[x_v_n].shape)
        print _n, 'Using variable for Y-coords <{0}>, of shape <{1}>'.format(y_v_n, Vars[y_v_n].shape)
        while coord_mode not in ['cartesian', 'geographic', 'c', 'g']:
            coord_mode = raw_input('\n'+_n+'Coord-type not understood. Choose cartesian or geographic. Type [c/g]:')
        if coord_mode == 'c': coord_mode = 'cartesian'
        if coord_mode == 'g': coord_mode = 'geographic'


    # now determine if it is rectangular grid or curvilinear....
    #
    # if x-coord and y-coord are 1d arrays
    # then the grid is rectangular but the cells may be of
    # any rectangular shape
    #
    # if x-y- coords are 2d arrays,
    # then the grid is curvilinear

    x = Vars[x_v_n]
    y = Vars[y_v_n]

    if len(x.shape) == 1 and len(y.shape) == 1:
        grid_type = 'rectangular'
    elif len(x.shape) == 2 and len(y.shape) == 2:
        grid_type = 'curvilinear'
    else:
        print _n, 'Using variable for X-coords <{0}>, of shape <{1}>'.format(x_v_n, x.shape)
        print _n, 'Using variable for Y-coords <{0}>, of shape <{1}>'.format(y_v_n, y.shape)
        raise ValueError('Grid type not understood. X and Y should be either two 1D arrays or two 2D arrays')


    # now determine if coords are at T (cell center) or X (cell nodes) points...
    if grid_type == 'rectangular':
        if x.shape[0] == b.shape[1] and y.shape[0] == b.shape[0]:  # same as bathymetry
            data_location = 'T_points'
        elif x.shape[0] == (b.shape[1]+1) and y.shape[0] == (b.shape[0]+1):  #+1 more
            data_location = 'X_points'
        else:
            print _n, 'Using bathymetry variable <{0}>, of shape <{1}>'.format(bathymetry_vname, b.shape)
            print _n, 'Using variable for X-coords <{0}>, of shape <{1}>'.format(x_v_n, x.shape)
            print _n, 'Using variable for Y-coords <{0}>, of shape <{1}>'.format(y_v_n, y.shape)
            raise ValueError('Invalid variable dimensions!\nLength of x- or y- dimension in X- or Y- coord-variable should be equal (T_points) or 1 more (X_points) than in bathymetry file')
    else:  #grid_type == 'curvilinear'
        if (x.shape[0] == b.shape[0]) and (x.shape[1] == b.shape[1]
            ) and (y.shape[0] == b.shape[0]) and (y.shape[1] == b.shape[1]):  # same as bathymetry
            data_location = 'T_points'
        elif (x.shape[0] == (b.shape[0]+1)) and (x.shape[1] == (b.shape[1]+1)
            ) and (y.shape[0] == (b.shape[0]+1)) and (y.shape[1] == (b.shape[1]+1)):
            data_location = 'X_points'
        else:
            print _n, 'Using bathymetry variable <{0}>, of shape <{1}>'.format(bathymetry_vname, b.shape)
            print _n, 'Using variable for X-coords <{0}>, of shape <{1}>'.format(x_v_n, x.shape)
            print _n, 'Using variable for Y-coords <{0}>, of shape <{1}>'.format(y_v_n, y.shape)
            raise ValueError('Invalid variable dimensions!\nLength of x- or y- dimension in X- or Y- coord-variable should be equal (T_points) or 1 more (X_points) than in bathymetry file')
    # now print summary....
    print _n, '-'*50
    print _n, 'Summary'
    print _n, '-'*50
    print _n, '\tbathymetry: <{0}>, of shape <{1}>'.format(bathymetry_vname, b.shape)
    print _n, '\tX-coords: <{0}>, of shape <{1}>, range <[{2:.2f}:{3:.2f}]>, units <{4}>'.format(
                x_v_n, x.shape, x[:].min(), x[:].max(), x.units if 'units' in x.ncattrs() else 'unknown')
    print _n, '\tY-coords: <{0}>, of shape <{1}>, range <[{2:.2f}:{3:.2f}]>, units <{4}>'.format(
                y_v_n, y.shape, y[:].min(), y[:].max(), y.units if 'units' in y.ncattrs() else 'unknown')
    print _n, '\tCoordinate mode: <{0}>'.format(coord_mode)
    print _n, '\tGrid_type: <{0}>'.format(grid_type)
    print _n, '\tData location: <{0}>'.format(data_location)
    print _n, '-'*50

    x = x[:]
    y = y[:]
    root_grp.close()
    return {'x' : x, 'y' : y, 'coord_mode' : coord_mode, 'grid_type' : grid_type, 'data_location' : data_location}



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




def read_mossco_nc_rawvar(filename, varname, mask=None):
    root_grp = Dataset(filename, mode='r')
    a = root_grp.variables[varname][:]
    root_grp.close()
    return a


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
    _n = 'get_number_of_depth_layer_from_mossco():'
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
        print _n, 'found vertical-layers:', nLayers, '\t(from dimension <{0}>)'.format(dimname)
        return nLayers
    else:
        print _n, 'Vertical layers not found. Variable <{0}> not found in files: {1}'.format(dimname, list_with_filenames)
        answer = None
        while answer not in ['y', 'Y', 'Yes', 'yes', 'n', 'N', 'no', 'No', 'NO']:
            answer = raw_input(_n+' Will continue with 2D. Explicitly set nLayers=1 [y/n]?:')
        if answer in ['y', 'Y', 'Yes', 'yes']:
            return 1
        else:
            sys.exit(1)


def get_davit_friendly_variables(filename, tdim=['time'], zdim=['getmGrid3D_getm_3'],
        ydim=['getmGrid3D_getm_2', 'getmGrid2D_getm_2', 'y', 'yx', 'y_x', 'yX', 'y_X', 'yt', 'y_t', 'y_T', 'yT' 'yc', 'y_c', 'y_C'
              'lat', 'latc', 'latx', 'latt', 'lat_c', 'lat_x', 'lat_t'],
        xdim=['getmGrid3D_getm_1', 'getmGrid2D_getm_1', 'x', 'xx', 'x_x', 'xX', 'x_X', 'xt', 'x_t', 'x_T', 'xT' 'xc', 'x_c', 'x_C'
              'lon', 'lonc', 'lonx', 'lont', 'lon_c', 'lon_x', 'lon_t'],
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


def get_sigma_coordinates(list_with_filenames, nLayers, varname='level', log=False):
    '''
    get sigma coordinates information. Information is read from the given variable,
    and compared first to the number of layers.

    Output is a 1D array of length <nLayers> or <nLayers+1>
        if the length is <nLayers> then the values represent the sigma coords of the layer centers
        if the length is <nLayers+1> then the values represent the sigma coords of the layer borders
        

    Example: sigma = [-1., -0.75, -0.5, -0.25, 0.], with nLayers=4 it means that we have 4 layers 25% each,
                    and the values describes the border coords. In this case output will look like...
    Output:
        [[-1., -0.75, -0.5, -0.25, 0.], 'border']

    '''
    
    sigma_found = False
    sigma_type = None
    
    for nc_file in list_with_filenames:
        root_grp = Dataset(nc_file, mode='r')
        if varname in root_grp.variables.keys():
            var_nc = root_grp.variables[varname]
            if nLayers == (var_nc.__len__()-1):  # we have sigma coordinates of the layer borders (N_borders = N_layers+1)
                sigma_found = True
                sigma_type = 'border'
                sigma = var_nc[:]
                root_grp.close()
                break
            
            elif nLayers == (var_nc.__len__()):  # we have sigma coordinates of the layer center (N_centers = N_layers)
                sigma_found = True
                sigma_type = 'center'
                sigma = var_nc[:]
                root_grp.close()
                break

        root_grp.close()

        
    if sigma_found:
        print 'get_sigma_coordinates(): Found sigma-coords:', sigma, '\n From variable <{0}>)'.format(varname)
        print 'get_sigma_coordinates(): Values represent coordinates of layer <{0}>'.format(sigma_type)
        print 'get_sigma_coordinates(): Number of layers <{0}>; Number of sigma coords found <{1}>'.format(nLayers, sigma.__len__())
        return [sigma, sigma_type]
    else:
        raise ValueError('Sigma coordinates not found. Variable <{0}> not found in files: {1}'.format(varname, list_with_filenames))



def get_water_level(list_with_filenames, varname='water_level', water_depth_at_soil_surface=None, bathymetry=None, log=False):
    '''
    get water_level information from variable, or generate it from....
        <water_level> = <water_depth_at_soil_surface> - <bathymetry>
    here we assume that <bathymetry> and <water_level> are with respect to the Mean Sea Level, and
    <water_depth_at_soil_surface> is always positive distance between sea bed and water surface

    flatten - if True, x,y dimensions of the array will be compressed into one single
    mask    - boolean 2d mask to ignore elements during flattening. see <process_mossco_netcdf.make_mask_array_from_mossco_bathymetry()>

    '''
    
    water_level = None
    
    for nc_file in list_with_filenames:
        root_grp = Dataset(nc_file, mode='r')
        if varname in root_grp.variables.keys():
            nc_file_found = nc_file
            water_level = root_grp.variables[varname][:]
            root_grp.close()
            break
    

    if not water_level:  # if <water_level> has not been found... generate it!
        if water_depth_at_soil_surface is None or bathymetry is None:  # check inputs
            msg = 'get_water_level(): Searching for water_level data. Variable <{0}> not found in none of the following files:'.format(varname)
            for nc_file in list_with_filenames:
                msg = msg+'\n\t'+nc_file
            msg = msg+'\n Therefore function tried to calculate water_level as <water_level(t,y,x)> = <water_depth_at_soil_surface(t,y,x)> - <bathymetry(y,x)>'
            msg = msg+'\n but inputs <water_depth_at_soil_surface(t,y,x)>, <bathymetry(y,x)> were not specified correctly. Pass proper arrays and restart the script'
            raise ValueError(msg)

        # inputs are correct, generate waterlevel!
        water_level = np.zeros(water_depth_at_soil_surface.shape)
        for t in xrange(water_depth_at_soil_surface.shape[0]):
            water_level[t, ...] = water_depth_at_soil_surface[t, ...] - bathymetry


            if log:
                msg = 'get_water_level(): Searching for water_level data. Variable <{0}> not found in none of the following files:'.format(varname)
                for nc_file in list_with_filenames:
                    msg = msg+'\nget_water_level():\t'+nc_file
                msg = msg+'\nget_water_level(): Therefore function calculated water_level as <water_level(t,y,x)> = <water_depth_at_soil_surface(t,y,x)> - <bathymetry(y,x)>'
                msg = msg+'\nget_water_level(): Returning array of shape {0}'.format(water_level.shape)
                print msg

    else:
        if log:
            print 'get_water_level(): Found <water_level> of shape {0}, from variable <{1}>)'.format(water_level.shape, varname)
            print 'get_water_level(): From file <{0}>'.format(nc_file_found)
        pass
  
    return water_level
    
