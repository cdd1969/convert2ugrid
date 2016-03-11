#!/usr/bin/python
#-*- coding: utf-8 -*-
# Copyright (C) 2015, Nikolai Chernikov <nikolai.chernikov.ru@gmail.com>
#
# "convert2ugrid" is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License v3+. "convert2ugrid" is distributed in the
# hope that it will be useful, but WITHOUT ANY WARRANTY.  Consult the file
# LICENSE.GPL or www.gnu.org/licenses/gpl-3.0.txt for the full license terms.


'''
This module contains functions for processing data stored in
MOSSCO output netcdf file
'''

from __future__ import division
from netCDF4 import Dataset
import numpy as np
import sys
from . import ui, sprint
import traceback


def read_mossco_nc_0d(filename, varname, log=False):
    nc = Dataset(filename, mode='r')
    Vars = nc.variables

    v = Vars[varname]
    a = v[...]
    if len(a.shape) != 0:
        raise TypeError('read_mossco_nc_0d(): Invalid array shape. Should be 0D. Received shape {0}, after squeezing {1}'.format(Vars[varname].shape, a.shape))

    nc.close()
    del nc
    if log: print 'read_mossco_nc_0d(): returning {0}'.format(a)
    return a


def read_mossco_nc_1d(filename, varname, log=False):
    nc = Dataset(filename, mode='r')
    Vars = nc.variables

    a = np.squeeze(Vars[varname][:])
    if len(a.shape) != 1:
        raise TypeError('read_mossco_nc_1d(): Invalid array shape. Should be 1D. Received shape {0}, after squeezing {1}'.format(Vars[varname].shape, a.shape))

    nc.close()
    del nc
    if log: print 'read_mossco_nc_1d(): returning {0} of shape {1}'.format(type(a), a.shape)
    return a


def read_mossco_nc_2d(filename, varname, flatten=False, mask=None, log=False):
    nc = Dataset(filename, mode='r')
    Vars = nc.variables

    a = np.squeeze(Vars[varname][:])
    if len(a.shape) != 2:
        raise TypeError('read_mossco_nc_2d(): Invalid array shape. Should be 2D. Received shape {0}, after squeezing {1}'.format(Vars[varname].shape, a.shape))
    if mask is not None and a.shape != mask.shape:
        raise ValueError('read_mossco_nc_2d(): Invalid data or mask shape. Should be equal. Received data shape(squeezed) - {0}, mask - {1}'.format(a.shape, mask.shape))

    if mask is None and flatten:
        a = a.flatten(order='F')

    elif mask is not None:
        a = np.ma.array(a, mask=mask)
        if flatten:
            a = a.flatten(order='F').compressed()

    nc.close()
    del nc
    if log: print 'read_mossco_nc_2d(): returning {0} of shape {1}'.format(type(a), a.shape)
    return  a


def read_mossco_nc_3d(filename, varname, flatten=False, mask=None, log=False):
    _n   = 'read_mossco_nc_3d(): '
    nc   = Dataset(filename, mode='r')
    Vars = nc.variables
    v    = Vars[varname]
    vs   = np.squeeze(v[:])
    
    if not flatten:
        a = vs
    
    if len(vs.shape) != 3:
        raise TypeError(_n+'Invalid array shape. Should be 3D. Received shape {0}, after squeezing {1}'.format(v.shape, vs.shape))
    if mask is not None and vs[0, :, :].shape != mask.shape:
        raise ValueError(_n+'Invalid data or mask shape. Last two dimensions should be equal. Received data shape(squeezed) - {0}, mask - {1}'.format(vs.shape, mask.shape))

    if mask is None and flatten:
        a = np.zeros(tuple([vs.shape[0], vs.shape[1]*vs.shape[2]]))
        for i in xrange(v.shape[0]):
            a[i, :] = vs[i, ...].flatten(1)

    elif mask is not None:  # IF MASKED ARRAY
        n_valid_2d = np.sum(np.invert(mask))  #number of valid elements in 2d grid. invert - because True(True=1) is an invalid element
        #print _n+'working with masked array. Mask shape {0}. N_valid elements {1}'.format(mask.shape, n_valid_2d)
        #print _n+'data shape {0}'.format(v.shape)
        if flatten:
            a = np.zeros(tuple([vs.shape[0], n_valid_2d]))
        for i in xrange(vs.shape[0]):
            var_masked = np.ma.array(vs[i, ...], mask=mask)
            if flatten:
                var_masked = var_masked.flatten(order='F').compressed()
                #print _n+'var_masked shape {0}'.format(var_masked.shape)
            a[i, ...] = var_masked
    nc.close()
    del nc
    if log: print _n+'returning {0} of shape {1}'.format(type(a), a.shape)
    return a


def read_mossco_nc_4d(filename, varname, flatten=False, mask=None, log=False):
    # this will work only for 4d variables which have dimensions (time, z, y, x)
    _n   = 'read_mossco_nc_4d(): '
    nc   = Dataset(filename, mode='r')
    Vars = nc.variables
    v    = Vars[varname]
    vs   = np.squeeze(v[:])
    if not flatten:
        a = vs

    if len(vs.shape) != 4:
        raise TypeError(_n+'Invalid array shape. Should be 4D. Received shape {0}, after squeezing {1}'.format(v.shape, vs.shape))
    if mask is not None and vs[0, 0, :, :].shape != mask.shape:
        raise ValueError(_n+'Invalid data or mask shape. Last two dimensions should be equal. Received data shape(squeezed) - {0}, mask - {1}'.format(vs.shape, mask.shape))

    if mask is None and flatten:  # IF NOT MASKED ARRAY
        a = np.zeros(tuple([vs.shape[0], vs.shape[1], vs.shape[2]*vs.shape[3]]))
        for t in xrange(vs.shape[0]):
            for z in xrange(vs.shape[1]):
                a[t, z, :] = vs[t, z, ...].flatten(1)


    elif mask is not None:  # IF MASKED ARRAY
        n_valid_2d = np.sum(np.invert(mask))  #number of valid elements in 2d part. invert - because True is an invalid element
        if flatten:
            a = np.zeros(tuple([v.shape[0], v.shape[1], n_valid_2d]))
        for t in xrange(v.shape[0]):
            for z in xrange(v.shape[1]):
                var_masked = np.ma.array(v[t, z, ...], mask=mask)
                if flatten:
                    var_masked = var_masked.flatten(order='F').compressed()
                a[t, z, :] = var_masked

    nc.close()
    del nc
    if log: print _n+'returning {0} of shape {1}'.format(type(a), a.shape)
    return a




def read_mossco_nc_rawvar(filename, varname):
    nc = Dataset(filename, mode='r')
    a = nc.variables[varname][:]
    nc.close()
    return a



def get_number_of_depth_layer_from_mossco(list_with_filenames, dimname='getmGrid3D_getm_3'):
    _n = 'get_number_of_depth_layer_from_mossco():'
    nLayers = False
    for nc_file in list_with_filenames:
        try:
            nc = Dataset(nc_file, mode='r')
            d_nc = nc.dimensions[dimname]
            nLayers = d_nc.__len__()
            nc.close()
            layer_fname = nc_file
        except KeyError:  # if dimension is not found in cuurent file > skip to next file
            pass
        except RuntimeError as err:
            print _n, 'filelist:', list_with_filenames
            print _n, 'reading file:', nc_file
            raise err
    if nLayers:
        print _n, 'found vertical-layers:', nLayers, '\t(from dimension <{0}>) in file <{1}>'.format(dimname, layer_fname)
        return nLayers, layer_fname
    else:
        print _n, 'Vertical layers not recognized. Dimension <{0}> not found in files: {1}'.format(dimname, list_with_filenames)
        nLayers = -1
        while (nLayers < 1):
            user_input = ui.promt(_n+' Set number of vertical layers manually. To continue with 2D, set 1 (nLayers integer >=1):', color='yellow', type=int, show_default=False)
            try:
                nLayers = int(user_input)
            except:
                nLayers = -1

        return nLayers, None




def get_davit_friendly_variables(filename, tdim=['time'], zdim=['getmGrid3D_getm_3'],
        ydim=['getmGrid3D_getm_2', 'getmGrid2D_getm_2', 'y', 'yx', 'y_x', 'yX', 'y_X', 'yt', 'y_t', 'y_T', 'yT' 'yc', 'y_c', 'y_C',
              'lat', 'latc', 'latx', 'latt', 'lat_c', 'lat_x', 'lat_t'],
        xdim=['getmGrid3D_getm_1', 'getmGrid2D_getm_1', 'x', 'xx', 'x_x', 'xX', 'x_X', 'xt', 'x_t', 'x_T', 'xT' 'xc', 'x_c', 'x_C',
              'lon', 'lonc', 'lonx', 'lont', 'lon_c', 'lon_x', 'lon_t'],
        log=True):
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
    _n = 'get_davit_friendly_variables(): '
    FRIENDLY_VAR_DICT = dict()
    _1D = list()
    _2D = list()
    _3D = list()
    _4D = list()

    if log: print _n, 'Scanning file :', filename
    root_grp = Dataset(filename, mode='r')
    for var_name, var in root_grp.variables.iteritems():
        dims = [str(d) for d in var.dimensions]
        sprint('Found: ', log=log, newline=False, indent='\t')
        sprint(var_name+' ', log=log, newline=False, mode='warning')
        sprint(str(dims)+' >>> ', log=log, newline=False)
        
        if (len(var.shape) == 1) and (dims[0] in tdim):
            _1D.append(var_name)
            sprint('added', log=log, mode='ok')
        if (len(var.shape) == 2) and (dims[0] in ydim) and (dims[1] in xdim):
            _2D.append(var_name)
            sprint('added', log=log, mode='ok')
        elif (len(var.shape) == 3) and (dims[0] in tdim) and (dims[1] in ydim) and (dims[2] in xdim):
            _3D.append(var_name)
            sprint('added', log=log, mode='ok')
        elif (len(var.shape) == 4) and (dims[0] in tdim) and (dims[1] in zdim) and (dims[2] in ydim) and (dims[3] in xdim):
            _4D.append(var_name)
            sprint('added', log=log, mode='ok')
        else:
            sprint('skipped', log=log, mode='fail')

    root_grp.close()
    FRIENDLY_VAR_DICT['1D'] = _1D
    FRIENDLY_VAR_DICT['2D'] = _2D
    FRIENDLY_VAR_DICT['3D'] = _3D
    FRIENDLY_VAR_DICT['4D'] = _4D
    
    return FRIENDLY_VAR_DICT


def get_water_level(list_with_filenames, wl_vname='water_level', water_depth_at_soil_surface_vname=None, bathymetry_vname=None, log=False, indent=''):
    '''Try to get water_level information from variable `wl_vname`.
    If fails, try to generate it as:
        water_level = `water_depth_at_soil_surface` - `bathymetry`

    Args:
    -----
        list_with_filenames (list(str)):
            list with names of netcdf files to work with

        wl_vname (str):
            name of the variable with water-level data. Water level is with respect to MSL.
            Default: 'water_level'

        water_depth_at_soil_surface_vname (None or str):
            IS USED ONLY IF VARIABLE <wl_vname> NOT FOUND
            None or name of the variable that represents water depth at soil surface.
            Distance between sea bed and water surface. Always positive
        
        bathymetry_vname (None or str):
            IS USED ONLY IF VARIABLE <wl_vname> NOT FOUND
            None or name of the variable that represents bathymetry values. Distance between
            sea bed and MSL. Always positive.
    
    Return:
    -------
        water_level (3D-array):
            3D-array of floats of (time, y, x) shape, representing water level with respect to MSL
    '''
    _n = 'get_water_level():'
    _i = indent+'\t'
    water_level = None
    
    for nc_file in list_with_filenames:
        root_grp = Dataset(nc_file, mode='r')
        if wl_vname in root_grp.variables.keys():
            nc_file_found = nc_file
            water_level = root_grp.variables[wl_vname][:]
            root_grp.close()
            break
    

    if not water_level:  # if <water_level> has not been found... generate it!
        sprint(_n, ' Searching for water_level data. Variable <{0}> not found in none of files: {1}\n'.format(wl_vname, list_with_filenames), log=log, indent=_i)
        sprint(_n, ' Now try to calculate water_level as <water_level(t,y,x)> = <water_depth_at_soil_surface(t,y,x)> - <bathymetry(y,x)>.\n {2}\tWith\n{2}\twater_depth_at_soil_surface={0}\n{2}\tbathymetry={1}'.format(water_depth_at_soil_surface_vname, bathymetry_vname, _i), log=log, indent=_i)

        # First find arrays based on varnames
        water_depth_at_soil_surface = None
        bathymetry = None
        wd_found = False
        ba_found = False
        for nc_file in list_with_filenames:
            root_grp = Dataset(nc_file, mode='r')

            if water_depth_at_soil_surface_vname in root_grp.variables.keys() and wd_found is False:
                water_depth_at_soil_surface = read_mossco_nc_rawvar(nc_file, water_depth_at_soil_surface_vname)
                wd_found = 'found in file {0}'.format(nc_file)

            if bathymetry_vname in root_grp.variables.keys() and ba_found is False:
                bathymetry = read_mossco_nc_rawvar(nc_file, bathymetry_vname)
                ba_found = 'found in file {0}'.format(nc_file)

            root_grp.close()
            
            if ba_found is not False and wd_found is not False:
                break

        if water_depth_at_soil_surface is None or bathymetry is None:
            raise ValueError('Cannot calculate `water_level`. When water-level is not found, in order to calculate it, both `water_depth_at_soil_surface` and `bathymetry` must be present. Variable for `water_depth_at_soil_surface`={0} {1}. Variable for `bathymetry`={2} {3}'.format(water_depth_at_soil_surface_vname, wd_found if wd_found else 'not found', bathymetry_vname, ba_found if ba_found else 'not found'))

        # arrays are found, generate waterlevel!
        water_level = np.zeros(water_depth_at_soil_surface.shape)
        for t in xrange(water_depth_at_soil_surface.shape[0]):
            water_level[t, ...] = water_depth_at_soil_surface[t, ...] - bathymetry
        sprint(_n, ' Variable <water_level> of shape {0}, calculated successfully based on : {1})'.format(water_level.shape, [water_depth_at_soil_surface_vname, bathymetry_vname]), log=log, indent=_i)

    else:
        sprint(_n+' Found <water_level> of shape {0}, from variable <{1}> from file {2})'.format(water_level.shape, wl_vname, nc_file_found), log=log, indent=_i)

    return water_level


def get_sigma_coordinates(list_with_filenames, nLayers, sigma_varname='level', waterdepth_varname='water_depth_at_soil_surface',
        layerdepth_varname='getmGrid3D_getm_layer', log=False, indent=''):
    ''' Get sigma coordinates of layer center/layer borders. This is done in two ways:

        1) Try to read specific variable given with `sigma_varname`

        2) If 1st approach fails, try to calculate sigma coordinates based on
            `waterdepth_varname` and `layerdepth_varname` (to get more info about
            these variable read docstring of function `caclulate_relative_layer_thickness()`)

    Args:
    -----
        list_with_filenames (list(str)):
            list with names of netcdf files to work with
        
        nLayers (int):
            number of vertical layers in the file
        
        sigma_varname (str):
            ONLY NEEDED FOR APPROACH #1
            name of the variable (in one of the netcdf files) with sigma layer information.
            The data stored within should be a 1D-array of <nLayers> (layer-centers) or
            <nLayers>+1 (layer-borders) length.

        waterdepth_varname (str):
            ONLY NEEDED FOR APPROACH #2
            name of the variable that represents water depth at soil surface.
            Units must be shared with `layerdepth_varname`. Always positive.
            3D-Array with (time, y, x) dimensions;
            Default: "water_depth_at_soil_surface"
        
        layerdepth_varname (str):
            ONLY NEEDED FOR APPROACH #2
            name of the variable in `nc_in` netcdf file, that represents
            layer depth below water surface at the element center. Units
            must be shared with `waterdepth_varname`. Always negative.
            3D-Array with (z, y, x) dimensions;
            Default: "getmGrid3D_getm_layer"

        log (bool):
            flag to print additional output
    
    Return:
    -------
        sigma (1D-array):
            1D array that represents sigma coordinates at cell centers/borders. Therefore it can
            has length of <nLayers> or <nLayers+1>.
        sigma_type ('center'|'border'):
            position of sigma coordinates: 'center' or 'border'
    
    Example:
    --------
        sigma = [-1., -0.75, -0.5, -0.25, 0.], with nLayers=4 it means that we have 4 layers 25% each,
                    and the values describes the border coords. In this case output will look like...
        Example Output:
            [[-1., -0.75, -0.5, -0.25, 0.], 'border']
            [[-0.875, -0.625, -0.375, -0.125], 'center']

    '''
    _i = indent+'\t'
    _n = 'get_sigma_coordinates():'
    sprint (_n, 'Try to get sigma coordinates', log=log, indent=indent, mode='bold')
    sigma_found = False
    sigma_type = None
    
    for nc_file in list_with_filenames:
        root_grp = Dataset(nc_file, mode='r')
        if sigma_varname in root_grp.variables.keys():
            var_nc = root_grp.variables[sigma_varname]
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
        sprint('Found sigma-coords:', sigma, '\n From variable <{0}>)'.format(sigma_varname), log=log, indent=_i)
        sprint('Values represent coordinates of layer <{0}>'.format(sigma_type), log=log, indent=_i)
        sprint('Number of layers <{0}>; Number of sigma coords found <{1}>'.format(nLayers, sigma.__len__()), log=log, indent=_i)
        return [sigma, sigma_type]
    else:
        kwargs = dict()
        kwargs['waterdepth_varname'] = waterdepth_varname
        kwargs['layerdepth_varname'] = layerdepth_varname

        sprint('Sigma coordinates not found. Variable <{0}> not found in files: {1}'.format(sigma_varname, list_with_filenames), log=log, indent=_i, mode='warning')
        ui.promt('{0}{2} Now i will try calculating relative layer thickness based on variables:\n{0}{1}\nPress Enter to continue'.format(indent, kwargs, _n), color='yellow', show_default=False)
        sigma = None
        for nc_file in list_with_filenames:
            try:
                sigma = get_mossco_relative_layer_thickness(nc_file, indent=_i+'\t', **kwargs)
            except Exception, err:
                #sprint(err, mode='fail')
                #traceback.print_exc()
                continue
            if sigma is not None:
                sigma_found = True
                sigma_type = 'center'
                break

        if sigma_found:
            sprint('Calculated sigma-coords:', sigma, '\n{1}From variables:\n{1}{0})'.format(kwargs, _i), log=log, indent=_i)
            sprint('Values represent coordinates of layer <{0}>'.format(sigma_type), log=log, indent=_i)
            sprint('Number of layers <{0}>; Number of sigma coords <{1}>'.format(nLayers, sigma.__len__()), log=log, indent=_i)
            return [sigma, sigma_type]
        else:
            raise ValueError('Sigma coordinates neither found (1) nor calculated (2).\n1): Variable <{0}> not found in files: {1}\n\n2): Sigma was not calculated based on variables {2}'.format(sigma_varname, list_with_filenames, kwargs))


def get_mossco_relative_layer_thickness(nc_in,
        waterdepth_varname='water_depth_at_soil_surface',
        layerdepth_varname='getmGrid3D_getm_layer',
        indent=''):
    ''' Function return 1D array with relative layer thickness (sigma-layers) of cell centers.
        See doc string of function `caclulate_relative_layer_thickness()`

    Return:
    -------
        rel_thick_1D (1D-array):
            1D array of sigma coordinates (either cell center of cell borders)
    '''
    nc = Dataset(nc_in , mode='r')
    
    #>>> Now get the relative layer thickness
    layer_relthickness = caclulate_relative_layer_thickness(nc.variables[layerdepth_varname][:], nc.variables[waterdepth_varname][:], include_time=False)

    valid_cell_ji = (layer_relthickness[0, :, :].nonzero()[0][0], layer_relthickness[0, :, :].nonzero()[1][0])  # see ISSUE #2
    rel_thick_1D = layer_relthickness[:, valid_cell_ji[0], valid_cell_ji[1]]  # see ISSUE #1 , #2
    print indent+'>>> layerdepth  ', nc.variables[layerdepth_varname][:, valid_cell_ji[0], valid_cell_ji[1]]  # see ISSUE #1 , #2
    print indent+'>>> waterdepth  ', nc.variables[waterdepth_varname][0, valid_cell_ji[0], valid_cell_ji[1]]  # see ISSUE #1 , #2
    print indent+'>>> relthickness', rel_thick_1D
    nc.close()
    return rel_thick_1D
    

def caclulate_relative_layer_thickness(layer_depth, water_depth, include_time=False):
    '''Calculate relative layer thickness for mutlidimensional array
    (t, z, y, x) at timestep t=0 at cell center
    Assumptions:
        1) z=0 is the closest to soil-surface layer, where `z` is the index
           of the vertical dimension of the matrix
        2) relative layer thickness is not time dependent !!!
    Args:
    -----
        layer_depth (3D numpy array):
            3D-Array with (z, y, x) dimensions; represents layer depth below
            water surface (0, level) at center of each z,y,x cell.
            Units are shared with `water_depth` array. Always negative.
        
        water_depth (3D numpy array):
            3D-Array with (time, y, x) dimensions; represents water depth at soil
            surface. Units are shared with `layer_depth` array. Always positive.
        include_time (Optional[bool]):
            Flag to control the shape of the output array. It does not affect the data
            that is stored in output array: it is completely identical.
            ---
            If `True` - the shape will be (time, z, y, x)
            If `False` - the shape will be (z, y, x)  (DEFAULT)
    Return:
    -------
        layer_relthickness (3D|4D numpy array):
            3D|4D-Array with (z, y, x)|(t, z, y, x) dimensions; represents relative layer thickness.
            Is dimensionless and is constant within simulation period (time). Values are at cell center
    '''
    # >>> Get dimensions
    z, y, x = layer_depth.shape
    t, y, x = water_depth.shape

    # >>> Get optional mask
    if isinstance(layer_depth, np.ma.MaskedArray):
        mask = layer_depth.mask
        # >>> Allocate memory for `layear_thickness` array, initialize it. This array represents layer thickness at timestep t=0. Values are always positive
        layer_thickness = np.ma.array(np.empty((z, y, x), dtype=float), mask=mask)
        # >>> Allocate memory for `relative_thcikness` array, initialize it. This array represents relaitve layer thickness at timestep t=0 with respect to total water-depth. Values are always positive, dimensionless
        layer_relthickness = np.ma.array(np.empty((z, y, x), dtype=float), mask=mask)
    elif isinstance(water_depth, np.ma.MaskedArray):
        # add new axis at 0-index position
        # and repeat the array (z,y,x) `t` times along 0-index position axis
        mask = water_depth[0, :, :].mask
        mask = mask.reshape(1, y, x).repeat(z, 0)
        layer_thickness = np.ma.array(np.empty((z, y, x), dtype=float), mask=mask)
        layer_relthickness = np.ma.array(np.empty((z, y, x), dtype=float), mask=mask)
    else:
        layer_thickness = np.empty((z, y, x), dtype=float)
        layer_relthickness = np.empty((z, y, x), dtype=float)
    


    # >>> Layer thickness of the near-bottom layer
    layer_thickness[0, :, :] = (water_depth[0, :, :] - (-layer_depth[0, :, :]) ) * 2.0
    # >>> Layer thickness of the rest layers
    for k in xrange(1, z):
        #layer_thickness[k, :, :] = (water_depth[0, :, :] - (-layer_depth[k, :, :]) - layer_thickness[k-1, :, :]) * 2.0
        layer_thickness[k, :, :] = (water_depth[0, :, :] - (-layer_depth[k, :, :]) - np.sum(layer_thickness[0:k, :, :], axis=0)) * 2.0
    # >>> Now calculate relative layer thickness
    for k in xrange(z):
        layer_relthickness[k, :, :] = abs(layer_thickness[k, :, :] / water_depth[0, :, :])
    # >>> Finally return the result
    if include_time is False:
        return layer_relthickness
    else:
        # add new axis at 0-index position
        # and repeat the array (z,y,x) `t` times along 0-index position axis
        return layer_relthickness.reshape(1, z, y, x).repeat(t, 0)
