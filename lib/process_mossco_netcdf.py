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
import ui
from Mesh2 import gridhelp




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

    if mask is None and flatten:
        a = a.T.flatten(order='F')

    elif mask is not None:
        a = np.ma.array(a, mask=mask.T)
        if flatten:
            a = a.T.flatten(order='F').compressed()

    nc.close()
    del nc
    if log: print 'read_mossco_nc_2d(): returning {0} of shape {1}'.format(type(a), a.shape)
    return  a


def read_mossco_nc_3d(filename, varname, flatten=False, mask=None, log=False):
    nc   = Dataset(filename, mode='r')
    Vars = nc.variables
    v    = Vars[varname]
    a    = np.squeeze(v[:])
    
    if len(a.shape) != 3:
        raise TypeError('read_mossco_nc_3d(): Invalid array shape. Should be 3D. Received shape {0}, after squeezing {1}'.format(v.shape, a.shape))

    if mask is None and flatten:
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
                if log: print a.shape
                if log: print layer_thicness0.shape
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

    nc.close()
    del nc
    if log: print 'read_mossco_nc_3d: returning {0} of shape {1}'.format(type(a), a.shape)
    return a


def read_mossco_nc_4d(filename, varname, mask=None, log=False):
    # this will work only for 4d variables which have dimensions (time, z, y, x)
    nc = Dataset(filename, mode='r')
    Vars = nc.variables

    v = Vars[varname]
    vv = np.squeeze(v[:])
    if len(v.shape) != 4:
        raise TypeError('read_mossco_nc_4d(): Invalid array shape. Should be 4D. Received shape {0}, after squeezing {1}'.format(Vars[varname].shape, vv.shape))

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

    nc.close()
    del nc
    #return a.flatten(1), dims_v
    if log: print 'read_mossco_nc_4d: returning {0} of shape {1}'.format(type(a), a.shape)
    return a




def read_mossco_nc_rawvar(filename, varname, mask=None):
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
            user_input = raw_input(_n+' Set number of vertical layers manually. To continue with 2D, set 1 (nLayers integer >=1):')
            try:
                nLayers = int(user_input)
            except:
                nLayers = -1

        return nLayers, None




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
        print '-'*50
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
    _n = 'get_sigma_coordinates():'
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
        print _n, 'Found sigma-coords:', sigma, '\n From variable <{0}>)'.format(varname)
        print _n, 'Values represent coordinates of layer <{0}>'.format(sigma_type)
        print _n, 'Number of layers <{0}>; Number of sigma coords found <{1}>'.format(nLayers, sigma.__len__())
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
    
