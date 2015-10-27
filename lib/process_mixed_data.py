#!/usr/bin/python
#-*- coding: utf-8 -*-
# Copyright (C) 2015, Nikolai Chernikov <nikolai.chernikov.ru@gmail.com>
#
# "convert2ugrid" is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License v3+. "convert2ugrid" is distributed in the
# hope that it will be useful, but WITHOUT ANY WARRANTY.  Consult the file
# LICENSE.GPL or www.gnu.org/licenses/gpl-3.0.txt for the full license terms.

import numpy as np
import os.path
import netCDF4

import ui
import process_cdl
import process_mossco_netcdf
from Mesh2 import gridhelp


class netcdfVariableReader(object):
    """ This is the basic class for inheritance of other classes,
        providing methods to access netCDF files"""
    
    def __init__(self):
        pass
    
    def check_dtype(self, _object, dtype, raise_error=True):
        """Compares type(_object) to "dtype" """
        if isinstance(_object, dtype):
            return True
        else:
            # here we treat unicode and simple string as same objects
            if isinstance(_object, unicode) and dtype is str:
                return True
            if raise_error: raise TypeError('<{0}> should be of type <{2}>. Is {1}'.format(_object, type(_object), str(dtype)))
            else: return False

    def get_string_with_netcdf_varnames(self, ncname, prefix='', withdims=True):
        """ Method returns a nice string of enumerated varnames within netcdf file
        Inputs:
            ncname      -   string, name of the netcdf datafile
            prefix      -   string, prefix for first line of the "s"
        Returns:
            s           -   string
        """
        self.check_dtype(ncname,   str)
        self.check_dtype(prefix,   str)
        self.check_dtype(withdims, bool)
        
        nc = netCDF4.Dataset(ncname, mode='r')
        s = prefix+' File <{0}> contains following variables:'.format(ncname)
        for i, varname in enumerate( sorted(nc.variables.keys()) ):
            if withdims: s += '\n\t{0:3d}) {1} {2} {3}'.format(i, varname, nc.variables[varname].dimensions, nc.variables[varname].shape)
            else: s += '\n\t{0:3d}) {1}'.format(i, varname)
        nc.close()
        del nc
        return s


    def variableIsFound(self, ncname, varname):
        """ Method tells user if variable has been found in netcdf file.
        Inputs:
            ncname      -   string, name of the netcdf datafile
            varname     -   string, name of the variable within netcdf file
        Returns:
            found       -   True/False
        """
        if varname is None: return False
        self.check_dtype(ncname,  str)
        self.check_dtype(varname, str)

        found = False
        
        nc = netCDF4.Dataset(ncname, mode='r')
        if varname in nc.variables.keys():
            found = True
        nc.close()
        del nc
        return found


    def promtVariableNameInput(self, ncname, promtInputMessage='Type the name of the variable to pick:', printAvailableVarNames=True):
        """ Method asks user to choose variable from file
        Inputs:
            ncname      -   string, name of the netcdf datafile
            promtInputMessage - string, message to ask for variable input. Is used only if "promtInput=True"
        Returns:
            varname     -   string, name of the variable within netcdf file, since it may not be equall to input
        """
        self.check_dtype(printAvailableVarNames, bool)
        self.check_dtype(promtInputMessage, str)
        self.check_dtype(ncname, str)
        
        nc = netCDF4.Dataset(ncname, mode='r')
        if printAvailableVarNames:
            print self.get_string_with_netcdf_varnames(ncname)
            print 'Select variable from file <{0}>. All available variables are listed above'.format(ncname)
        else:
            pass
            #print 'Select variable from file <{0}>'.format(ncname)
        varname = None
        while varname not in nc.variables.keys():
            varname = raw_input(promtInputMessage)
        nc.close()
        del nc
        return varname


    def read_netcdfVarMetadata(self, ncname, varname, raise_error=True, promtInput=True, **kwargs):
        """ Method reads given variable from netcdf file and returns metadata
        Inputs:
            ncname      -   string, name of the netcdf datafile
            varname     -   string, name of the variable within netcdf file
            raise_error -   True/False, if True will raise error when file not found
            promtInput  -   True/False,If True, will not return False if variable not found
                            not found in netcdf file. Instead ask user to type var-name manually,
                            using "promtInputMessage" in **kwargs
            **kwargs    -   will be passed to promtVariableNameInput()
        Returns:
            metadata    -   dictionary or None
        """
        self.check_dtype(ncname,  str)
        self.check_dtype(varname, str)

        if not self.variableIsFound(ncname, varname):
            if promtInput is True:
                varname = self.promtVariableNameInput(**kwargs)
            else:
                if raise_error:
                    raise IndexError('No such variable <0> in file <{1}>'.format(varname, ncname))
                else:
                    return None

        nc = netCDF4.Dataset(ncname, mode='r')
        metadata = dict()
        var = nc.variables[varname]
        metadata['fname'] = ncname     # source
        metadata['name']  = varname    # source
        metadata['dims']  = var.dimensions
        metadata['dtype'] = var.dtype
        metadata['shape'] = var.shape  # will be <()> if single value
        metadata['size']  = var.size   # number of elements
        # for some reason not working, giving "attribute 'mask' not found"
        #metadata['mask']  = var.mask   # True/False flag
        metadata['fillvalue']    = var._FillValue if '_FillValue' in var.ncattrs() else None
        metadata['nNonOneDims'] = len(filter(lambda a: a != 1, var.shape))  # length of the array.shape excluding single dimensions, for <()> will give 0
        metadata['attrs'] = dict()
        if var.ncattrs():  # if not empty list
            for attr_n in sorted(var.ncattrs()):
                metadata['attrs'][attr_n] = var.getncattr(attr_n)
        nc.close()
        del nc
        return metadata


    def read_netcdfVarData(self, ncname, varname, squeeze=False, raise_error=True, promtInput=True, **kwargs):
        """ Method reads given variable from netcdf file and returns array
        Inputs:
            ncname      -   string, name of the netcdf datafile
            varname     -   string, name of the variable within netcdf file
            squeeze     -   True/False. If True, remove single-dimensional entries from the shape of an array.
            raise_error -   True/False, if True will raise error when file not found
            promtInput  -   True/False,If True, will not return False if variable not found
                            not found in netcdf file. Instead ask user to type var-name manually,
                            using "promtInputMessage" in **kwargs
            **kwargs    -   will be passed to variableIsFound()
        Returns:
            data        -   ndarray or None
        """
        self.check_dtype(ncname,  str)
        self.check_dtype(varname, str)
        self.check_dtype(squeeze, bool)

        if not self.variableIsFound(ncname, varname):
            if promtInput is True:
                varname = self.promtVariableNameInput(**kwargs)
            else:
                if raise_error:
                    raise IndexError('No such variable <0> in file <{1}>'.format(varname, ncname))
                else:
                    return None
        nc = netCDF4.Dataset(ncname, mode='r')
        data = nc.variables[varname][...]  # with this syntax we catch also null-dimensional arrays (i.e. scalars)
        if squeeze:
            # note here errors may arise: if scalar is used, it will be converted to ndarray of shape <()>
            data = np.squeeze(data)
        nc.close()
        del nc
        return data








class cdlVariableExt(netcdfVariableReader):
    """this class is an extension of class cdlVariable()
        It is considered as an addon to already existing metadata variable.
        
        Inputs:
            parent  - instance of type "cdlVariable".

        This object, in contrast to its parent, can process netcdf data-file and read
        meta-data from them. This is exactly hat is being done: gathering meta-info from
        the source netcdf-variable of the passed parent"""
    
    def __init__(self, parent):
        super(cdlVariableExt, self).__init__()
        if not isinstance(parent, process_cdl.cdlVariable):
            raise TypeError('Invalid type. Should be <cdlVariable>, received {0}'.format(type(parent)))
        self.__parent = parent

        self._init_constants()

        # do not know if it is good idea to check at the init stage
        # whether netcdf file is valid, and var is whithin
        self.__source = self._search_source_metadata()


    def _init_constants(self):
        # overwriting method of parent
        self.__source_attrs = ['_source_filename', '_source_varname']


    def _search_source_metadata(self):
        parent_attrs = self.__parent.get_attrs()
        if all(attr in parent_attrs.keys() for attr in self.__source_attrs):
            # if the file exists
            if os.path.isfile( parent_attrs[self.__source_attrs[0]] ):
                return self.read_netcdfVarMetadata(parent_attrs[self.__source_attrs[0]], parent_attrs[self.__source_attrs[1]])

            else:
                raise IOError('Special attribute <_source_filename>: No such file <{0}>'.format(parent_attrs[self.__source_attrs[0]]))

        else:
            raise ValueError('Add missing attributes to know where to search for data. Both these attributess should exist {0}'.format(self.__source_attrs))


    def get_parent(self):
        return self.__parent

    def get_source_metadata(self):
        #return self._search_source_metadata()
        return self.__source


















class containerForGridGeneration(netcdfVariableReader):
    """this class is a container for information, required to generate
        grid."""
    
    def __init__(self, topofile, log=False):
        super(containerForGridGeneration, self).__init__()
        self.check_dtype(topofile, str)
        self.__topofile = topofile
        self._init_constants()

        self._init_coordinates_and_bathymetry_metadata(self.__topofile, log=log)
        if log: print self.get_string_metadata_summary()

    def _init_constants(self):
        self.__bathymetry = dict()
        self.__x_coords   = dict()
        self.__y_coords   = dict()
        self.__grid_info  = dict()

        self.__constants = dict()
        self.__constants['x_cartesian_vnames'] = ['x', 'xx', 'x_x' , 'xX', 'x_X', 'xt', 'x_t', 'xT', 'x_T']
        self.__constants['x_geographi_vnames'] = ['lon', 'lonx', 'lon_x' , 'lonX', 'lon_X', 'lont', 'lon_t', 'lonT', 'lon_T']
        self.__constants['y_cartesian_vnames'] = ['y', 'yx', 'y_x' , 'yX', 'y_X', 'yt', 'y_t', 'yT', 'y_T']
        self.__constants['y_geographi_vnames'] = ['lat', 'latx', 'lat_x' , 'latX', 'lat_X', 'latt', 'lat_t', 'latT', 'lat_T']
        self.__constants['bathymetry_vnames']  = ['bathymetry']
        self.__constants['fv_attr_namelist']   = ['_FillValue', 'missing_value']



    def get_constants(self):
        return self.__constants

    def get_metadata(self):
        return {'x': self.__x_coords['meta'],
                'y': self.__y_coords['meta'],
                'b': self.__bathymetry['meta'],
                'grid_info': self.__grid_info}
    
    def get_data(self, **kwargs):
        return {'x': self.read_netcdfVarData(self.__topofile, self.get_metadata()['x']['name'], **kwargs),
                'y': self.read_netcdfVarData(self.__topofile, self.get_metadata()['y']['name'], **kwargs),
                'b': self.read_netcdfVarData(self.__topofile, self.get_metadata()['b']['name'], **kwargs)}

    def get_mask(self, transpose=True, **kwargs):
        return self._generate_mask(transpose=transpose, **kwargs)

    def get_string_metadata_summary(self):
        s = ''
        s += '\n'+'-------------------------------------------------------------------------------'
        s += '\n'+'------------------------- Grid Detection Summary ------------------------------'
        s += '\n'+'-------------------------------------------------------------------------------'
        s += '\n'+'X-coords  :'
        s += '\n'+'            variable   <{0}>'.format(self.get_metadata()['x']['name'])
        s += '\n'+'            dimensions <{0}>'.format(self.get_metadata()['x']['dims'])
        s += '\n'+'            shape      <{0}>'.format(self.get_metadata()['x']['shape'])
        s += '\n'+'            location   <{0}>'.format(self.get_metadata()['x']['points_location'])
        s += '\n'+'            type       <{0}>'.format(self.get_metadata()['x']['coordinate_type'])
        s += '\n'+'            units      <{0}>'.format(self.get_metadata()['x']['attrs']['units'] if 'units' in self.get_metadata()['x']['attrs'] else 'unknown')
        s += '\n'+'Y-coords  :'
        s += '\n'+'            variable   <{0}>'.format(self.get_metadata()['y']['name'])
        s += '\n'+'            dimensions <{0}>'.format(self.get_metadata()['y']['dims'])
        s += '\n'+'            shape      <{0}>'.format(self.get_metadata()['y']['shape'])
        s += '\n'+'            location   <{0}>'.format(self.get_metadata()['y']['points_location'])
        s += '\n'+'            type       <{0}>'.format(self.get_metadata()['y']['coordinate_type'])
        s += '\n'+'            units      <{0}>'.format(self.get_metadata()['y']['attrs']['units'] if 'units' in self.get_metadata()['y']['attrs'] else 'unknown')
        s += '\n'+'Bathymetry:'
        s += '\n'+'            variable   <{0}>'.format(self.get_metadata()['b']['name'])
        s += '\n'+'            dimensions <{0}>'.format(self.get_metadata()['b']['dims'])
        s += '\n'+'            shape      <{0}>'.format(self.get_metadata()['b']['shape'])
        s += '\n'+'            location   <{0}>'.format(self.get_metadata()['b']['points_location'])
        s += '\n'+'Grid Info :'
        s += '\n'+'            type       <{0}>'.format(self.get_metadata()['grid_info']['type'])
        s += '\n'+'-------------------------------------------------------------------------------'
        s += '\n'+'-------------------------------------------------------------------------------'
        s += '\n'+'-------------------------------------------------------------------------------'
        return s

    def _init_coordinates_and_bathymetry_metadata(self, ncname, bathymetry_vname=None, x_vname=None, y_vname=None, log=False):
        self.check_dtype(ncname,  str)
        long_var_list_not_shot = True
        if bathymetry_vname is not None: self.check_dtype(bathymetry_vname, str)
        if x_vname is not None: self.check_dtype(x_vname, str)
        if y_vname is not None: self.check_dtype(y_vname, str)


        # searching for bathymetry
        # ----------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------
        if bathymetry_vname and not self.variableIsFound(ncname, bathymetry_vname) :
            raw_input('Bathymetry: User-defined bathymetry-variable name <{0}> not found in file <{1}>. I will try to match default names. Press Enter to continue'.format(
                            bathymetry_vname, ncname))
        
        if not bathymetry_vname:
            if log: print 'Bathymetry: Will try to match bathymetry-name from default list...'
            for bathymetry_vname in self.get_constants()['bathymetry_vnames']:
                found = self.variableIsFound(ncname, bathymetry_vname)
                if found:  # if var is found
                    if log: print 'Bathymetry: Variable found: <{0}>'.format(bathymetry_vname)
                    break
            if not found:  # we cycled through whole loop and havent found any var
                print 'Bathymetry: Bathymetry not found: No variable with name from default namelist found.'
                print 'Bathymetry: Default bathymetry namelist: {0}'.format(self.get_constants()['bathymetry_vnames'])
                if long_var_list_not_shot: raw_input('Bathymetry: Choose manually. Press Enter to list all available variables in current file')
                ndim = -1
                while ndim != 2:
                    bathymetry_vname = self.promtVariableNameInput(ncname, promtInputMessage='Bathymetry: Type the name of the bathymetry-variable to pick (should be 2D):', printAvailableVarNames=long_var_list_not_shot)
                    ndim = self.read_netcdfVarMetadata(ncname, bathymetry_vname, raise_error=False, promtInput=False)['nNonOneDims']
                    if long_var_list_not_shot: long_var_list_not_shot ^= long_var_list_not_shot

        # now saving bathymetry meta-data
        self.__bathymetry['meta'] = self.read_netcdfVarMetadata(ncname, bathymetry_vname, raise_error=True, promtInput=False)
        if log: print 'Bathymetry: Picking bathymetry variable: <{0}> , dimensions {1}, shape {2}'.format(
                        self.__bathymetry['meta']['name'], self.__bathymetry['meta']['dims'], self.__bathymetry['meta']['shape'])
        # ----------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------

        
        # now picking X and Y coords based on bathymetry
        # ----------------------------------------------
        # 1) try to use passed if any
        #       if not passed
        #           go to (2)
        #       elif not found
        #           go to (2)
        # 2) try to use from bathymetry dim-name
        #       if not found
        #           go to (3)
        # 3) try to use from default list
        #       if nothing has been found
        #           ask user to type manually
        #       elif one var has been found for x and for y
        #           take them
        #       elif more then one var has been found for x and for y
        #           ask user to choose one
        # ----------------------------------------------------
        
        # X coords...
        # ----------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------
        if x_vname and not self.variableIsFound(ncname, x_vname):
                raw_input('X coords: User-defined X-coord-variable name <{0}> not found in file <{1}>. Press Enter to continue auto-search'.format(x_vname, ncname))

        if not x_vname:
            dim_bath = self.__bathymetry['meta']['dims']
            if log: print 'X coords: Trying to find variable based on bathymetry dimensions <{0}>. Searching for variable <{1}>'.format(dim_bath, dim_bath[-1] )
            
            if not self.variableIsFound(ncname, dim_bath[-1]):
                if log: print 'X coords: variable name <{0}> not found in file <{1}>.'.format(dim_bath[-1], ncname)
                x_default_list = self.get_constants()['x_cartesian_vnames']+self.get_constants()['x_geographi_vnames']
                if log: print 'X coords: Searching for variables with name from default list {0}'.format(x_default_list)
                
                x_found = list()
                for x_n in x_default_list:
                    if self.variableIsFound(ncname, x_n):
                        x_found.append(x_n)
                # X case (1)
                if len(x_found) == 0:  # nothing found, ask user to type name
                    print 'X coords: Nothing found. Choose manually.'
                    if long_var_list_not_shot: raw_input('X coords: Press Enter to list all available variables in current file')
                    ndim = -1
                    while ndim not in [1, 2]:
                        x_vname = self.promtVariableNameInput(ncname, promtInputMessage='X coords: Type the name of the X-coords-variable to pick (should be 1D or 2D):', printAvailableVarNames=long_var_list_not_shot)
                        ndim = self.read_netcdfVarMetadata(ncname, x_vname, raise_error=False, promtInput=False)['nNonOneDims']
                        if long_var_list_not_shot: long_var_list_not_shot ^= long_var_list_not_shot
            
                # X case (2)
                elif len(x_found) == 1:
                    x_vname = x_found[0]
            
                # X case (3)
                else:  # len(x_found) > 1
                    print 'X coords: More than 1 X-coords variable exist; variables found:'
                    for x_v_n in x_found:
                        print '\t', x_v_n
                    x_vname = self.promtVariableNameInput(ncname, promtInputMessage='X coords: Type the name of the X-coords-variable to pick:', printAvailableVarNames=False)
            
            else:
                x_vname = dim_bath[-1]
                
        self.__x_coords['meta'] = self.read_netcdfVarMetadata(ncname, x_vname, raise_error=True, promtInput=False)
        if log: print 'X coords: Picking Y-coords variable: <{0}> , dimensions {1}, shape {2}'.format(
                        self.__x_coords['meta']['name'], self.__x_coords['meta']['dims'], self.__x_coords['meta']['shape'])
        # ----------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------



        # Y coords...
        # ----------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------
        if y_vname and not self.variableIsFound(ncname, y_vname):
                raw_input('Y coords: User-defined Y-coord-variable name <{0}> not found in file <{1}>. Press Enter to continue auto-search'.format(y_vname, ncname))

        if not y_vname:
            dim_bath = self.__bathymetry['meta']['dims']
            if log: print 'Y coords: Trying to find variable based on bathymetry dimensions <{0}>. Searching for variable <{1}>'.format(dim_bath, dim_bath[-2] )
            
            if not self.variableIsFound(ncname, dim_bath[-2]):
                if log: print 'Y coords: variable name <{0}> not found in file <{1}>.'.format(dim_bath[-2], ncname)
                y_default_list = self.get_constants()['y_cartesian_vnames']+self.get_constants()['y_geographi_vnames']
                if log: print 'Y coords: Searching for variables with name from default list {0}'.format(y_default_list)
                
                y_found = list()
                for y_n in y_default_list:
                    if self.variableIsFound(ncname, y_n):
                        y_found.append(y_n)
                # Y case (1)
                if len(y_found) == 0:  # nothing found, ask user to type name
                    print 'Y coords: Nothing found. Choose manually.'
                    if long_var_list_not_shot: raw_input('Y coords: Press Enter to list all available variables in current file')
                    ndim = -1
                    while ndim not in [1, 2]:
                        y_vname = self.promtVariableNameInput(ncname, promtInputMessage='Y coords: Type the name of the Y-coords-variable to pick (should be 1D or 2D):', printAvailableVarNames=long_var_list_not_shot)
                        ndim = self.read_netcdfVarMetadata(ncname, y_vname, raise_error=False, promtInput=False)['nNonOneDims']
                        if long_var_list_not_shot: long_var_list_not_shot ^= long_var_list_not_shot
            
                # Y case (2)
                elif len(y_found) == 1:
                    y_vname = y_found[0]
            
                # Y case (3)
                else:  # len(y_found) > 1
                    print 'Y coords: More than 1 Y-coords variable exist; variables found:'
                    for y_v_n in y_found:
                        print '\t', y_v_n
                    y_vname = self.promtVariableNameInput(ncname, promtInputMessage='Y coords: Type the name of the Y-coords-variable to pick:', printAvailableVarNames=False)
            
            else:
                y_vname = dim_bath[-2]
                
        self.__y_coords['meta'] = self.read_netcdfVarMetadata(ncname, y_vname, raise_error=True, promtInput=False)
        if log: print 'Y coords: Picking Y-coords variable: <{0}> , dimensions {1}, shape {2}'.format(
                        self.__y_coords['meta']['name'], self.__y_coords['meta']['dims'], self.__y_coords['meta']['shape'])
        # ----------------------------------------------------------------------------------

        # -------------------------------------------------------
        # now we know var-names for X and Y
        # -------------------------------------------------------
        # now determine if it is rectangular grid or curvilinear....
        #
        # if x-coord and y-coord are 1d arrays
        # then the grid is rectangular but the cells may be of
        # any rectangular shape
        #
        # if x-y- coords are 2d arrays,
        # then the grid is curvilinear
        # NOTE: here someone can paste rectangular grid with two 2d arrays... but it will still work


        if   len(self.__x_coords['meta']['shape']) == 1 and len(self.__y_coords['meta']['shape']) == 1:
            self.__grid_info['type'] = 'rectangular'
        
        elif len(self.__x_coords['meta']['shape']) == 2 and len(self.__y_coords['meta']['shape']) == 2:
            self.__grid_info['type'] = 'curvilinear'
        else:
            print 'GridType: Using variable for X-coords <{0}>, of shape <{1}>'.format(self.__x_coords['meta']['name'], self.__x_coords['meta']['shape'])
            print 'GridType: Using variable for Y-coords <{0}>, of shape <{1}>'.format(self.__y_coords['meta']['name'], self.__y_coords['meta']['shape'])
            raise ValueError('Grid type not understood. X and Y should be either two 1D arrays or two 2D arrays'+'\n'+gridhelp())


        # now determine if coords are at T (cell center) or X (cell nodes) points...GridType:
        print 'Data location: X-coords   variable: <{0}> , dimensions {1}, shape {2}'.format(self.__x_coords['meta']['name'], self.__x_coords['meta']['dims'], self.__x_coords['meta']['shape'])
        print 'Data location: Y-coords   variable: <{0}> , dimensions {1}, shape {2}'.format(self.__y_coords['meta']['name'], self.__y_coords['meta']['dims'], self.__y_coords['meta']['shape'])
        print 'Data location: Bathymetry variable: <{0}> , dimensions {1}, shape {2}'.format(self.__bathymetry['meta']['name'], self.__bathymetry['meta']['dims'], self.__bathymetry['meta']['shape'])
        xy_location = None
        while xy_location not in ['x', 'X', 't', 'T']:
            xy_location = raw_input('Data location: select whether origin XY-COORDS data is located at nodes(X_points) or at cell centers(T_points) [X/T]:')
        
        if xy_location in ['x', 'X']:
            self.__x_coords['meta']['points_location'] = 'X_points'
            self.__y_coords['meta']['points_location'] = 'X_points'
        else:
            self.__x_coords['meta']['points_location'] = 'T_points'
            self.__y_coords['meta']['points_location'] = 'T_points'
        
        bt_location = None
        while bt_location not in ['x', 'X', 't', 'T']:
            bt_location = raw_input('Data location: select whether origin BATHYMETRY data is located at nodes(X_points) or at cell centers(T_points) [X/T]:')
        if bt_location in ['x', 'X']:
            self.__bathymetry['meta']['points_location'] = 'X_points'
        else:
            self.__bathymetry['meta']['points_location'] = 'T_points'

        
        # determine mode.... (Geographic or Cartesian)
        if   self.__x_coords['meta']['name'] in self.get_constants()['x_cartesian_vnames'] and self.__y_coords['meta']['name'] in self.get_constants()['y_cartesian_vnames']:
            self.__x_coords['meta']['coordinate_type'] = 'cartesian'
            self.__y_coords['meta']['coordinate_type'] = 'cartesian'
        elif self.__x_coords['meta']['name'] in self.get_constants()['x_geographi_vnames'] and self.__y_coords['meta']['name'] in self.get_constants()['y_geographi_vnames']:
            self.__x_coords['meta']['coordinate_type'] = 'geographic'
            self.__y_coords['meta']['coordinate_type'] = 'geographic'
        else:
            # now show user found vars
            print 'Coords type: X-coords   variable: <{0}> , dimensions {1}, shape {2}'.format(self.__x_coords['meta']['name'], self.__x_coords['meta']['dims'], self.__x_coords['meta']['shape'])
            print 'Coords type: Y-coords   variable: <{0}> , dimensions {1}, shape {2}'.format(self.__y_coords['meta']['name'], self.__y_coords['meta']['dims'], self.__y_coords['meta']['shape'])
            coord_mode = None
            while coord_mode not in ['c', 'g']:
                coord_mode = raw_input('Coords type: coord-type not understood. Choose cartesian or geographic. Type [c/g]:' )
            if coord_mode == 'c':
                self.__x_coords['meta']['coordinate_type'] = 'cartesian'
                self.__y_coords['meta']['coordinate_type'] = 'cartesian'
            if coord_mode == 'g':
                self.__x_coords['meta']['coordinate_type'] = 'geographic'
                self.__y_coords['meta']['coordinate_type'] = 'geographic'




    def _generate_mask(self, maskvalue=None, transpose=False, log=False):
        ''' Function reads an 2D array (in MOSSCO - 'bathymetry') and generates a boolean mask
        array (True means that value is masked), based on 'fillvalue'

        We have here 2 possibilities depending on data - location:
            1) bathymetry at X_points
            2) bathymetry at T_points

        
                    x-----x-----x-----x-----x
                    |     |     |     |     |
                    |  T  |  T  |  T  |  T  |
                    |     |     |     |     |
                    x-----x-----x-----x-----x
                    |     |     |     |     |
                    |  T  |  T  |  T  |  T  |
                    |     |     |     |     |
                    x-----x-----x-----x-----x
                    |     |     |     |     |
                    |  T  |  T  |  T  |  T  |
                    |     |     |     |     |
                    x-----x-----x-----x-----x
        input:
            filename        - string to filename (relative path)
            varname         - string containing name of the variable
            fillvalue       - float indicating the values to be masked (mask=True). By default (fillvalue=None)
                            gets _FillValue automatically
            transpose       - a boolean flag to return transposed array or not
        '''
        _n = 'generate_mask(): '
        meta = self.get_metadata()
        if log: print _n, 'Working with variable <{0}> of shape {2} from file <{1}>'.format(meta['b']['name'], meta['b']['fname'], meta['b']['shape'])
        
        array = self.get_data(squeeze=True)['b']
        
        if len(array.shape) == 2:
            pass
        else:
            raise ValueError(_n+" Array should be two dimensional. Received {0}-dimensional".format(len(array.shape)))

        location = meta['b']['points_location']

        mask_in = None
        # if it is already masked array, simply take the mask
        if isinstance(array, np.ma.MaskedArray):
            mask_in = np.ma.getmaskarray(array)
            if not ui.promtYesNo(_n+'Array is already of type <numpy.ma.MaskedArray>\n' + _n +
                    'Default mask of shape {0} found. Use it? If "no" i will try to build mask based on _FillValue.'.format(mask_in.shape)):
                mask_in = None

        # ... or create own mask based on fv
        if mask_in is None:
            if maskvalue is None:  # nothing passed by user
                for attr in self.get_constants()['fv_attr_namelist']:
                    if attr in meta['b']['attrs'].keys():
                        maskvalue = float(meta['b']['attrs'][attr])
                        if log: print _n, 'Found attribute <{0}={1}>. I will use this value to generate mask'.format(attr, maskvalue)
                        break
                # if nothing found in <for> , here maskvalue=None
                if not maskvalue:
                    raw_input(_n+'Current virable does not have any of the following attributes {0}. I will continue without mask. Press ENTER to continue'.format(
                                self.get_constants()['fv_attr_namelist']))
                    return None
            masked_arr = np.ma.masked_values(array, float(maskvalue))
            mask_in = np.ma.getmaskarray(masked_arr)
            del masked_arr

        # so... at this point we should have <mask_in> - the input mask
        
        if location == 'X_points':
            # create new face map ( 2D array of T_points), it will have one less index in each dimension

            mask_out = np.zeros(tuple([array.shape[0]-1, array.shape[1]-1]), dtype=bool)  # by default all cells are valid (mask=False)
            # loop over all face-indexes. If any of the surrounding nodes is invalid (mask=True)
            # then the whole face is considered to be invalid as well
            for j in xrange(mask_out.shape[0]):
                for i in xrange(mask_out.shape[1]):
                    tl = mask_in[j  , i  ]  #topleft node
                    tr = mask_in[j  , i+1]  #topright node
                    br = mask_in[j+1, i+1]  #bottomright node
                    bl = mask_in[j+1, i  ]  #bottomleft node
                    if any(node is True for node in [tl, tr, br, bl]):
                        mask_out[j, i] = True

        elif location == 'T_points' :  # location == 'T_points'
            # from numpy docs: "We must keep in mind that a True entry in the mask indicates an invalid data"
                mask_out = mask_in
        else:
            raise ValueError('Invalid data location. Should be <T_points> or <X_points>, received {0}'.format(location))

        # at this point mask is created....
        if any(mask_out.ravel()) > 0:  #True is equal to 1. So if array has at least one masked element, sum will be more than 0
            if transpose: mask_out = mask_out.T
            if log: print _n, 'Created boolean mask of shape {0} (created from array of shape{1})'.format(mask_out.shape, array.shape)
        else:
            if log: print _n, 'Array has no invalid entries. Mask is not needed. Returning <None>'
            mask_out = None  # overwriting nomask value
        
        return mask_out



























        

def create_magnitude_variable_from_x_y_component(VARS, varname, varval, mask=None, log=False):
    '''
    in:
    ----
        VARS     - dictionary generated by function process_cdl.read_file_with_only_variables()
        varname  - string, item-key of a dictionary VARS, should correspond to varval
        varval   - item-value of a dictionary VARS, should correspond to varname
                    varval[0] >>> datatype
                    varval[1] >>> [dim1,dim2,...]
                    varval[2] >>> dict(attributes)
                    Note: all values are stored as strings
    out:
    ----
        magnitude - False, or Ndimensional numpy array (maybe masked, depending on the inputs).
                    It will usually be 1D array (<...> in selections such as x[t, ..., f] is made
                    to increase flexibility of accepted variables. Nevertheless mostly it will be x[t, f])

    '''
    magnitude = None
    # -----------------------------------------------------------------------------------------------
    # Create auto variables
    # -----------------------------------------------------------------------------------------------
    if '_auto_creation' in varval[2].keys():
        if log: print 'Autocreation-variable'
        _fnx = VARS[ varval[2]['_auto_creation'].split(',')[0].strip() ] [2]['_mossco_filename']
        _fny = VARS[ varval[2]['_auto_creation'].split(',')[1].strip() ] [2]['_mossco_filename']
        _vnx = VARS[ varval[2]['_auto_creation'].split(',')[0].strip() ] [2]['_mossco_varname']
        _vny = VARS[ varval[2]['_auto_creation'].split(',')[1].strip() ] [2]['_mossco_varname']

        # -----------------------------------------------------------------------------------------------
        # now check which function (read_mossco_nc_3d, read_mossco_nc_4d) to use to get data....
        # -----------------------------------------------------------------------------------------------
        if varname.endswith('_2d'):
            x, _ = process_mossco_netcdf.read_mossco_nc_3d(_fnx, _vnx, mask=mask)
            y, _ = process_mossco_netcdf.read_mossco_nc_3d(_fny, _vny, mask=mask)
            magnitude = np.zeros(x.shape)
            for t in xrange(x.shape[0]):
                for f in xrange(x.shape[-1]):
                    _x = x[t, ..., f]
                    _y = y[t, ..., f]
                    _magnitude = (_x**2 + _y**2)**(1./2.)
                    magnitude[t, ..., f] = _magnitude

        elif varname.endswith('_3d'):
            x, _ = process_mossco_netcdf.read_mossco_nc_4d(_fnx, _vnx, mask=mask)
            y, _ = process_mossco_netcdf.read_mossco_nc_4d(_fny, _vny, mask=mask)
            magnitude = np.zeros(x.shape)
            for t in xrange(x.shape[0]):
                for z in xrange(x.shape[1]):
                    for f in xrange(x.shape[-1]):
                        _x = x[t, z, ..., f]
                        _y = y[t, z, ..., f]
                        _magnitude = (_x**2 + _y**2)**(1./2.)
                        magnitude[t, ..., f] = _magnitude

    return magnitude




def create_layer_elevation_from_sigma_coords(eta, sigma, depth, flatten=False, mask=None, log=False):
    '''
    create elevations <z> with respect to MSL of passed sigma-coordinates.

    Calculations are performed in accordance with CF-conventions (CF 1.6)
    for variable "Ocean Sigma Coordinate" as...
        z(n,k,j,i) = eta(n,j,i) + sigma(k)*(depth(j,i)+eta(n,j,i))

        where:
            z, eta, sigma, depth - numpy arrays
            n   - integer, timesteps
            k   - integer, layers
            j,i - integer, y,x indices

    flatten - if True, x,y dimensions of the array will be compressed into one single
    mask    - boolean 2d mask to ignore elements during flattening. see <process_mossco_netcdf.make_mask_array_from_mossco_bathymetry()>
    '''
    _name = 'create_layer_elevation_from_sigma_coords():'
    if log:
        print _name, 'Shapes of the inputs...'
        print _name, 'eta: <{0}>, sigma: <{1}>, depth: <{2}>'.format(eta.shape, sigma.shape, depth.shape)
    
    elev         = np.zeros((eta.shape[0], len(sigma), eta.shape[1], eta.shape[2]))
    elev_borders = np.zeros((2, eta.shape[0], len(sigma), eta.shape[1], eta.shape[2]))

    # create sigma coordinates of the borders
    sigma_borders = np.zeros(len(sigma)+1)
    sigma_borders[0] = -1.  # coordinate of the very bottom
    for z in xrange(len(sigma)):
        sigma_borders[z+1] = sigma_borders[z] - (sigma_borders[z] - sigma[z])*2
    

    if abs(sigma_borders[-1]) > 0.005:
        print _name, 'sigma layer centers', sigma
        print _name, 'sigma layer borders', sigma_borders
        raise ValueError('Sigma values for layer-borders calculated not correctly')
    else:
        sigma_borders[-1] = 0.  # corrdinate of the very top

    if log: print _name, 'sigma layer centers', sigma
    if log: print _name, 'sigma layer borders', sigma_borders
    


    for t in xrange(elev.shape[0]):
        for z in xrange(elev.shape[1]):
            elev[t, z, ...] = eta[t, ...] + sigma[z]*(depth + eta[t, ...])
            for border in xrange(2):
                elev_borders[border, t, z, ...] = eta[t, ...] + sigma_borders[z+border]*(depth + eta[t, ...])
            #if log: print 't=', t, 'z=', z, 'elev:', elev[t, z, 12, 12]
            #if log: print 't=', t, 'z=', z, 'elev_bnb [0]:', elev_borders[0, t, z, 12, 12]
            #if log: print 't=', t, 'z=', z, 'elev_bnb [1]:', elev_borders[1, t, z, 12, 12]
    



    if log: print _name, 'Elevation array created of shape <{0}>'.format(elev.shape)
    if log: print _name, 'Elevation border array created of shape <{0}>'.format(elev_borders.shape)
    return elev, elev_borders



def create_sigma_coords_of_layer_center(sigma_border):
    '''
        creates arrays of sigma-coordinates for layer centers, when a
        corresponding array is given for the layer borders
    '''

    sigma_center = np.zeros(len(sigma_border)-1)

    for z in xrange(sigma_center.__len__()):
        sigma_center[z] = .5*(sigma_border[z] + sigma_border[z+1])
    return sigma_center






def flatten_xy_data(data, mask=None):
    data = np.squeeze(data)
    if mask is None:
        
        dims = list(data.shape[:-2])
        dims.append(data.shape[-1]*data.shape[-2])

        a = np.zeros(dims)
        if len(dims) == 3:
            for t in xrange(data.shape[0]):
                for z in xrange(data.shape[1]):
                    a[t, z, :] = data[t, z, ...].flatten(order='F')
        elif len(dims) == 2:
            for t in xrange(data.shape[0]):
                a[t, :] = data[t, ...].flatten(order='F')
        elif len(dims) == 1:
            a[:] = data[...].flatten(order='F')

        elif len(dims) == 4 and dims[0] == 2:  # if we have boundary var (i.e. Mesh2_face_bnd(two, t, z, face))
            del a
            dims.append(2)
            a = np.zeros(dims[1::])  # make the dimension two appear at the end... (two, t, z, face) => (t, z, face, two)

            for t in xrange(data.shape[1]):
                for z in xrange(data.shape[2]):
                    for bnd in xrange(data.shape[0]):
                        a[t, z, :, bnd] = data[bnd, t, z, ...].flatten(order='F')
        else:
            raise ValueError('Number of array dimensions <{0}> is not supported.'.format(len(dims)))
    else:
        n_valid_2d = np.sum(np.invert(mask))  #number of valid elements in 2d part. invert - because True is an invalid element
        
        dims = list(data.shape[:-2])
        dims.append(n_valid_2d)

        a = np.zeros(dims)

        if len(dims) == 3:
            for t in xrange(data.shape[0]):
                for z in xrange(data.shape[1]):
                    var_masked = np.ma.array(data[t, z, ...], mask=mask)
                    var_masked = var_masked.flatten(order='F').compressed()
                    a[t, z, :] = var_masked
        elif len(dims) == 2:
            for t in xrange(data.shape[0]):
                var_masked = np.ma.array(data[t, ...], mask=mask)
                var_masked = var_masked.flatten(order='F').compressed()
                a[t, :] = var_masked
        elif len(dims) == 1:
            var_masked = np.ma.array(data[...], mask=mask)
            var_masked = var_masked.flatten(order='F').compressed()
            a[:] = var_masked
        elif len(dims) == 4 and dims[0] == 2:  # if we have boundary var (i.e. Mesh2_face_bnd(two, t, z, face))
            del a
            dims.append(2)
            a = np.zeros(dims[1::])  # make the dimension two appear at the end... (two, t, z, face) => (t, z, face, two)
            for t in xrange(data.shape[1]):
                for z in xrange(data.shape[2]):
                    for bnd in xrange(data.shape[0]):
                        var_masked = np.ma.array(data[bnd, t, z, ...], mask=mask)
                        var_masked = var_masked.flatten(order='F').compressed()
                        a[t, z, :, bnd] = var_masked
        else:
            raise ValueError('Number of array dimensions <{0}> is not supported.'.format(len(dims)))

    return a
