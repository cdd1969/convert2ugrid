#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# This file is part of "convert2ugrid" tool
#
# Author: Nikolai Chernikov, nikolai.chernikov.ru@gmail.com


'''
This module contains functions relevant for processing CDL format file, which is being used
for further creation of netcdf file compatible with DAVIT
'''
from __future__ import division
import numpy as np
import sys
import os
import re

import process_mossco_netcdf


class cdlVariable(object):
    '''
        class describes meta-data of variables from cdl file format.
        It does not support storage of any large-arrays (i.e. data).
        Data stored within the object can be nicely printed in cdl syntax
        by <print> command:
                >>> a = cdlVariable()
                >>> print a

    '''
    def __init__(self, _id=-1):
        self.__id    = int(_id)  # currenty not used, but in future may be useful
        self.__name  = None      # inits var-name
        self.__dtype = None      # inits var d-type
        self.__dims  = None      # inits var dims
        self.__fv    = None      # inits var fill-value
        self.__attrs = dict()    # inits var attributes

        self._init_constants()



    def _init_constants(self):
        # dictionary <self.__dtypes_cdl2pync> connects valid datatype names
        #       from python-netcdf syntax: <http://unidata.github.io/netcdf4-python/#netCDF4.Dataset.createVariable>
        #       with names in CDL  syntax: <https://www.unidata.ucar.edu/software/netcdf/docs/cdl_data_types.html>
        # dict.key = name in CDL. dict.value = name in python-netcdf API
        # NOTE! here CompondType and VLType not implemented!
        self.__dtypes_cdl2pync = dict()
        self.__dtypes_cdl2pync['char']   = 'c'
        self.__dtypes_cdl2pync['byte']   = 'b'
        self.__dtypes_cdl2pync['short']  = 'i2'
        self.__dtypes_cdl2pync['int']    = 'i4'
        self.__dtypes_cdl2pync['long']   = 'i4'      # deprecated, synonymous with int
        self.__dtypes_cdl2pync['float']  = 'f4'
        self.__dtypes_cdl2pync['real']   = 'f4'      # synonymous with real
        self.__dtypes_cdl2pync['double'] = 'f8'
        self.__dtypes_cdl2pync['ubyte']  = 'u1'      # only supported in netCDF4
        self.__dtypes_cdl2pync['ushort'] = 'u2'      # only supported in netCDF4
        self.__dtypes_cdl2pync['uint']   = 'u4'      # only supported in netCDF4
        self.__dtypes_cdl2pync['int64']  = 'int64'   # only supported in netCDF4
        self.__dtypes_cdl2pync['uint64'] = 'uint64'  # only supported in netCDF4
        self.__dtypes_cdl2pync['string'] = str       # only supported in netCDF4. Python built-in type <str>

        # __dtypes_list is list with datatypes in CDL syntax supported by both netcdf3 and netcdf4
        self.__dtypes_list = self.__dtypes_cdl2pync.keys()


    def check_dtype(self, _object, dtype, raise_error=True):
        # compare type(<_object>) to <dtype>
        if isinstance(_object, dtype):
            return True
        else:
            # here we treat unicode and simple string as same objects
            if isinstance(_object, unicode) and dtype is str:
                return True
            if raise_error: raise TypeError('<{0}> should be of type <{2}>. Is {1}'.format(_object, type(_object), str(dtype)))
            else: return False

    def _isInt(self, val):
        # currently we dont care if it is i32 or i64
        if type(val) in [int, long, np.int_, np.intc, np.intp, np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64]:
            return True
        else:
            return False
    
    def _isFloat(self, val):
        # currently we dont care if it is f32 or f64
        if type(val) in [float, np.float16, np.float32, np.float64]:
            return True
        else:
            return False

    def set_name(self, name):
        # name - string
        self.check_dtype(name, str)
        self.__name = name

    def set_dtype(self, dtype):
        # dtype - string
        self.check_dtype(dtype, str)
        if dtype in self.__dtypes_list:
            self.__dtype = dtype
        else:
            raise TypeError('Invalid datatype. Received <{0}> not found in default list {1}'.format(dtype, self.__dtypes_list))

    def set_fillvalue(self, fv):
        # fv - None, integer or float. If not any of them will try to set to float by default
        if hasattr(fv, '__iter__') and len(fv) == 1:
            self.set_fillvalue(fv[0])  #if it is a one-element list, recurse it
        else:
            if fv is None:          self.__fv = fv
            elif self._isInt(fv):   self.__fv = int(fv)
            elif self._isFloat(fv): self.__fv = float(fv)
            else:
                try:    self.__fv = float(fv)
                except: raise TypeError('Invalid fillvalue type. Received <{0}> is of type {1}. Cannot convert to float'.format(fv, type(fv)))
            # now if no errors occured
            self.set_attr('_FillValue', fv)

    def set_dims(self, dims):
        # dims list of strings
        if not hasattr(dims, '__iter__'):
            raise TypeError('Invalid dimension list type. Passed parameter <{0}> is not iterable. Use list() or tuple()'.format(dims))
        lst = list()
        for d in dims:
            self.check_dtype(d, str)
            lst.append(d)
        self.__dims = tuple(lst)



    def set_attr(self, name, val):
        # name - string, val - any type
        self.check_dtype(name, str)
        if name in self.get_attrs().keys():
            del self.__attrs[name]
        if name == '_FillValue' and val != self.get_fillvalue():
            self.set_fillvalue(val)

        self.__attrs[name] = val

    def is_valid_variable(self):
        if any(prop is None for prop in [self.get_name(), self.get_dtype(), self.get_dims()]):
            return False
        else:
            return True


    def get_name(self):
        return self.__name
    
    def get_dims(self):
        return self.__dims

    def get_fillvalue(self):
        return self.__fv
    
    def get_attrs(self):
        return self.__attrs
    
    def get_dtype(self, syntax='cdl'):
        """ Return variable datatype
                either in cdl syntax >>> string
                or in python-netcdf  >>> string or <str> built-in-type
            Input:
                syntax='cdl'  (default)
                syntax='python-netcdf'
        """
        if self.__dtype is None: return None
        else:
            if syntax == 'cdl':
                return self.__dtype
            elif syntax == 'python-netcdf':
                return self.__dtypes_cdl2pync[self.__dtype]
            else:
                raise ValueError('Invalid input for parameter <syntax>. Received <{0}>\n{1}'.format(syntax, self.get_dtype().__doc__))

    def get_string(self):
        """ Here we create nice string to print all info ala CDL
        """
        if not self.is_valid_variable():
            return 'Invalid cdlVariable(): (name, dtype, dims)=({0}, {1}, {2})'.format(self.get_name(), self.get_dtype(), self.get_dims())
        attrs_str = str()
        # join dict if there are entries
        if  self.get_attrs():
            for attr_n in sorted(self.get_attrs().keys()):
                # if string - add apostrofs around str >>> "str"
                if isinstance(self.get_attrs()[attr_n], str):
                    a = '\t{0}:{1} = {2} ;\n'.format(self.get_name(), attr_n, '"'+self.get_attrs()[attr_n]+'"')
                # if iterable - show without brackets [1, 2, 3] >>> 1, 2, 3
                elif hasattr(self.get_attrs()[attr_n], '__iter__'):
                    a = '\t{0}:{1} = {2} ;\n'.format(self.get_name(), attr_n, ', '.join([str(item) for item in self.get_attrs()[attr_n]]))
                # everything else
                else:
                    a = '\t{0}:{1} = {2} ;\n'.format(self.get_name(), attr_n, str(self.get_attrs()[attr_n]))
                attrs_str += a
        # no entries in dict
        else:
            attrs_str = ''


        s = '{0} {1}{2};\n{3}'.format(
            self.get_dtype(),
            self.get_name(),
            '('+', '.join(self.get_dims())+') ' if self.get_dims() else ' ',  # if dims is None or len(dims)==0
            attrs_str
            )
        return s
    

    def get_dtypes_list(self):
        """ These are possible names of variable data-type
        """
        return self.__dtypes_list
    

    def __str__ (self):
        return self.get_string()






















def remove_comments(string):
    pattern = r"(\".*?\"|\'.*?\')|(//[^\n]*$)"
    # first group captures quoted strings (double or single)
    # second group captures comments (//single-line)
    regex = re.compile(pattern, re.DOTALL)
    def _replacer(match):
        # if the 2nd group (capturing comments) is not None,
        # it means we have captured a non-quoted (real) comment string.
        if match.group(2) is not None:
            return ""  # so we will return empty to remove the comment
        else:  # otherwise, we will return the 1st group
            return match.group(1)  # captured quoted-string
    string_fixed = regex.sub(_replacer, string)
    #print '<', string
    #print '>', string_fixed
    #print '-'*5
    return string_fixed


def read_file_delete_comments(path, comments=None):
    '''
    function read a file line by line, excluding comment and blank lines and then combines list of lines

    input:
        path     [str] - full path to file to be read
        comments [str] - a string whith which comment lines start, default is None
    out:
        file_content_str [list with strings]
    '''

    with open(path, 'r') as f:
        file_content = f.readlines()
        f.close()
    if comments in [u'//', '//']:  # remove comment lines that start with "//"
        file_content = [remove_comments(line).strip() for line in file_content if remove_comments(line).strip() not in ['', '\n', '\r']]
        return file_content
    else:
        raise KeyError('Unsupported comment character <{0}>. Should be <//>'.format(comments))


def read_file_with_only_variables(content_in_lines, log=False):
    '''
        function seraches within a given string for variables (see regex) and returns a dictionary with their list

        Simplifications (!!!):
            - content_in_lines should consist of stripped lines.
            - Each line should be terminated at the end(!!!) with semicolon ";". Multiline defenitions are not allowed
    
        input:
            content_in_lines [list] - list with strings representing lines, note that each line should be stripped, to successfully apply .startswith()
        out:
            variables [dict] - variables that have been found, (key - name of the variable, value = object of type <cdlVariable>
    '''
    _n = "read_file_with_only_variables():"
    VARIABLES = dict()
    
    # init variable to be read
    buffer_var = cdlVariable()

    for line in content_in_lines:
        #print line
        line_contains_var_definition = False
        
        for dt in buffer_var.get_dtypes_list():
            if line.startswith(dt):  # if this is a line with definition a variable
                line_contains_var_definition = True
                break
        
        if line_contains_var_definition:
                # save previously gathered information, if we have it
                if buffer_var.get_name() is not None:
                    VARIABLES[buffer_var.get_name()] = buffer_var

                # clear buffer, after we have saved previous info
                del buffer_var
                buffer_var = cdlVariable()
                
                buffer_var.set_name( line.split()[1].strip() )
                buffer_var.set_dtype( line.split()[0].strip() )

                # save dims
                if '(' in line:  # if variable has dimensions (because they are declared in parenthesis "()")
                    buffer_var.set_name( line.split('(')[0].split()[1].strip() )  # save var-name
                    vardims = re.match('.*?\((.*?)\).*', line).group(1).strip()
                    if len(vardims) == 0:  #checking length of the string! number of symbols! if we have ()
                        buffer_var.set_dims( tuple() )  # empty dims
                    elif ',' not in vardims:  # we have something like (dim1)
                        buffer_var.set_dims( [vardims] )
                    else:  # we have something like (dim1, dim2, dim3)
                        vardims = vardims.split(',')
                        buffer_var.set_dims( [vardim.strip() for vardim in vardims] )  # real dims
                else:
                    buffer_var.set_dims( tuple() )
        
        elif not line_contains_var_definition:
            if line.startswith(buffer_var.get_name()):  # if this is a line with attribute
                attr_line = ':'.join(line.split(':')[1::])
                attr_name = attr_line.split('=', 2)[0].strip()

                attr_val_str = '='.join(attr_line.split('=')[1::]) [0:-1].strip()  #removing last symbol, since the line is stripped it is semicolon ';'

                #now convert attr_val_str to a proper datatype
                if '"' in attr_val_str or "'" in attr_val_str:
                    # attribute value is a string
                    attr_val_proper_type = attr_val_str.strip()[1:-1].strip()
                elif ',' in attr_val_str:  # not a string but comma is present ==> list
                    # attribute value is a list
                    attr_val_str = attr_val_str.split(',')
                    attr_val_str = [item.strip() for item in attr_val_str]  # ensure that there are no whitespaces
                else:
                    attr_val_str = [attr_val_str]

                # now we have an attribute either a string of list of values
                
                if isinstance(attr_val_str, list):  # if we have list
                    attr_val_proper_type = []
                    for item in attr_val_str:
                        if ('.' in item) or ('e' in item) or ('f' in item):
                            try:
                                item = np.float32(item)
                            except Exception as err1:
                                if 'f' in item:
                                    try:
                                        item = np.float32(re.sub('f', '', item))
                                    except Exception as err2:
                                        print 'WARNING! : line <{0}> not understood. Could not convert <{1}> to float()'.format(line, item)
                                        print 'WARNING! : therefore I tried deleting "f" symbols. Still cannot convert <{1}> to float()'.format(re.sub('f', '', item))
                                        raise err2
                                else:
                                    print 'WARNING! : line <{0}> not understood. Cannot convert <{1}> to float()'.format(line, item)
                                    raise err1
                        else:
                            try:
                                item = int(item)
                            except Exception as err:
                                print 'WARNING! : line <{0}> not understood. Cannot convert <{1}> to int()'.format(line, item)
                                raise err
                        attr_val_proper_type.append(item)

                buffer_var.set_attr(attr_name, attr_val_proper_type)

            elif line not in ['\n', '']:
                raise IOError('{1}\n{2}\nline below not understood (should either start with <{3}> or <datatype_name(i.e float, double)>):\n{0}\n{2}\n'.format(line, _n, '-'*50, buffer_var['name']))

        #add last variable
        if line is content_in_lines[-1]:  #if we have finished proceeding the last line
            if buffer_var is not None:  # if we have variable in buffer
                VARIABLES[buffer_var.get_name()] = buffer_var
    
    return VARIABLES



def read_baw_mossco_varname_dictionary(fname, log=False):
    f = read_file_delete_comments(fname, comments='//')
    BAW_MOSSCO_VARNAMES = dict()
    for l in f:
        try:
            baw_vn = l.split('>>>')[0].strip()
            baw_vn = re.sub('[\"\']', '', baw_vn)
            
            mossco_vn = l.split('>>>')[1].strip()
            mossco_vn = re.sub('[\"\']', '', mossco_vn)
        
            BAW_MOSSCO_VARNAMES[baw_vn] = [mossco_vn]
        except:
            if log: print 'read_baw_mossco_varname_dictionary(): Skipping line "{0}"'.format(l)
            pass
    return BAW_MOSSCO_VARNAMES


def create_txt_mossco_baw(list_with_ncfnames, output_fname, baw_mossco_varname_dictionary, log=False):
    """
    BAW_MOSSCO_VARNAMES = dict()
    BAW_MOSSCO_VARNAMES['Mesh2_edge_Stroemungsgeschwindigkeit_x_2d'] = ['depth_averaged_x_velocity_in_water']
    BAW_MOSSCO_VARNAMES['Mesh2_edge_Stroemungsgeschwindigkeit_y_2d'] = ['depth_averaged_y_velocity_in_water']
    BAW_MOSSCO_VARNAMES['Mesh2_Temperatur_3d'] = ['temperature_in_water']
    BAW_MOSSCO_VARNAMES['Mesh2_face_depth_2d'] = ['bathymetry']
    BAW_MOSSCO_VARNAMES['Mesh2_Wasserstand_2d'] = ['water_depth_at_soil_surface']
    BAW_MOSSCO_VARNAMES['Mesh2_WaveHeight_2d'] = ['wave_height']
    BAW_MOSSCO_VARNAMES['nMesh2_data_time'] = ['time']
    """

    BAW_MOSSCO_VARNAMES = read_baw_mossco_varname_dictionary(baw_mossco_varname_dictionary, log=log)
    
    if log:
        print '-'*50
        print 'BAW_MOSSCO_DICT'
        for i, v in BAW_MOSSCO_VARNAMES.iteritems():
            print i, '>>>', v
        print 'BAW_MOSSCO_DICT'
        print '-'*50

    with open(output_fname, 'w+') as f:
        f.write('// file describes corelation between mossco-output-netcdf variable names and\n')
        f.write('// baw-format variable names\n')
        f.write('//\n')
        f.write('// format: "filename", "mossco variable name" >>> "corresponding baw format variable name"\n')
        f.write('// format: "filename", "mossco variable name" >>> NOT_INCLUDED\n')
        f.write('// spaces and tabs may be used freely for readability. Comments may follow after "//"\n')
        f.write('//\n')
        f.write('// Displaying Davit-friendly variables found in MOSSCO output netcdf file\n')
        f.write('// Let Davit-friendly variables be those, which have dimensions:\n// 1D: (tdim)\n// 2D: (ydim, xdim)\n// 3D: (tdim, ydim, xdim)\n// 4D: (tdim, zdim, ydim, xdim)\n')
        f.write('//\n')
        f.write('// WARNING! Currently these dimensions are !HARD-CODED! in function get_davit_friendly_variables()\n')

        f.write('// '+'-'*100+'\n')

        for nc in list_with_ncfnames:
            if log: print 'searching for davit-friendly variables in file:', nc
            f.write('\n// '+'-'*100+'\n')
            nc_friendly_vars = process_mossco_netcdf.get_davit_friendly_variables(nc, log=log)
            for typ in ['1D', '2D', '3D', '4D']:
                f.write('//\t'+typ+'\n')
                for v in nc_friendly_vars[typ]:
                    for BAW_VN, MSC_VN in BAW_MOSSCO_VARNAMES.iteritems():
                        if v in MSC_VN:
                            found_baw_vn = '"'+BAW_VN+'"'
                            break
                        else:
                            found_baw_vn = 'NOT_INCLUDED'
                    f.write('   "{0}",\t{1:70s} >>> {2}\n'.format(os.path.abspath(nc), '"'+v+'"', found_baw_vn))
        
        f.close()



def read_txt_mossco_baw(txt_mossco_baw, log=False):
    '''
    in:
        txt_mossco_baw - string, path to dictionary 2 (see documentation)
    out:
        dictionary {'baw_variable_name': ['path_to_mossco_netcdf', 'mossco_variable_name']}
    '''
    VARS = dict()
    f = read_file_delete_comments(txt_mossco_baw, comments='//')
    for i, l in enumerate(f):
        if log: print l
        try:
            baw_vn = l.split('>>>')[1].strip()
            if baw_vn != 'NOT_INCLUDED':
                baw_vn = re.sub('[\"\']', '', baw_vn)
                mossco_vn = re.sub('[\"\']', '', l.split(',')[1].split(">>>")[0].strip())
                mossco_nc = re.sub('[\"\']', '', l.split(',')[0].strip())
                VARS[baw_vn] = [mossco_nc, mossco_vn]
        except Exception, err:
            print 'Error reading line <{0}>: {1}'.format(i, l)
            print(sys.exc_info())
            raise err
    return VARS


def create_dictionary4(cdl_baw, txt_mossco_baw, output_fname, log=False):
    """
        script processes information and creates a CDL sample file with description of variables.

        input:
            cdl_baw             - string, containing name of CDL file, where BAW standart variables are described
            txt_mossco_baw      - string, containing name of TXT file, correlation between BAW standart variable
                                  names and MOSSCO variable names is specified
            output_fname        - string containing name of output file
    """
    cdl_baw, txt_mossco_baw = cdl_baw, txt_mossco_baw
    MOSSCO_VARS_TO_INCLUDE = read_txt_mossco_baw(txt_mossco_baw, log=log)
    #MOSSCO_VARS_TO_INCLUDE <<< [dict of strings]  (key - name of the variable, value = [mossco_filepath, mossco_varname, davit_dims])
    BAW_VARS_DESCRIPTION = read_file_with_only_variables(read_file_delete_comments(cdl_baw, comments='//'))

    with open(output_fname, mode='w+') as f:
        f.write('// This file has been generated by script, which couples BAW-synoptic data standart format with MOSSCO variable names\n')
        f.write('// This file is needed for further conversion of MOSSCO-netcdf file to DAVIT-friendly-netcdf file, and serves as an \n')
        f.write('// input for another processing script\n')
        f.write('//\n')
        f.write('// Below follows the description of synoptic data, presented in a form of NETCDF variables.\n')
        f.write('// Only those variables, listed below will be included in DAVIT-friendly netcdf.\n')
        f.write('// Note, that every variable has two additional non-standart string attributes: "_mossco_filename", "_mossco_varname".\n')
        f.write('// These attributes will not be included in DAVIT-friendly-netcdf, while only indicating the exact location of\n')
        f.write('// data for processing script. For example:\n')
        f.write('// \t\t _mossco_filename = "D:\data\mossco_output.nc"\n')
        f.write('// \t\t _mossco_varname = "depth_averaged_x_velocity_in_water"\n')
        f.write('// The standart BAW attributes, dimensions, and names are well described here:\n')
        f.write('// \thttp://www.baw.de/methoden/index.php5/NetCDF_Synoptische_Daten_im_Dreiecksgitter\n')
        f.write('// '+'-'*100+'\n\n')



        #iterate over variables to include
        for vn, vv in MOSSCO_VARS_TO_INCLUDE.iteritems():
            if vn in BAW_VARS_DESCRIPTION.keys():
                # write variable definition
                var = BAW_VARS_DESCRIPTION[vn]
                var.set_attr('_source_filename', vv[0])
                var.set_attr('_source_varname',  vv[1])
                f.write( BAW_VARS_DESCRIPTION[vn].get_string() )
                
            else:
                # (nMesh2_data_time)
                # (nMesh2_data_time, nMesh2_face)
                # (nMesh2_data_time, nMesh2_layer_2d, nMesh2_face)
                # (nMesh2_data_time, nMesh2_layer_3d, nMesh2_face)
                # (nMesh2_data_time, nMesh2_suspension_classes, nMesh2_layer_2d, nMesh2_face)
                # (nMesh2_data_time, nMesh2_suspension_classes, nMesh2_layer_3d, nMesh2_face)
                raise EOFError('a variable "{0}" not found in:\n{1}\nbut is present in:\n{2}\n'.format(vn, cdl_baw, txt_mossco_baw))

        f.close()




if __name__ == '__main__':

    a = cdlVariable()
    a.set_name("Mesh2_face_test")
    a.set_dtype("float")
    a.set_dims(['time', 'z', 'y_t', 'x_t'])
    a.set_fillvalue(-999)
    a.set_attr('test_attr1', 12.32)
    a.set_attr('test_attr3', 'help i need somebody help not just anybody')
    a.set_attr('test_attr2', 1.3e31)
    a.set_attr('missing_value', -999)
    a.set_attr('test_attr4', True)
    a.set_attr('test_attr5', [1, 2])

    b = cdlVariable()
    b.set_name("Mesh2_face_grid_type")
    b.set_dtype("int")
    b.set_dims([])
    #b.set_attr('long_name', 'zaebalo progat!')


    c = cdlVariable()
    c.set_name("my_var_name")
    print a
    print b
    print c
    

    with open('TES1T', mode='r') as f:
        f.write(a.get_string())
        f.write(b.get_string())
        f.close()
