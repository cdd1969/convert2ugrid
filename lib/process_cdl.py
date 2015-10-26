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
from netCDF4 import Dataset
import numpy as np
import time
import sys
import os
import re
import process_mossco_netcdf
import copy


class cdlVariable(object):
    '''
        class describes variables from cdl file format
    '''
    def __init__(self):
        self.__name  = None  # inits var-name
        self.__dtype = None  # inits var d-type
        self.__dims  = None  # inits var dims
        self.__fv    = None  # inits var fill-value
        self.__attrs = None  # inits var attributes

        self._init_constants()

    def _init_constants(self):
        self.__dtypes_list = ["char", "byte", "short", "int", "long", "float", "real", "double",
                    "ubyte", "ushort", "uint", "int64", "uint64", "string"]

    def set_name(self, name):
        # name - string
        if self._check_dtype(name, str):
            self.__name = name

    def set_dtype(self, dtype):
        # dtype - string
        if self._check_dtype(dtype, str):
            if dtype in self.__dtypes_list:
                self.__dtype = dtype
            else:
                raise TypeError('Invalid datatype. Received <{0}> not found in default list {1}'.format(dtype, self.__dtypes_list))

    def _check_dtype(self, _object, dtype):
        # compare type(<_object>) to <dtype>
        if isinstance(_object, dtype):
            return True
        else:
            raise TypeError('<{0}> should be of type <{2}>. Is {1}'.format(_object, type(_object), str(dtype)))




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
            variables [dict] - variables that have been found, (key - name of the variable, value = [datatype, [dim1,dim2,...], dict(attributes)])
                                Note: all values are stored as strings
    '''
    _n = "read_file_with_only_variables():"
    VARIABLES = dict()
    # this is list with datatypes both supported by netcdf3(first row) and netcdf4(second row)
    dtype_list = ["char", "byte", "short", "int", "long", "float", "real", "double", ]
    #                "ubyte", "ushort", "uint", "int64", "uint64", "string"]

    # init variable to be read
    buffer_var = None

    for line in content_in_lines:
        #print line
        line_contains_var_definition = False
        
        for dt in dtype_list:
            if line.startswith(dt):  # if this is a line with definition a variable
                line_contains_var_definition = True
                
                if buffer_var is not None:  #save previously gathered information, if we have
                    VARIABLES[buffer_var['name']] = [buffer_var['dtype'], buffer_var['dims'], buffer_var['attrs']]

                # clear buffer
                buffer_var          = dict()  # initialize
                buffer_var['attrs'] = dict()  # initialize
                
                # save var-name
                buffer_var['name']  = line.split()[1].strip()
                
                # save data-type
                buffer_var['dtype'] = line.split()[0].strip()

                # save dims
                if '(' in line:  # if variable has dimensions (because they are declared in parenthesis "()")
                    buffer_var['name']  = line.split('(')[0].split()[1].strip()  # save var-name
                    vardims = re.match('.*?\((.*?)\).*', line).group(1).strip()
                    if len(vardims) == 0:  #checking length of the string! number of symbols! if we have ()
                        buffer_var['dims'] = list()  # empty dims
                    elif ',' not in vardims:  # we have something like (dim1)
                        buffer_var['dims'] = [vardims]
                    else:  # we have something like (dim1, dim2, dim3)
                        vardims = vardims.split(',')
                        buffer_var['dims'] = [vardim.strip() for vardim in vardims]  # real dims
                else:
                    buffer_var['dims'] = list()  # empty dims
        
        if not line_contains_var_definition:
            if line.startswith(buffer_var['name']):  # if this is a line with attribute
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

                buffer_var['attrs'][attr_name] = attr_val_proper_type

            elif line not in ['\n', '']:
                raise IOError('{1}\n{2}\nline below not understood (should either start with <{3}> or <datatype_name(i.e float, double)>):\n{0}\n{2}\n'.format(line, _n, '-'*50, buffer_var['name']))

        #add last variable
        if line is content_in_lines[-1]:  #if we have finished proceeding the last line
            if buffer_var is not None:  # if we have variable in buffer
                VARIABLES[buffer_var['name']] = [buffer_var['dtype'], buffer_var['dims'], buffer_var['attrs']]
    
    return VARIABLES


def print_nicely_variables(vars_dict):
    for vn, vv in vars_dict.iteritems():
        print '\n'
        print re.sub(r'[\[\]\u\']', '', '{0} {1} ({2}) ;'.format(vv[0], vn, vv[1]))
        if vv[2]:
            for an, av in vv[2].iteritems():
                if isinstance(av, unicode) or isinstance(av, str):
                    print '\t{0}: {1} = "{2}" ;'.format(vn, an, av)
                elif isinstance(av, list):
                    if isinstance(av[0], float):
                        string = '\t{0}: {1} = '.format(vn, an)
                        for item in av:
                            #
                            # HERE IMPLEMENT   Number. This is the same as 'g', except that it uses the
                            # current locale setting to insert the appropriate
                            # number separator characters.
                            #
                            # https://www.python.org/dev/peps/pep-3101/
                            item_str = '{:g}'.format(item)
                            if 'e' in item_str:
                                item_str = re.sub('e', '.e', item_str)
                            elif '.' not in item_str:
                                item_str += '.'
                            string += '{0}, '.format(item_str)
                        string = string[0:-2]+' ;'

                    else:
                        string = '\t{0}: {1} = {2} ;'.format(vn, an, av)
                    print re.sub(r'[\[\]]', '', string)
                else:
                    print '\t{0}: {1} = {2} ;'.format(vn, an, av)



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
            f.write('\n')
            f.write('// '+'-'*100+'\n')
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


def create_cdl_file(cdl_baw, txt_mossco_baw, output_fname, log=False):
    """
        script processes information and creates a CDL sample file with description of variables.

        input:
            cdl_baw             - string, containing name of CDL file, where BAW standart variables are described
            txt_mossco_baw      - string, containing name of TXT file, correlation between BAW standart variable
                                  names and MOSSCO variable names is specified
            output_fname        - string containing name of output file
    """
    cdl_baw, txt_mossco_baw = os.path.abspath(cdl_baw), os.path.abspath(txt_mossco_baw)
    MOSSCO_VARS_TO_INCLUDE = read_txt_mossco_baw(txt_mossco_baw, log=log)
    #MOSSCO_VARS_TO_INCLUDE <<< [dict of strings]  (key - name of the variable, value = [mossco_filepath, mossco_varname, davit_dims])
    BAW_VARS_DESCRIPTION = read_file_with_only_variables(read_file_delete_comments(cdl_baw, comments='//'))
    #BAW_VARS_DESCRIPTION <<< [dict of strings]  (key - name of the variable, value = [datatype, [dim1,dim2,...], dict(attributes)])

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


        # auto-velocity magnitude
        # IT IS IMPORTANT TO USE 2d or 3d at the end!!! this will be used to choose proper function for reading data.. 2d or 3d
        auto_vars = dict()
        auto_vars['Mesh2_face_WindVelocityAt10m_m_2d']       = ['Mesh2_face_WindVelocityAt10m_x_2d', 'Mesh2_face_WindVelocityAt10m_y_2d']
        auto_vars['Mesh2_face_WaterVelocityAtSoil_m_2d']     = ['Mesh2_face_WaterVelocityAtSoil_x_2d', 'Mesh2_face_WaterVelocityAtSoil_y_2d']
        auto_vars['Mesh2_face_Stroemungsgeschwindigkeit_m_2d'] = ['Mesh2_face_Stroemungsgeschwindigkeit_x_2d', 'Mesh2_face_Stroemungsgeschwindigkeit_y_2d']
        auto_vars['Mesh2_face_Stroemungsgeschwindigkeit_m_3d'] = ['Mesh2_face_Stroemungsgeschwindigkeit_x_3d', 'Mesh2_face_Stroemungsgeschwindigkeit_y_3d']
        
        _AUTO_VARS = copy.deepcopy(auto_vars)


        #iterate over variables to include
        for vn, vv in MOSSCO_VARS_TO_INCLUDE.iteritems():

            # removing parents
            for auto_vn, auto_vv in auto_vars.iteritems():
                if vn in auto_vv:
                    auto_vars[auto_vn].remove(vn)

            # look up corresponding varname
            if vn in BAW_VARS_DESCRIPTION.keys():
                
                # write variable definition
                if vv[1]:  # if variable with dimensions...
                    f.write(re.sub(r'[\[\]\']', '', '{0} {1} ({2}) ;\n'.format(BAW_VARS_DESCRIPTION[vn][0], vn, BAW_VARS_DESCRIPTION[vn][1])) )
                else:  # if variable without dimensions...
                    f.write('{0} {1} ;\n'.format(BAW_VARS_DESCRIPTION[vn][0], vn) )
                
                # write out all standart attributes specified in baw cdl
                if BAW_VARS_DESCRIPTION[vn][2]:
                    for an, av in BAW_VARS_DESCRIPTION[vn][2].iteritems():
                        # specially treat _FillValue...
                        if an == '_FillValue':
                            # if we have autovariable, it is not existing... so we have to use <try>
                            try:
                                nc = Dataset(vv[0], mode='r')
                                av = nc.variables[vv[1]]._FillValue
                                nc.close()
                                f.write('\t{0}: {1} = {2} ;\n'.format(vn, an, av) )
                            except:
                                if isinstance(av, list):
                                    if isinstance(av[0], float):
                                        string = '\t{0}: {1} = '.format(vn, an)
                                        for item in av:
                                            item_str = '{:g}'.format(item)
                                            if 'e' in item_str:
                                                item_str = re.sub('e', '.e', item_str)
                                            elif '.' not in item_str:
                                                item_str += '.'
                                            string += '{0}, '.format(item_str)
                                        string = string[0:-2]+' ;'

                                    else:
                                        string = '\t{0}: {1} = {2} ;'.format(vn, an, av)
                                    f.write( re.sub(r'[\[\]]', '', string)+'\n' )
                                else:
                                    f.write(  '\t{0}: {1} = {2} ;\n'.format(vn, an, av) )


                        # if the current attribute is not _FillValue    
                        else:
                            if isinstance(av, unicode) or isinstance(av, str):
                                f.write('\t{0}: {1} = "{2}" ;\n'.format(vn, an, av) )
                            elif isinstance(av, list):
                                if isinstance(av[0], float):
                                    string = '\t{0}: {1} = '.format(vn, an)
                                    for item in av:
                                        item_str = '{:g}'.format(item)
                                        if 'e' in item_str:
                                            item_str = re.sub('e', '.e', item_str)
                                        elif '.' not in item_str:
                                            item_str += '.'
                                        string += '{0}, '.format(item_str)
                                    string = string[0:-2]+' ;'

                                else:
                                    string = '\t{0}: {1} = {2} ;'.format(vn, an, av)
                                f.write( re.sub(r'[\[\]]', '', string)+'\n' )
                            else:
                                f.write(  '\t{0}: {1} = {2} ;\n'.format(vn, an, av) )

            else:
                # (nMesh2_data_time)
                # (nMesh2_data_time, nMesh2_face)
                # (nMesh2_data_time, nMesh2_layer_2d, nMesh2_face)
                # (nMesh2_data_time, nMesh2_layer_3d, nMesh2_face)
                # (nMesh2_data_time, nMesh2_suspension_classes, nMesh2_layer_2d, nMesh2_face)
                # (nMesh2_data_time, nMesh2_suspension_classes, nMesh2_layer_3d, nMesh2_face)
                raise EOFError('a variable "{0}" not found in:\n{1}\nbut is present in:\n{2}\n'.format(vn, cdl_baw, txt_mossco_baw))
                """
                if vv[2]:  # if variable with dimensions...
                    f.write(re.sub(r'[\[\]\u\']', '', '{0} {1} ({2}) ;\n'.format('float', vn, vv[2])) )
                else:
                    f.write(re.sub(r'[\[\]\u\']', '', '{0} {1} ;\n'.format('float', vn)) )
                """
            

            #now write additional attributes
            f.write(  '\t{0}: _mossco_filename = "{1}" ;\n'.format(vn, vv[0]) )
            f.write(  '\t{0}: _mossco_varname = "{1}" ;\n'.format(vn, vv[1]) )
            f.write(  '\n\n')



        # CREATING AUTOMATICALLY VARIABLES
        #
        #
        #
        #
        #if skip = [] ==> create automatically variable
        for vn, skip in auto_vars.iteritems():
            if not skip:
                if vn in BAW_VARS_DESCRIPTION.keys():
                    
                    # write variable definition
                    if vv[1]:  # if variable with dimensions...
                        f.write(re.sub(r'[\[\]\']', '', '{0} {1} ({2}) ;\n'.format(BAW_VARS_DESCRIPTION[vn][0], vn, BAW_VARS_DESCRIPTION[vn][1])) )
                    else:  # if variable without dimensions...
                        f.write('{0} {1} ;\n'.format(BAW_VARS_DESCRIPTION[vn][0], vn) )

                    if BAW_VARS_DESCRIPTION[vn][2]:
                        for an, av in BAW_VARS_DESCRIPTION[vn][2].iteritems():
                            if isinstance(av, unicode) or isinstance(av, str):
                                f.write('\t{0}: {1} = "{2}" ;\n'.format(vn, an, av) )
                            elif isinstance(av, list):
                                if isinstance(av[0], float):
                                    string = '\t{0}: {1} = '.format(vn, an)
                                    for item in av:
                                        item_str = '{:g}'.format(item)
                                        if 'e' in item_str:
                                            item_str = re.sub('e', '.e', item_str)
                                        elif '.' not in item_str:
                                            item_str += '.'
                                        string += '{0}, '.format(item_str)
                                    string = string[0:-2]+' ;'
                                else:
                                    string = '\t{0}: {1} = {2} ;'.format(vn, an, av)
                                f.write( re.sub(r'[\[\]]', '', string)+'\n' )
                            else:
                                f.write(  '\t{0}: {1} = {2} ;\n'.format(vn, an, av) )
                    
                    f.write('\t{0}: _auto_creation = "{1}, {2}" ;\n'.format(vn, _AUTO_VARS[vn][0], _AUTO_VARS[vn][1]) )
                    f.write('\n\n')
                else:
                    raise EOFError('a variable "{0}" not found in:\n{1}\nbut is present in:\n{2}\n'.format(vn, cdl_baw, txt_mossco_baw))
        f.close()




if __name__ == '__main__':

    NCfname = [os.path.join(os.path.dirname(sys.argv[0]), '../data/NSBS/'+name) for name in ['netcdf_reference_3d.nc'] ]
    
    TXTfname = os.path.join(os.path.dirname(sys.argv[0]), '../user_input/test_mossco_baw.txt')
    CDL_outputfname = os.path.join(os.path.dirname(sys.argv[0]), '../output/test_mossco_to_davit.cdl')
    CDL_inputfname = os.path.join(os.path.dirname(sys.argv[0]), '../user_input/baw_standart.cdl')
    create_txt_mossco_baw(NCfname, TXTfname)
    create_cdl_file(CDL_inputfname, TXTfname, CDL_outputfname)
    