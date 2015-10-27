#!/usr/bin/python
#-*- coding: utf-8 -*-
# Copyright (C) 2015, Nikolai Chernikov <nikolai.chernikov.ru@gmail.com>
#
# "convert2ugrid" is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License v3+. "convert2ugrid" is distributed in the
# hope that it will be useful, but WITHOUT ANY WARRANTY.  Consult the file
# LICENSE.GPL or www.gnu.org/licenses/gpl-3.0.txt for the full license terms.



import os
import sys
import re


import make_grid
import process_davit_ncdf
import process_mossco_netcdf
import process_cdl
import process_mixed_data
import ui



def step_1(list_with_synoptic_nc, dictionary_1, dictionary_2, log=False):
    '''
    Script looks through passed MOSSCO netcdf files and searches for those variables
    that can be theoretically converted (valid dimensions); then it cycles over the
    found variables and tries to match them with names from "Dictionary 1". As a
    result of this step "Dictionary 2" is created.

    
    INPUT
        list_with_synoptic_nc - list of strings, paths to netcdf files with synoptic data
        dictionary_1          - string, path to ASCII file with dictionary to suggest
                                standard mossco-baw variable name correlation
    OUTPUT
        dictionary_2          - string, path to ASCII file after scanning variables
    '''
    if log: print ' '*20+'+'*30
    if log: print ' '*20+'+'*30
    if log: print ' '*20+'+'*10+'  STEP 1  '+'+'*10
    if log: print ' '*20+'+'*30
    if log: print ' '*20+'+'*30
    if log: print 'Script is going to scan following NC files:\n{0}\n'.format(list_with_synoptic_nc)
    process_cdl.create_txt_mossco_baw(list_with_synoptic_nc, dictionary_2, dictionary_1, log=log)
    if log: print 'Following file has been created:\n{0}\nIt shows the relation between MOSSCO and DAVIT variable names.'.format(dictionary_2)
    

def step_2(dictionary_2, dictionary_3, dictionary_4, log=False):
    '''
    Scripts goes through "Dictionary 2" and for each variable looks up proper description
    in "Dictionary 3". At the end, joined information is stored in "Dictionary 4".

    INPUT
        dictionary_2          - string, path to ASCII file after scanning variables
        dictionary_3          - string, path to CDL file with standard variables

    OUTPUT
        dictionary_4          - string, path to CDL file to be created
    '''
    if log: print ' '*20+'+'*30
    if log: print ' '*20+'+'*30
    if log: print ' '*20+'+'*10+'  STEP 2  '+'+'*10
    if log: print ' '*20+'+'*30
    if log: print ' '*20+'+'*30
    process_cdl.create_dictionary4(dictionary_3, dictionary_2, dictionary_4, log=log)
    if log: print 'Following file has been created:\n{0}\nIt shows the synoptic data which will be added to NetCDF.'.format(dictionary_4)


def step_3(topo_nc, list_with_synoptic_nc, dictionary_4, nc_out, create_davit_netcdf=True, log=False):
    '''
    INPUT
        topo_nc                 - string, path to netcdf file with x,y vectors and "bathymetry" variable for mask
        list_with_synoptic_nc   - list of strings, paths to netcdf with synoptic data. At first position should be the main data-file
        dictionary_1            - string, path to txt file with dictionary to suggest standart mossco-baw variable name correlation
        dictionary_2            - string, path to txt file after scanning variables
        dictionary_3            - string, path to cdl file with standard variables
        dictionary_4            - string, path to cdl file to be created

    OUTPUT
        nc_out                  - string, path to netcdf to be created
    '''
    if log: print ' '*20+'+'*30
    if log: print ' '*20+'+'*30
    if log: print ' '*20+'+'*10+'  STEP 3  '+'+'*10
    if log: print ' '*20+'+'*30
    if log: print ' '*20+'+'*30
    
    # --------------------------------------------------
    # 1) Read, x any y vectors from the netcdf
    # --------------------------------------------------
    if log: print 'searching x, y, bathymetry ...'
    #coords = process_mossco_netcdf.find_coordinate_vars(topo_nc, log=log)
    structGrid = process_mixed_data.containerForGridGeneration(topo_nc, log=log)
    coords = structGrid.get_data()
    meta   = structGrid.get_metadata()

    # --------------------------------------------------
    # 2) Create mask
    # --------------------------------------------------
    if log: print 'creating mask...'
    #m = process_mossco_netcdf.make_mask_array_from_mossco_bathymetry(topo_nc, varname=coords['bName'], fillvalue=None, transpose=True, log=log)
    m = structGrid.get_mask(transpose=False, log=log)

    

    # --------------------------------------------------
    # 3) Create grid, and unpack values
    # --------------------------------------------------
    start_index = 0
    print 'Generating uGrid... (be patient, this may take a while)'
    dims, topo, nodes, edges, faces, bounds = make_grid.make_2d_qudratic_grid_or_curvilinear(coords['x'], coords['y'],
                                                data_location=meta['x']['points_location'], mask=m, log=log, startingindex=start_index)

    nMesh2_node = dims[0]
    nMesh2_edge = dims[1]
    nMesh2_face = dims[2]
    nMaxMesh2_face_nodes = dims[3]

    Mesh2_edge_nodes = topo[0]
    Mesh2_edge_faces = topo[1]
    Mesh2_face_nodes = topo[2]
    Mesh2_face_edges = topo[3]

    Mesh2_node_x = nodes[0]
    Mesh2_node_y = nodes[1]

    Mesh2_edge_x = edges[0]
    Mesh2_edge_y = edges[1]

    Mesh2_face_x = faces[0]
    Mesh2_face_y = faces[1]
    Mesh2_face_center_x = faces[2]
    Mesh2_face_center_y = faces[3]
    Mesh2_face_area = faces[4]

    Mesh2_edge_x_bnd = bounds[0]
    Mesh2_edge_y_bnd = bounds[1]
    Mesh2_face_x_bnd = bounds[2]
    Mesh2_face_y_bnd = bounds[3]

    # --------------------------------------------------
    # 4) Get bumber of vertical layers, dimensions from MOSSCO
    # --------------------------------------------------
    nLayers, layer_fname = process_mossco_netcdf.get_number_of_depth_layer_from_mossco(list_with_synoptic_nc, dimname='getmGrid3D_getm_3')

    # --------------------------------------------------
    # 5) Create netcdf
    # --------------------------------------------------
    if not create_davit_netcdf:
        print "\nAbortig here...\nNetcdf has not been created, switch flag 'create_davit_netcdf' to True."
        sys.exit(0)

    process_davit_ncdf.create_uGrid_ncdf(nc_out,
                        nMesh2_node, nMesh2_edge, nMesh2_face, nMaxMesh2_face_nodes,
                        Mesh2_edge_nodes, Mesh2_edge_faces, Mesh2_face_nodes, Mesh2_face_edges,
                        Mesh2_node_x, Mesh2_node_y, Mesh2_edge_x, Mesh2_edge_y,
                        Mesh2_face_x, Mesh2_face_y, Mesh2_face_center_x, Mesh2_face_center_y,
                        Mesh2_face_area=Mesh2_face_area,
                        Mesh2_edge_x_bnd=Mesh2_edge_x_bnd, Mesh2_edge_y_bnd=Mesh2_edge_y_bnd,
                        Mesh2_face_x_bnd=Mesh2_face_x_bnd, Mesh2_face_y_bnd=Mesh2_face_y_bnd,
                        coord_type=meta['x']['coordinate_type'],
                        dim_nMesh2_layer2d=1, dim_nMesh2_layer3d=nLayers, dim_nMesh2_class_names_strlen=20, dim_nMesh2_suspension_classes=1,
                        start_index=start_index,
                        log=log)
    if log: print 'grid created', '\n', '-'*100, '\n', 'Now adding data', '\n', '-'*100

    
    # --------------------------------------------------
    # 6) fill netcdf with time variables (TIME and DATATIME)
    # --------------------------------------------------
    if log: print '-'*100
    time_added = False
    for nc_file in list_with_synoptic_nc:
        try:
            process_davit_ncdf.append_Time_andDatetime_to_netcdf(nc_out, nc_file, time_var_name='time', log=log)
            if log: print 'added "nMesh2_data_time" from file ', nc_file
            time_added = True
            break
        except KeyError:  # if var is not found in current file => skip to next file
            pass
    if not time_added:
        print 'WARNING! added dummy "nMesh2_data_time"'
        process_davit_ncdf.append_Time_andDatetime_to_netcdf(nc_out, dummy_values=True, log=log)


    if log: print '-'*100
    # ---------------------------------------------------------------------------
    # 7) Layer thicness
    # ---------------------------------------------------------------------------
    
    # now moved to the very bottom, since we want to use dictionary VARS
    # that is generated by processing Dict4

    # --------------------------------------------------
    # 8) fill netcdf with SYNOPTIC data!!!
    # --------------------------------------------------
    if log: print '-'*100+'\n+Now adding data based on Dictionary 4.\n'+'-'*100
    # --------------------------------------------------
    # -- 8.1) fill netcdf with SYNOPTIC data!!!
    # --------------------------------------------------
    VARS = process_cdl.read_file_with_only_variables(process_cdl.read_file_delete_comments(dictionary_4, comments='//'))

    # -----------------------------------------------------------------------------------------------
    # -- 8.2) cycle through found variables
    # -----------------------------------------------------------------------------------------------
    for var_name, var in VARS.iteritems():
        if log: print 'Variable read from Dictionary 4:', var.get_name()
        varExt = process_mixed_data.cdlVariableExt(var)  #var is an instance of type <cdlVariable>
        source = varExt.get_source_metadata()

        # -----------------------------------------------------------------------------------------------
        # 8.2.4) add data
        # -----------------------------------------------------------------------------------------------

        if source['nNonOneDims'] not in [0, 1, 2, 3, 4]:
            raw_input('WARNING! Skipping variable: <{0}> with dimensions <{1}> of shape <{2}>. It has <{3}> non-one dimensions. Currently <=4 is supported. Press ENTER to continue'.format(
                        source['name'], source['dims'], source['shape']), source['nNonOneDims'])
            break
        raw_data = process_mossco_netcdf.read_mossco_nc_rawvar(source['fname'], var.get_name())
        var_data = process_mixed_data.flatten_xy_data(raw_data, mask=m)


        # -----------------------------------------------------------------------------------------------
        # 8.2.6) append data variable to nc_out...
        # -----------------------------------------------------------------------------------------------
        process_davit_ncdf.append_VariableData_to_netcdf(nc_out, var, var_data, fv=source['fillvalue'],  log=log)
        if log: print '-'*100
        del varExt
    


    # -----------------------------------------------------------------------------------------------
    # 9) Layer thickness
    # -----------------------------------------------------------------------------------------------
    if log: print '-'*100
    if log: print 'Now adding vertical layer information'
    if log: print '-'*100
    vertical_coord_mode = 'sigma'

    if nLayers > 1:  #if a real 3d is here
        try:
            if vertical_coord_mode == 'sigma':
                add_eta   = True
                add_depth = True

                if 'Mesh2_face_Wasserstand_2d' in VARS.keys():
                    add_eta = False
                if 'Mesh2_face_depth_2d' in VARS.keys():
                    add_depth = False
                process_davit_ncdf.append_sigma_vertical_coord_vars(list_with_synoptic_nc, nLayers, nc_out, add_eta=add_eta, add_depth=add_depth, mask=m, log=True)
        except Exception as err:
            try:
                print err
                raw_input('Now i will try to add dummy vertical data. Press ENTER to see info about these values')
                print process_davit_ncdf.append_test_Mesh2_face_z_3d_and_Mesh2_face_z_3d_bnd.__doc__
                if ui.promtYesNo('Do you want to proceed?', quitonno=True):
                    process_davit_ncdf.append_test_Mesh2_face_z_3d_and_Mesh2_face_z_3d_bnd(nc_out, nx=meta['b']['shape'][-1],
                                    ny=meta['b']['shape'][-2], mask=m, log=log)
            except Exception as err:
                print err
                raw_input('Failed to find any vertical-layer information. Will proceed without any. Press ENTER')
                pass

                
    if log: print '-'*100


    if create_davit_netcdf:
        print 'File created:', os.path.abspath(nc_out)



















def create_davit_friendly_netcdf(topo_nc=None, list_with_synoptic_nc=[None], nc_out=None, dictionary_1=None, dictionary_3=None,
                                            start_from_step=1, dictionary_2=None, dictionary_4=None,
                                            create_davit_netcdf=True, log=False, overwrite=False):
    """
        topo_nc                 - string, path to netcdf file with x,y vectors and "bathymetry" variable for mask
        list_with_synoptic_nc   - list of strings, paths to netcdf with synoptic data
        nc_out                  - string, path to netcdf to be created
        dictionary_3            - string, path to cdl file with standard variables
        dictionary_4            - string, path to cdl file to be created, NOW CREATED AUTOMATICALLY
        dictionary_2            - string, path to txt file after scanning variables, NOW CREATED AUTOMATICALLY
        dictionary_1            - string, path to txt file with dictionary to suggest standart mossco-baw variable name correlation

        start_from_step         - integer, (1,2,3) to indicate from which step to start
        overwrite               - True/False , flag to use force overwrite existing files
    """
    
    # It can happen that the user will specify as output parameter an input dictionary
    # In order to protect it (it may be owerwritten), we rename the files...

    nc_out = rename_existing_file(nc_out, force_overwrite=overwrite, log=False)



    # now start program
    if start_from_step == 1:
        if all([topo_nc, list_with_synoptic_nc[0], nc_out, dictionary_1, dictionary_2, dictionary_3, dictionary_4]):
            dictionary_2 = rename_existing_file(dictionary_2, force_overwrite=overwrite, log=False)
            dictionary_4 = rename_existing_file(dictionary_4, force_overwrite=overwrite, log=False)
            
            step_1(list_with_synoptic_nc, dictionary_1, dictionary_2, log=log)
            step_2(dictionary_2, dictionary_3, dictionary_4, log=log)
            step_3(topo_nc, list_with_synoptic_nc, dictionary_4, nc_out, create_davit_netcdf=create_davit_netcdf, log=log)
        else:
            raise ValueError('Fill missing inputs...')
    
    elif start_from_step == 2:
        if all([topo_nc, list_with_synoptic_nc[0], nc_out, dictionary_2, dictionary_3, dictionary_4]):
            dictionary_4 = rename_existing_file(dictionary_4, force_overwrite=overwrite, log=False)
            
            step_2(dictionary_2, dictionary_3, dictionary_4, log=log)
            step_3(topo_nc, list_with_synoptic_nc, dictionary_4, nc_out, create_davit_netcdf=create_davit_netcdf, log=log)
        else:
            raise ValueError('Fill missing inputs...')
    elif start_from_step == 3:
        if all([topo_nc, list_with_synoptic_nc, nc_out, dictionary_4]):
            step_3(topo_nc, list_with_synoptic_nc, dictionary_4, nc_out, create_davit_netcdf=create_davit_netcdf, log=log)
        else:
            raise ValueError('Fill missing inputs...')
    else:
        raise ValueError('Invalid "start_from_step" passed. Choose integer [1, 2, 3]')




def generate_dict_name(dict_name, name_of_the_ncout):
    dirname = os.path.dirname(name_of_the_ncout)
    return os.path.join(dirname, dict_name)




def rename_existing_file(filename, force_overwrite=False, log=False):
    '''
    Input:
        filename - string, path to file (fullpath)
        force_overwrite - True/False, if True, skip this function and return <filename>
    Output:
        new_filename - string with filename to be saved (fullpath)
    ----------------------------------------------------------------------

    the function check whether the passed argument is an existing
    file. If yes, it promts user to confirm overwriting or ...
    generates a new name a new name:
        <oldfilename>(Copy_1).extension
    If copy already existits it will create a new copy (see examples below)
    
    Example 1:
        Lets imagine we have following directory with file test.txt
                ../
                ./
                test.txt
        After running scripts passing fullpath to the file <text.txt>, function
        will return
                path/to/file/test(Copy_1).txt

    Example 2:
        Lets imagine we have following directory with 3 file test.txt
                ../
                ./
                test.txt
                test(Copy_1).txt
                test(Copy_2).txt
        After running scripts passing fullpath to the file <text.txt>, function
        will return
                path/to/file/test(Copy_3).txt
    ------------------------------------------------------------------------
    '''
    _n = 'rename_existing_file():'
    if force_overwrite:
        return filename

    new_filename = filename

    while os.path.isfile(new_filename):
        if ui.promtYesNo('WARNING! File already exists:<{0}> Overwrite?'.format(new_filename)):
            print 'WARNING! I will overwrite the file: <{0}>'.format(filename)
            break
        else:
            dirname = os.path.dirname(new_filename)
            basename = os.path.basename(new_filename)
            raw_name, extension = os.path.splitext(basename)
            if log: print 'filename  :', new_filename
            if log: print 'dirname   :', dirname
            if log: print 'basename  :', basename
            if log: print 'raw_name  :', raw_name
            if log: print 'extension :', extension

            # we dont care about overwriting NetCDF files
            # UPD: since we can overwrite existing data file.... let us always create new names
            #if extension.endswith('nc'):
            #    break

            # check if it is already a copy
            copy_str = re.findall('.*\(Copy\_(\d+)\).*', raw_name)
            
            if copy_str:  # if something was found...
                copy_n = int(copy_str[0])+1  #new number is one more than the old one
                new_raw_name = re.sub('\(Copy\_\d+\)', '(Copy_{0})'.format(copy_n), raw_name)
                new_fname = '{0}{1}'.format(new_raw_name, extension)
            else:  # if no copy exists
                new_fname = '{0}(Copy_1){1}'.format(raw_name, extension)
            
            new_filename = os.path.join(dirname, new_fname)

    if new_filename != filename:
        print _n, 'WARNING! Since <old> file exists, i will create a <new> one and work with it.'
        print _n, '<old>: {0}'.format(filename)
        print _n, '<new>: {0}'.format(new_filename)
    return os.path.abspath(new_filename)







if __name__ == '__main__':

    fname = '/net/widar/home/ak2stud/Nick/python_scripts/dev/uGrid/code/readme.txt'
    print rename_existing_file(fname)
    


