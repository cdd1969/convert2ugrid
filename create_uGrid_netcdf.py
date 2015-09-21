#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# This file is part of "convert2ugrid" tool
#
# Author: Nikolai Chernikov, nikolai.chernikov.ru@gmail.com
#
# version = 0.1


import os
import sys
import inspect

# use this if you want to include modules from a subfolder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(
                    inspect.getfile( inspect.currentframe() ))[0], "lib")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import make_grid
import process_davit_ncdf
import process_mossco_netcdf
import process_cdl
import process_mixed_data


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
    process_cdl.create_cdl_file(dictionary_3, dictionary_2, dictionary_4)
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
    print 'searching x,y vector...'
    try:
        coord_mode = 'local'
        x_vector, _ = process_mossco_netcdf.read_mossco_nc_1d(topo_nc, 'x')
        y_vector, _ = process_mossco_netcdf.read_mossco_nc_1d(topo_nc, 'y')
    except KeyError:
        coord_mode = 'geographic'
        x_vector, _ = process_mossco_netcdf.read_mossco_nc_1d(topo_nc, 'lon')
        y_vector, _ = process_mossco_netcdf.read_mossco_nc_1d(topo_nc, 'lat')
    print 'working with... {0} coordinates'.format(coord_mode)

    # --------------------------------------------------
    # 2) Create mask
    # --------------------------------------------------
    print 'creating mask...'
    m = process_mossco_netcdf.make_mask_array_from_mossco_bathymetry(topo_nc, varname='bathymetry', fillvalue=None, transpose=True)

    # --------------------------------------------------
    # 3) Create grid, and unpack values
    # --------------------------------------------------
    dims, topo, nodes, edges, faces, bounds = make_grid.make_2d_rectangular_grid(x_vector, y_vector, mask=m, log=True, startingindex=0)

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
    nLayers = process_mossco_netcdf.get_number_of_depth_layer_from_mossco(list_with_synoptic_nc, dimname='getmGrid3D_getm_3')

    # --------------------------------------------------
    # 5) Create netcdf
    # --------------------------------------------------
    if not create_davit_netcdf:
        print "\nExiting...\nNetcdf has not been created, switch flag 'create_davit_netcdf' to True."
        sys.exit(0)

    process_davit_ncdf.create_uGrid_ncdf(nc_out,
                        nMesh2_node, nMesh2_edge, nMesh2_face, nMaxMesh2_face_nodes,
                        Mesh2_edge_nodes, Mesh2_edge_faces, Mesh2_face_nodes, Mesh2_face_edges,
                        Mesh2_node_x, Mesh2_node_y, Mesh2_edge_x, Mesh2_edge_y,
                        Mesh2_face_x, Mesh2_face_y, Mesh2_face_center_x, Mesh2_face_center_y,
                        Mesh2_face_area=Mesh2_face_area,
                        Mesh2_edge_x_bnd=Mesh2_edge_x_bnd, Mesh2_edge_y_bnd=Mesh2_edge_y_bnd,
                        Mesh2_face_x_bnd=Mesh2_face_x_bnd, Mesh2_face_y_bnd=Mesh2_face_y_bnd,
                        coord_mode=coord_mode,
                        dim_nMesh2_layer2d=1, dim_nMesh2_layer3d=nLayers, dim_nMesh2_class_names_strlen=20, dim_nMesh2_suspension_classes=1)
    print 'grid created', '\n', '-'*100, '\n', 'Now adding data', '\n', '-'*100

    
    # --------------------------------------------------
    # 6) fill netcdf with time variables (TIME and DATATIME)
    # --------------------------------------------------
    for nc_file in list_with_synoptic_nc:
        try:
            process_davit_ncdf.append_Time_andDatetime_to_netcdf(nc_out, nc_file, time_var_name='time')
            print 'added "nMesh2_data_time" from file ', nc_file
        except KeyError:  # if var is not found in current file => skip to next file
            pass


    # ---------------------------------------------------------------------------
    # 7) Layer thicness
    # ---------------------------------------------------------------------------
    
    if nLayers > 1:  #if a real 3d is here
        #pass
        process_davit_ncdf.append_test_Mesh2_face_z_3d_and_Mesh2_face_z_3d_bnd(nc_out, list_with_synoptic_nc[0], mask=m, log=False)


    # --------------------------------------------------
    # 8) fill netcdf with SYNOPTIC data!!!
    # --------------------------------------------------
    print '-'*100
    print 'Now adding data based on CDL file'
    print '-'*100
    # --------------------------------------------------
    # -- 8.1) fill netcdf with SYNOPTIC data!!!
    # --------------------------------------------------
    VARS = process_cdl.read_file_with_only_variables(process_cdl.read_file_delete_comments(dictionary_4, comments='//'))

    # -----------------------------------------------------------------------------------------------
    # -- 8.2) cycle through found variables
    # -----------------------------------------------------------------------------------------------
    for VN, VV in VARS.iteritems():
        print 'Working with variable read from CDL:', VN

        # -----------------------------------------------------------------------------------------------
        # -- 8.2.1) add basic info
        # -----------------------------------------------------------------------------------------------
        var_to_add = dict()
        var_to_add['vname'] = VN
        if VV[0] == 'float' : VV[0] = 'f4'
        if VV[0] == 'double': VV[0] = 'f8'
        var_to_add['dtype'] = VV[0]  # here there might be a problem, cause datatypes are specified in CDL format (i.e. "double", "float", "int", etc.)
        var_to_add['dims'] = tuple(VV[1])

        # -----------------------------------------------------------------------------------------------
        # -- 8.2.2) adding fillvalue
        # -----------------------------------------------------------------------------------------------
        if '_FillValue' in VV[2].keys():
            var_to_add['_FillValue'] = VV[2]['_FillValue']
        else:
            var_to_add['_FillValue'] = False
        
        # -----------------------------------------------------------------------------------------------
        # -- 8.2.3) adding attributes
        # -----------------------------------------------------------------------------------------------
        attrs = dict()
        for ATR_N, ATR_V in VV[2].iteritems():
            if not ATR_N.startswith('_'):
                attrs[ATR_N] = ATR_V
        var_to_add['attributes'] = attrs

        # -----------------------------------------------------------------------------------------------
        # -- 8.2.4) check if we have data to add...
        # -----------------------------------------------------------------------------------------------
        if '_mossco_filename' in VV[2].keys():
            fname = VV[2]['_mossco_filename']
            varname = VV[2]['_mossco_varname']
        else:
            fname, varname, var_to_add['data'] = None, None, None
            # -----------------------------------------------------------------------------------------------
            # Create auto variables
            # -----------------------------------------------------------------------------------------------
            if '_auto_creation' in VV[2].keys():
                var_to_add['data'] = process_mixed_data.create_magnitude_variable_from_x_y_component(VARS,
                                        VN, VV, mask=m, log=True)

        # -----------------------------------------------------------------------------------------------
        # 8.2.5) add data
        # -----------------------------------------------------------------------------------------------
        if fname:
            # if 0D
            if var_to_add['dims'] == tuple([]):
                var_to_add['data'] = None
            
            # if 1D
            elif var_to_add['dims'] in [tuple(['nMesh2_data_time'])
                                        ]:
                var_to_add['data'], dim_shape = process_mossco_netcdf.read_mossco_nc_1d(fname, varname)
            
            # if 2D
            elif var_to_add['dims'] in [tuple(['nMesh2_time', 'nMesh2_face']),
                                        tuple(['nMesh2_time', 'nMesh2_layer_2d', 'nMesh2_face']),
                                        tuple(['nMesh2_time', 'nMesh2_suspension_classes', 'nMesh2_layer_2d', 'nMesh2_face'])
                                        ]:
                var_to_add['data'], dim_shape = process_mossco_netcdf.read_mossco_nc_2d(fname, varname, mask=m)
            
            # if 3D
            elif var_to_add['dims'] in [tuple(['nMesh2_data_time', 'nMesh2_face']),
                                        tuple(['nMesh2_data_time', 'nMesh2_layer_2d', 'nMesh2_face']),
                                        tuple(['nMesh2_data_time', 'nMesh2_suspension_classes', 'nMesh2_layer_2d', 'nMesh2_face']),
                                        tuple(['nMesh2_time', 'nMesh2_layer_3d', 'nMesh2_face']),
                                        tuple(['nMesh2_time', 'nMesh2_suspension_classes', 'nMesh2_layer_3d', 'nMesh2_face'])
                                        ]:
                var_to_add['data'], dim_shape = process_mossco_netcdf.read_mossco_nc_3d(fname, varname, mask=m)

            # if 4D
            elif var_to_add['dims'] in [tuple(['nMesh2_data_time', 'nMesh2_layer_3d', 'nMesh2_face']),
                                        tuple(['nMesh2_data_time', 'nMesh2_suspension_classes', 'nMesh2_layer_3d', 'nMesh2_face'])
                                        ]:
                var_to_add['data'], dim_shape = process_mossco_netcdf.read_mossco_nc_4d(fname, varname, mask=m)

            else:
                raise KeyError('Skipping variable: {1}\nDimensions "{0}" not recognised.\nCheck for hardcoded solution'.format(var_to_add['dims'], VN) )
                print 'Skipping variable: {1}\nDimensions "{0}" not recognised.\nCheck for hardcoded solution'.format(var_to_add['dims'], VN)
                break
        # -----------------------------------------------------------------------------------------------
        # 8.2.6) append data variable to nc_out...
        # -----------------------------------------------------------------------------------------------
        process_davit_ncdf.append_VariableData_to_netcdf(nc_out, var_to_add)
        print '-'*100
    
    if create_davit_netcdf and log:
        print '\n\nFile created:', os.path.abspath(nc_out)






def create_davit_friendly_netcdf(topo_nc=None, list_with_synoptic_nc=None, nc_out=None, dictionary_1=None, dictionary_2=None, dictionary_3=None,
                                            dictionary_4=None, start_from_step=1,
                                            create_davit_netcdf=True, log=False):
    """
        topo_nc                 - string, path to netcdf file with x,y vectors and "bathymetry" variable for mask
        list_with_synoptic_nc   - list of strings, paths to netcdf with synoptic data
        nc_out                  - string, path to netcdf to be created
        dictionary_3            - string, path to cdl file with standard variables
        dictionary_4            - string, path to cdl file to be created
        dictionary_2            - string, path to txt file after scanning variables
        dictionary_1            - string, path to txt file with dictionary to suggest standart mossco-baw variable name correlation

        start_from_step         - integer, (1,2,3) to indicate from which step to start
    """
    if start_from_step == 1:
        if all([topo_nc, list_with_synoptic_nc, nc_out, dictionary_1, dictionary_2, dictionary_3, dictionary_4]):
            step_1(list_with_synoptic_nc, dictionary_1, dictionary_2, log=log)
            step_2(dictionary_2, dictionary_3, dictionary_4, log=log)
            step_3(topo_nc, list_with_synoptic_nc, dictionary_4, nc_out, create_davit_netcdf=create_davit_netcdf, log=log)
        else:
            raise ValueError('Fill missing inputs...')
    
    elif start_from_step == 2:
        if all([topo_nc, list_with_synoptic_nc, nc_out, dictionary_2, dictionary_3, dictionary_4]):
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











if __name__ == '__main__':

    # setting paths: INPUT FILES...
    currentPath = os.path.dirname(sys.argv[0])
    dict1 = os.path.join(currentPath, 'user_input/dictionary1.txt')  # see description of dictionaries in documentation
    dict3 = os.path.join(currentPath, 'user_input/dictionary3.cdl')  # see description of dictionaries in documentation
    
    setup_path = '/net/widar/home/ak2stud/Nick/python_scripts/dev/uGrid/data/NSBS'
    
    topo_nc     = os.path.join(setup_path, 'topo.nc')                 #topo-file with bathymetry and grid
    synoptic_nc = os.path.join(setup_path, 'netcdf_reference_3d.nc')  #netcdf with simulation data (may be more than one, join them in a list below)
    list_with_synoptic_nc = [synoptic_nc, topo_nc]                                 #join files in a list

    # setting paths: OUTPUT FILES....
    dict2  = os.path.join(currentPath, '../data/NSBS/out/tmp', 'dictionary2.txt')  # see description of dictionaries in documentation
    dict4  = os.path.join(currentPath, '../data/NSBS/out/tmp', 'dictionary4.cdl')  # see description of dictionaries in documentation
    nc_out = os.path.join(currentPath, '../data/NSBS/out/tmp', 'nsbs_davit.nc')    # file to be created

    # running script...
    create_davit_friendly_netcdf(topo_nc=topo_nc, list_with_synoptic_nc=list_with_synoptic_nc, nc_out=nc_out,
                    dictionary_1=dict1, dictionary_2=dict2, dictionary_3=dict3, dictionary_4=dict4,
                    start_from_step=1, create_davit_netcdf=True, log=True)
