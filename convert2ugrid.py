#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# This file is part of "convert2ugrid" tool
#
# Author: Nikolai Chernikov, nikolai.chernikov.ru@gmail.com
#
# version = 0.3


import os
import sys
import inspect

# use this if you want to include modules from a subfolder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(
                    inspect.getfile( inspect.currentframe() ))[0], "lib")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import create_uGrid_netcdf
import ui

if __name__ == '__main__':
    """
        topo_nc                 - string, path to netcdf file with x,y vectors and "bathymetry" variable for mask
        nc_out                  - string, path to netcdf to be created
        dict1                   - string, path to txt file with dictionary to suggest standart mossco-baw variable name correlation
        dict2                   - string, path to txt file after scanning variables
        dict3                   - string, path to cdl file with standard variables
        dict4                   - string, path to cdl file to be created


        parameter in function create_davit_friendly_netcdf():
        start_from_step         - integer, (1,2,3) to indicate from which step to start (refer to documentation)
                                # start_from_step=1   => start from the very beginning, will create dict2, dict4 and nc_out
                                # start_from_step=2   => start from the second step, will create dict4 and nc_out. Dict2 should be given
                                # start_from_step=3   => start from the third step, will create only nc_out. Dict2 and dict4 should be given
    """
    

    dict1 = os.path.join(os.path.dirname(sys.argv[0]), 'lib', 'defaults', 'dictionary_1')
    dict3 = os.path.join(os.path.dirname(sys.argv[0]), 'lib', 'defaults', 'dictionary_3')
        
    setup_path = '//net/widar/data/nick/to_do/033_UGRID_topo_SNS'
    topo_nc     = os.path.join(setup_path, 'topo.nc')
    #synoptic_nc = os.path.join(setup_path, 'mossco_gfsen.nc')
    #synoptic_nc1 = os.path.join(setup_path, 'shallow_lake-1x1-erosion.3d.0000.nc')
    list_with_synoptic_nc = [topo_nc]

    # setting paths: OUTPUT FILES....
    nc_out = os.path.join(setup_path, 'output/topo_davit1.nc')


    # read command line args!
    params = ui.commandline_support()
    if 'nc_out' in params.keys():
        nc_out = params['nc_out']
    if 'nc_in' in params.keys():
        list_with_synoptic_nc = params['nc_in']
        topo_nc = params['nc_in'][0]
    if 'dict1' in params.keys():
        dict1 = params['dict1']
    if 'dict3' in params.keys():
        dict1 = params['dict3']


    # running script...
    create_uGrid_netcdf.create_davit_friendly_netcdf(topo_nc=topo_nc, list_with_synoptic_nc=list_with_synoptic_nc, nc_out=nc_out,
                    dictionary_1=dict1, dictionary_3=dict3,
                    start_from_step=params['step'], create_davit_netcdf=True, log=params['log'], overwrite=params['overwrite'])

