#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# This file is part of "convert2ugrid" tool
#
# Author: Nikolai Chernikov, nikolai.chernikov.ru@gmail.com
#
# version = 0.2


import os
import sys
import inspect

# use this if you want to include modules from a subfolder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(
                    inspect.getfile( inspect.currentframe() ))[0], "lib")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import create_uGrid_netcdf


if __name__ == '__main__':

    # setting paths: INPUT FILES...
    currentPath = os.path.dirname(sys.argv[0])
    dict1 = os.path.join(currentPath, 'user_input/dictionary1.txt')  # see description of dictionaries in documentation
    dict3 = os.path.join(currentPath, 'user_input/dictionary3.cdl')  # see description of dictionaries in documentation
    
    setup_path = '/net/widar/home/ak2stud/Nick/python_scripts/dev/uGrid/data/NSBS'
    
    topo_nc     = os.path.join(setup_path, 'topo.nc')                 #topo-file with bathymetry and grid
    synoptic_nc = os.path.join(setup_path, 'netcdf_reference_3d.nc')  #netcdf with simulation data (may be more than one, join them in a list below)
    list_with_synoptic_nc = [synoptic_nc, topo_nc]                    #join files in a list

    # setting paths: OUTPUT FILES....
    dict2  = os.path.join(currentPath, '../data/NSBS/out/tmp', 'dictionary2.txt')  # see description of dictionaries in documentation
    dict4  = os.path.join(currentPath, '../data/NSBS/out/tmp', 'dictionary4.cdl')  # see description of dictionaries in documentation
    nc_out = os.path.join(currentPath, '../data/NSBS/out/tmp', 'nsbs_davit.nc')    # file to be created

    # running script...
    create_uGrid_netcdf.create_davit_friendly_netcdf(topo_nc=topo_nc, list_with_synoptic_nc=list_with_synoptic_nc, nc_out=nc_out,
                    dictionary_1=dict1, dictionary_2=dict2, dictionary_3=dict3, dictionary_4=dict4,
                    start_from_step=1, create_davit_netcdf=True, log=True)
