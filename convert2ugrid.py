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
    setup_path = './info/tutorial'
    dict1 = os.path.join(setup_path, 'user_input/dictionary1.txt')
    dict3 = os.path.join(setup_path, 'user_input/dictionary3.cdl')
        
    topo_nc     = os.path.join(setup_path, 'data/topo.nc')
    synoptic_nc = os.path.join(setup_path, 'data/netcdf_reference_3d.nc')
    list_with_synoptic_nc = [synoptic_nc, topo_nc]

    # setting paths: OUTPUT FILES....
    dict2  = os.path.join(setup_path,  'output/dictionary2.txt')
    dict4  = os.path.join(setup_path, 'output/dictionary4.cdl')
    nc_out = os.path.join(setup_path, 'output/nsbs_davit_1.nc')

    # running script...
    create_uGrid_netcdf.create_davit_friendly_netcdf(topo_nc=topo_nc, list_with_synoptic_nc=list_with_synoptic_nc, nc_out=nc_out,
                    dictionary_1=dict1, dictionary_2=dict2, dictionary_3=dict3, dictionary_4=dict4,
                    start_from_step=1, create_davit_netcdf=True, log=True)

