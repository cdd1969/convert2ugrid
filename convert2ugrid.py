#!/usr/bin/python
#-*- coding: utf-8 -*-
# Copyright (C) 2015, Nikolai Chernikov <nikolai.chernikov.ru@gmail.com>
#
# "convert2ugrid" is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License v3+. "convert2ugrid" is distributed in the
# hope that it will be useful, but WITHOUT ANY WARRANTY.  Consult the file
# LICENSE.GPL or www.gnu.org/licenses/gpl-3.0.txt for the full license terms.

__author__  = 'Nikolai Chernikov <nikolai.chernikov.ru@gmail.com>'
__version__ = '0.3.1'


import os
#try:
#    import click
#    import lib.ui2 as ui
#except:
#    import lib.ui as ui
import lib.ui as ui  # this line should be deleted after ui2 will be tested (see imports above)
import lib.create_uGrid_netcdf


if __name__ == '__main__':
    """
    Script can convert structured grid netcdf files to unstructured grid format. Output files then can be viewed with DAVIT, QUICKPLOT or NCVIEW
    Example usage (from console):
     
     (variant 1)$python convert2ugrid.py -c
     (variant 1) >>> will use input parameters specified in the file <convert2ugrid.py>
     
     (variant 2)$python convert2ugrid.py [-s n] [-i s1,s2,...] [--dict1=s] [--dict3=s] [-o s] [-f] [-m] [-h]
     (variant 2) >>> see description below

        The inputs hard-coded here will only be used if you run the script with (variant 1), otherwise
        they will be ignored

    Args:
    -----
        topo_nc (str):
            path to netcdf file with x,y vectors and "bathymetry" variable for mask
        nc_out (str):
            path to netcdf to be created
        dict1 (str):
            path to txt file with dictionary to suggest standart mossco-baw variable name correlation
        dict2 (str):
            path to txt file after scanning variables
        dict3 (str):
            path to cdl file with standard variables
        dict4 (str):
            path to cdl file to be created


        parameter in function create_davit_friendly_netcdf():
        start_from_step         - integer, (1,2,3) to indicate from which step to start (refer to documentation)
                                # start_from_step=1   => start from the very beginning, will create dict2, dict4 and nc_out
                                # start_from_step=2   => start from the second step, will create dict4 and nc_out. Dict2 should be given
                                # start_from_step=3   => start from the third step, will create only nc_out. Dict2 and dict4 should be given
    """
    # ---------------------
    #     USER SETTINGS
    #  (only for variant1)
    # ---------------------
    show_more_logs  = True
    start_from_step = 1
    setup_path = '//net/widar/data/nick/to_do/033_UGRID_topo_SNS'
    
    topo_nc     = os.path.join(setup_path, 'topo.nc')
    synoptic_nc = os.path.join(setup_path, 'mossco_gfsen.nc')
    synoptic_nc1 = os.path.join(setup_path, 'shallow_lake-1x1-erosion.3d.0000.nc')
    list_with_synoptic_nc = [topo_nc, synoptic_nc, synoptic_nc1]
    

    dict1 = os.path.join(setup_path, 'user_input/dictionary1.txt')
    dict3 = os.path.join(setup_path, 'user_input/dictionary3.cdl')

    # setting paths: OUTPUT FILES....
    dict2  = os.path.join(setup_path, 'output/dictionary2.txt')
    dict4  = os.path.join(setup_path, 'output/dictionary4.cdl')
    nc_out = os.path.join(setup_path, 'output/nsbs_davit_1.nc')
    # ---------------------
    # END USER SETTINGS END
    #  (only for variant1)
    # ---------------------

    # read command line args, and set defaults
    params = ui.commandline_support()


    # running script...
    if params['use_code_inputs']:
        lib.create_uGrid_netcdf.create_davit_friendly_netcdf(topo_nc=topo_nc, list_with_synoptic_nc=list_with_synoptic_nc, nc_out=nc_out,
                    dictionary_1=dict1, dictionary_2=dict2, dictionary_3=dict3, dictionary_4=dict4,
                    start_from_step=start_from_step, create_davit_netcdf=True, log=show_more_logs)
    else:
        lib.create_uGrid_netcdf.create_davit_friendly_netcdf(topo_nc=params['nc_in'][0], list_with_synoptic_nc=params['nc_in'], nc_out=params['nc_out'],
                    dictionary_1=params['dict1'], dictionary_3=params['dict3'], dictionary_2=params['dict2'], dictionary_4=params['dict4'],
                    start_from_step=params['step'], create_davit_netcdf=True, log=params['log'], overwrite=params['overwrite'])
