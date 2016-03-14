#!/usr/bin/python
#-*- coding: utf-8 -*-
# Copyright (C) 2015, Nikolai Chernikov <nikolai.chernikov.ru@gmail.com>
#
# "convert2ugrid" is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License v3+. "convert2ugrid" is distributed in the
# hope that it will be useful, but WITHOUT ANY WARRANTY.  Consult the file
# LICENSE.GPL or www.gnu.org/licenses/gpl-3.0.txt for the full license terms.

__author__  = 'Nikolai Chernikov <nikolai.chernikov.ru@gmail.com>'

import os
import sys
import netCDF4
import lib.create_uGrid_netcdf
from lib import CLICK_IS_HERE
if CLICK_IS_HERE:
    import click
else:
    import lib.ui as ui


def main_old():
    """
    Script can convert structured grid netcdf files to unstructured grid format. Output files then can be viewed with DAVIT, QUICKPLOT or NCVIEW
    Example usage (from console):

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


    # read command line args, and set defaults
    params = ui.commandline_support()

    # running script...
    lib.create_uGrid_netcdf.create_davit_friendly_netcdf(topo_nc=params['nc_in'][0], list_with_synoptic_nc=params['nc_in'], nc_out=params['nc_out'],
        dictionary_1=params['dict1'], dictionary_3=params['dict3'], dictionary_2=params['dict2'], dictionary_4=params['dict4'],
        start_from_step=params['step'], create_davit_netcdf=True, log=params['log'], overwrite=params['overwrite'])














if __name__ == '__main__':
    if CLICK_IS_HERE:
        CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
        @click.command(context_settings=CONTEXT_SETTINGS)
        #@click.group(chain=True, invoke_without_command=True)
        @click.argument('nc_in', nargs=-1
            )
        @click.option('--nc_out', '-o', type=click.Path(exists=False, dir_okay=False), default='ugrid.nc',
                        help='Name of the output netcdf file(-s) with results. Default: `ugrid.nc`'
            )
        @click.option('--topo', '-t', type=click.Path(exists=True, dir_okay=False), default='topo.nc',
                        help='Name of input topo-netcdf file with coordinates and bathymetry. This file will be used for ugrid-generation. Default: `topo.nc`'
            )
        @click.option('--step', '-s', type=click.IntRange(min=1, max=3), default=1,
            help='Start script execution from this step (see DOCS) Default: 1.\n 1 => create dict2, dict4, nc_out\n 2 => create dict4, nc_out; dict2 must be already created\n 3 => create only nc_out; dict2 and dict4 must be already created'
            )
        @click.option('--dict1', '--d1', type=click.Path(exists=True, dir_okay=False), default=os.path.join(os.path.dirname(sys.argv[0]), 'lib', 'defaults', 'dictionary_1'),
            help='Name of the input Dictionary1. Default: `./lib/defaults/dictionary_1`'
            )
        @click.option('--dict2', '--d2', type=click.Path(exists=False, dir_okay=False), default='dictionary_2',
            help='Name of the input Dictionary2. Default: `dictionary_2`'
            )
        @click.option('--dict3', '--d3', type=click.Path(exists=True, dir_okay=False), default=os.path.join(os.path.dirname(sys.argv[0]), 'lib', 'defaults', 'dictionary_3'),
            help='Name of the input Dictionary3. Default: `./lib/defaults/dictionary_3`'
            )
        @click.option('--dict4', '--d4', type=click.Path(exists=False, dir_okay=False), default='dictionary_4',
            help='Name of the input Dictionary4. Default: `dictionary_4`'
            )
        @click.option('--force', '-f', is_flag=True, default=False,
            help='Flag to force overwrite. Do not promt user before overwriting existing files. Default: False'
            )
        @click.option('--mute', '-m', is_flag=True, default=False,
            help='Flag to run in silent mode; minimum logs. Default: False'
            )
        @click.option('--hardcoded', '-c', is_flag=True, default=False,
            help='Flag to use input params hard-coded in file `convert2ugrid.py`. Toggling this flag will IGNORE ALL other command-line inputs. Default: False'
            )
        def main_new(nc_in, nc_out, topo, step, d1, d2, d3, d4, force, mute, hardcoded):
            '''\nconvert2ugrid
            Script can convert structured grid netcdf files to unstructured grid format. Output files then can be viewed with DAVIT, QUICKPLOT or NCVIEW
            
            \b
            Example usage:
              $ python convert2ugrid.py -c
               >>> will use input parameters hard-coded in the file <convert2ugrid.py>\n
              $ python convert2ugrid.py [-s n] [--dict1 s] [--dict2 s] [--dict3 s] [--dict4 s] [-t s] [-o s] [-f] [-m] [-h] nc_in1 nc_in2
               >>> see description below'
            '''
            try:
                nc = netCDF4.Dataset(topo , mode='r')
                nc.close()
            except Exception, err:
                raise click.BadParameter('( {1} ) Can not read NetCDF file {0}'.format(topo, err), param_hint=['topo'])

            for fname_in in nc_in:
                try:
                    nc = netCDF4.Dataset(fname_in , mode='r')
                    nc.close()
                except Exception, err:
                    raise click.BadParameter('( {1} ) Can not read NetCDF file {0}'.format(fname_in, err), param_hint=['nc_in'])
            if d2 != 'dictionary_2':
                if not os.path.isfile(d2):
                    raise click.BadParameter('File {0} does not exist'.format(d2), param_hint=['d2'])
            if d4 != 'dictionary_4':
                if not os.path.isfile(d4):
                    raise click.BadParameter('File {0} does not exist'.format(d4), param_hint=['dict4'])
            #click.echo(click.style('Processing: <{0}>'.format(fname_in), fg='yellow', bold=True))

            lib.create_uGrid_netcdf.create_davit_friendly_netcdf(topo_nc=topo, list_with_synoptic_nc=[topo]+list(nc_in), nc_out=nc_out,
                dictionary_1=d1, dictionary_3=d3, dictionary_2=d2, dictionary_4=d4,
                start_from_step=step, create_davit_netcdf=True, log=not mute, overwrite=force)
            
            click.secho('Finished: <{0}> created successfully.'.format(nc_out), fg='green', bold=True)
        

        main_new()
    else:
        main_old()
