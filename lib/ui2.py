import click
import os
import sys
import netCDF4

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

def commandline_support(nc_in, nc_out, topo, step, d1, d2, d3, d4, force, mute):
    '''\nconvert2ugrid
    Script can convert structured grid netcdf files to unstructured grid format. Output files then can be viewed with DAVIT, QUICKPLOT or NCVIEW
    
    \b
    Example usage:
      $ python convert2ugrid.py [-s n] [--dict1 s] [--dict2 s] [--dict3 s] [--dict4 s] [-t s] [-o s] [-f] [-m] nc_in1 nc_in2
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
    #click.echo(click.style('Finished: <{0}> created successfully.'.format(fname_out), fg='green', bold=True))
    P = dict()
    P['log'] = not mute
    P['step'] = step
    P['overwrite'] = force
    P['nc_in']  = [topo]+list(nc_in)
    P['dict1']  = d1
    P['dict2']  = d2
    P['dict3']  = d3
    P['dict4']  = d4
    P['nc_out'] = nc_out

    print 'P', P
    return click.echo(P)


def promtYesNo(question='', quitonno=False):
    '''Promt a yes-no question in console and catch answer

    Inputs:
        question    -  string, will show this text before [y/n].
        quitonno    -  True/False. If True will force exit application if "No" answered

    Output:
        True        - if answer is Yes
        False       - if answer is No
        sys.exit(2) - if answer is No and quitonno=True
    '''
    answer = click.confirm(click.style(question, fg='yellow'), abort=quitonno)
    return answer


def promt(*args, **kwargs):
    if 'pause' in kwargs and kwargs['pause'] is True:
        click.pause()
        return
    if 'color' in kwargs:
        color = kwargs.pop('color')
        return click.prompt(click.style(args[0], fg=color), args[1::], **kwargs)
    else:
        return click.prompt(*args, **kwargs)

if __name__ == '__main__':
    commandline_support()
