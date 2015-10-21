#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# This file is part of "convert2ugrid" tool
#
# Author: Nikolai Chernikov, nikolai.chernikov.ru@gmail.com

import sys
import os
import getopt


class FileError(Exception):
    """Exception raised if file does nto exist.

    Attributes:
        filenamre -- string, abspath of the file
    """
    def __init__(self, filename):
        self.fn = filename
    def __str__(self):
        print ': File not found: {0}'.format(repr(self.fn))


class InputError(Exception):
    """Exception raised for errors in the input.

    Attributes:
        expr -- input expression in which the error occurred
        msg  -- explanation of the error
    """

    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg
    def __str__(self):
        print ': '+self.msg
        print 'You have given: <{0}>'.format(self.expr)
        print 'For help type: $python convert2ugrid -h'


def usage():
    print '\n\n\nconvert2ugrid'
    print 'Script can convert structured grid netcdf files to unstructured grid format. Output files then can be viewed with DAVIT, QUICKPLOT or NCVIEW'
    print 'Example usage:'
    print ' (variant 1)$python convert2ugrid.py -c'
    print ' (variant 1) >>> will use input parameters hard-coded in the file <convert2ugrid.py>\n'
    print ' (variant 2)$python convert2ugrid.py [-s n] [-i s1,s2,...] [--dict1=s] [--dict3=s] [-o s] [-f] [-m] [-h]'
    print ' (variant 2) >>> see description below'
    print '*'*120
    print ' [-c; --code]              >  Flag to use input params hard-coded in file <convert2ugrid.py>.'
    print '                              Toggling this flag will IGNORE ALL other command-line inputs.'
    print ''
    print ' [-s n; --step=n]          >  Start from step <n>, integer (refer to documentation).'
    print '                              If not set, [-s 1] is used.'
    print '                              n=1 => create dict2, dict4 and netcdf.'
    print '                              n=2 => create dict4 and netcdf, dict2 should be already created.'
    print '                              n=3 => create only netcdf, dict2 and dict4 should be already created.'
    print ''
    print ' [-i s1,s2; --input=s1,s2] >  <s1,s2,...> name(-s) of the input file(-s), comma-separated string.'
    print '                              First file should always be topo-file!'
    print '                              Example reading 1 file : [-i topo.nc]'
    print '                              Example reading 3 files: [-i topo.nc,data1.nc,data2.nc]'
    print ''
    print ' [-o s; --output=s]        >  <s> is the name of the outputfile to be created, string.'
    print '                              If not set, will create file "ugrid.nc" in the current working directory.'
    print ''
    print ' [--dict1=s]               >  <s> name of the input dictionary1, string.'
    print '                              If not set, will read from <./lib/defaults/dictionary_1>.'
    print ''
    print ' [--dict2=s]               >  <s> name of the input/output dictionary2, string.'
    print '                              If not set, will read/write within directory of the output netcdf file (see -o)'
    print '                              under default name "dictionary_2". Read/write status depends on the step (see -s).'
    print '                              This file is only read when [-s 2] or [-s 3], and is created with [-s 1].'
    print ''
    print ' [--dict3=s]               >  <s> name of the input dictionary3, string.'
    print '                              If not set, will read from <./lib/defaults/dictionary_3>.'
    print ''
    print ' [--dict4=s]               >  <s> name of the input/output dictionary4, string.'
    print '                              If not set, will read/write within directory of the output netcdf file (see -o)'
    print '                              under default name "dictionary_4". Read/write status depends on the step (see -s).'
    print '                              This file is only read when [-s 3], and is created with [-s 1] or [-s 2].'
    print ''
    print ' [-f; --force]             >  Flag to force overwrite. Do not promt user before overwriting existing files.'
    print ''
    print ' [-m; --mute]              >  Flag to run in silent mode; minimum logs.'
    print ''
    print ' [-h; --help]              >  Show help, exit.'
    print '*'*120





def check_inputs(params):
    # using hard-coded inputs
    if params['use_code_inputs'] is True:
        print '\nWARNING! [-c] argument passed. I ignore all the command-line inputs. Use those hard-coded in file <convert2ugrid.py>'
        return
    # using command-line inputs
    else:
        if params['nc_in'][0] is None:
            print 'No inputs for netcdf files. Pass at least one input file with [-i n]. For help parse [-h]'
            sys.exit(2)

        try:
            params['step'] = int(params['step'])
        except:
            msg = 'Parameter [-s; --step] should be integer. It may have one of the following values [1, 2, 3]'
            raise InputError(params['step'], msg)
            
        if params['step'] not in [1, 2, 3]:
            msg = 'Parameter [-s; --step] should be integer. It may have one of the following values [1, 2, 3]'
            raise InputError(params['step'], msg)

        
        else:
            for fn in params['nc_in']:
                if not os.path.isfile(fn):
                    raise FileError(fn)
            for pn in ['dict1', 'dict2', 'dict3', 'dict4']:
                if (pn == 'dict2' and params['dict2_default'] is True) or (
                    pn == 'dict4' and params['dict4_default'] is True):
                    # dict has not been passed. It may not exist yet, so we use default names
                    continue

                if not os.path.isfile(params[pn]):
                    raise FileError(params[pn])
        # now we know all the files exist!


def commandline_support():
    P = dict()
    # defaults....
    P['log'] = True
    P['step'] = 1
    P['overwrite'] = False
    P['nc_in']  = [None]
    P['dict1']  = os.path.join(os.path.dirname(sys.argv[0]), 'lib', 'defaults', 'dictionary_1')
    P['dict2']  = os.path.join(os.getcwd(), 'dictionary_2')
    P['dict2_default']  = True
    P['dict3']  = os.path.join(os.path.dirname(sys.argv[0]), 'lib', 'defaults', 'dictionary_3')
    P['dict4']  = os.path.join(os.getcwd(), 'dictionary_4')
    P['dict4_default']  = True
    P['nc_out'] = os.path.join(os.getcwd(), 'ugrid.nc')
    P['use_code_inputs'] = False
    # dont do anything if running without params
    if len(sys.argv) == 1:
        usage()
        sys.exit(2)
    else:
        # if more than 1 arg has been passed...  here the fun starts
        try:
            opts, args = getopt.getopt(sys.argv[1:], "s:o:i:fmhc",
                ['code', 'step=', 'input=', 'dict1=', 'dict2=', 'dict3=', 'dict4=' 'output=', 'force', "help", "mute", ])
        except getopt.GetoptError as err:
            print str(err)
            usage()
            sys.exit(2)

        for opt, arg in opts:
            if opt in ['-m', '--mute']:
                P['log'] = False
            elif opt in ['-s', '--step']:
                P['step'] = arg
            elif opt in ['-i', '--input']:
                P['nc_in'] = str(arg).split(',')
                for i, nc_name in enumerate(P['nc_in']):
                    P['nc_in'][i] = os.path.abspath(nc_name)
            elif opt in ['--dict1']:
                P['dict1'] = os.path.abspath(str(arg))
            elif opt in ['--dict2']:
                P['dict2'] = os.path.abspath(str(arg))
                P['dict2_default']  = False
            elif opt in ['--dict3']:
                P['dict3'] = os.path.abspath(str(arg))
            elif opt in ['--dict4']:
                P['dict4'] = os.path.abspath(str(arg))
                P['dict4_default']  = False
            elif opt in ['-o', '--output']:
                P['nc_out'] = os.path.abspath(str(arg))
            elif opt in ['-f', '--force']:
                P['overwrite'] = True
            elif opt in ['-c', '--code']:
                P['use_code_inputs'] = True
            elif opt in ['-h', '--help']:
                usage()
                sys.exit(2)

    check_inputs(P)
    return P
