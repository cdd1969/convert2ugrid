#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# This file is part of "convert2ugrid" tool
#
# Author: Nikolai Chernikov, nikolai.chernikov.ru@gmail.com

import sys
import os
import getopt

try:
    # This will modify the behavior of raw_input() so that it behaves more like the python interactive shell
    # in terms of history and line editing.
    import readline
except:
    pass #readline not available



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
        print 'You have typed: <{0}>'.format(self.expr)
        print 'For help type: $python convert2ugrid.py -h'


def usage():
    print '\n\n\nconvert2ugrid'
    print 'Script can convert structured grid netcdf files to unstructured grid format. Output files then can be viewed with DAVIT, QUICKPLOT or NCVIEW'
    print 'Example usage:'
    print ' (variant 1)$python convert2ugrid.py -c'
    print ' (variant 1) >>> will use input parameters hard-coded in the file <convert2ugrid.py>\n'
    print ' (variant 2)$python convert2ugrid.py [-s n] [--dict1=s] [--dict2=s] [--dict3=s] [--dict4=s] [-t s] [-o s] [-f] [-m] [-h] s1'
    print ' (variant 2) >>> see description below'
    print '*'*120
    print ' [-c; --code]   >  (OPTIONAL) Flag to use input params hard-coded in file <convert2ugrid.py>.'
    print '                   Toggling this flag will IGNORE ALL other command-line inputs.'
    print '                                              '
    print ' [-t s]         >  (MANDATORY) <s> is the name of the input topo-netcdf file with coordinates and bathymetry, string.'
    print '                   This file will be used for ugrid-generation'
    print '                                              '
    print ' s1             >  (OPTIONAL) Is the argument following the other opts. String. Should be at the end'
    print '                   Names of the netcdf files to be converted. Note: you dont have to write here name'
    print '                   of the mandatory topo-file specified with [-t s], since it will be included anyway.'
    print '                   You may parse unlimited number of these arguments <s1 s2 s3 ...>'
    print '                   Example: <$convert2ugrid.py -t topo.nc -o newfile.nc file1.nc file2.nc file3.nc>'
    print '                            <$convert2ugrid.py -t topo.nc -o newfile.nc>'
    print '                                              '
    print ' [-s n]         >  (OPTIONAL) Start from step <n>, integer (refer to documentation).'
    print '                   If not set, [-s 1] is used.'
    print '                   n=1 => create dict2, dict4 and netcdf.'
    print '                   n=2 => create dict4 and netcdf, dict2 should be already created.'
    print '                   n=3 => create only netcdf, dict2 and dict4 should be already created.'
    print '                                              '
    print ' [-o s]         >  (OPTIONAL) <s> is the name of the outputfile to be created, string.'
    print '                   If not set, will create file "ugrid.nc" in the current working directory.'
    print '                                              '
    print ' [--dict1=s]    >  (OPTIONAL) <s> name of the input dictionary1, string.'
    print '                   If not set, will read from <./lib/defaults/dictionary_1>.'
    print '                                              '
    print ' [--dict2=s]    >  (OPTIONAL) <s> name of the input/output dictionary2, string.'
    print '                   If not set, will read/write within directory of the output netcdf file (see -o)'
    print '                   under default name "dictionary_2". Read/write status depends on the step (see -s).'
    print '                   This file is only read when [-s 2] or [-s 3], and is created with [-s 1].'
    print '                                              '
    print ' [--dict3=s]    >  (OPTIONAL) <s> name of the input dictionary3, string.'
    print '                   If not set, will read from <./lib/defaults/dictionary_3>.'
    print '                                              '
    print ' [--dict4=s]    >  (OPTIONAL) <s> name of the input/output dictionary4, string.'
    print '                   If not set, will read/write within directory of the output netcdf file (see -o)'
    print '                   under default name "dictionary_4". Read/write status depends on the step (see -s).'
    print '                   This file is only read when [-s 3], and is created with [-s 1] or [-s 2].'
    print '                                              '
    print ' [-f; --force]  >  (OPTIONAL) Flag to force overwrite. Do not promt user before overwriting existing files.'
    print '                                              '
    print ' [-m; --mute]   >  (OPTIONAL) Flag to run in silent mode; minimum logs.'
    print '                                              '
    print ' [-h; --help]   >  Show help, exit.'
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


def commandline_support(log=False):
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
            opts, args = getopt.getopt(sys.argv[1:], "s:o:t:fmhc",
                ['dict1=', 'dict2=', 'dict3=', 'dict4=', 'force', 'help', 'mute', 'code'])
        except getopt.GetoptError as err:
            print str(err)
            usage()
            sys.exit(2)
        if log: print 'detected options:  ', opts
        if log: print 'detected arguments:', args
        for opt, arg in opts:
            if opt in ['-m', '--mute']:
                P['log'] = False
            elif opt in ['-s']:
                P['step'] = arg
            elif opt in ['-t']:
                P['nc_in'] = [os.path.abspath(str(arg))]
                if len(args) > 0:
                    P['nc_in'] = [os.path.abspath(str(arg))]+[os.path.abspath(nc_name) for nc_name in args]
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
            elif opt in ['-o']:
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



def promtYesNo(question='', quitonno=False):
    '''
    function promts a yes-no question in console and catches answer,
    based on raw_input()
    ! WARNING! raw_input() may not work properly in SublimeText texteditor

    Inputs:
        question    -  string, will show this text before [y/n].
        quitonno    -  True/False. If True will force exit application if "No" answered

    Output:
        True        - if answer is Yes
        False       - if answer is No
        sys.exit(2) - if answer is No and quitonno=True
    '''
    yes_answers = ['y', 'Y', 'yes', 'Yes', 'YES']
    no_answers  = ['n', 'N', 'no', 'No', 'NO']
    a = None  # answer

    question += ' [y/n]'
    
    while a not in yes_answers+no_answers:
        a = raw_input(question)
        if a in yes_answers:
            return True
        elif a in no_answers:
            if quitonno:
                print 'Aborting... Good bye'
                sys.exit(2)
            return False
