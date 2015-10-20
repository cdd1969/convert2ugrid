import sys
import os
import getopt


def usage():
    print '\n\n\nconvert2ugrid.py\nWill produce 3 files "dictionary_2", "dictionary_4" and user defined "nc_out"'
    print 'Both dictionaries will be created at the same directory where "nc_out" is specified'
    print '\tExample usage:'
    print ' $python convert2ugrid.py [-s n] [-i s1,s2,...] [--dict1=s] [--dict3=s] [-o s] [-f] [-m] [-h]'
    print '*'*100
    print ' [-s n; --step=n]          >  start from step <n>, integer (refer to documentation)'
    print '                              n=1 => create dict2, dict4 and NC'
    print '                              n=2 => create dict4 and NC. Dict2 should be already created'
    print '                              n=3 => create only NC. dict2 and dict4 should be already created'
    print ' [-i s1,s2; --input=s1,s2] >  <s1,s2,...> are the name(-s) of the input file(-s), comma-separated string'
    print '                              First file should always be topo-file (!!!)'
    print '                              example "read 1 file" : <-i topo.nc>'
    print '                              example "read more then one file": <-i topo.nc,data1.nc,data2.nc>'
    print ' [--dict1=s]               >  <s> is the name of the input dictionary1, string'
    print '                              If not set, will use default under ./lib/defaults/dictionary_1'
    print ' [--dict3=s]               >  <s> is the name of the input dictionary3, string'
    print '                              If not set, will use default under ./lib/defaults/dictionary_3'
    print ' [-o s; --output=s]        >  <s> is the name of the outputfile to be created, string'
    print ' [-f; --force]             >  flag to force overwrite. Do not promt user before overwriting existing files'
    print ' [-m; --mute]              >  flag to run in silent mode, minimum logs'
    print ' [-h; --help]              >  show help, exit'


def commandline_support():
    P = dict()
    # defaults....
    P['log'] = True
    P['step'] = 1
    P['overwrite'] = False
    # dont do anything if running without params
    if len(sys.argv) == 1:
        #use defaults
        pass
    else:
        # if more than 1 arg has been passed...  here the fun starts
        try:
            opts, args = getopt.getopt(sys.argv[1:], "s:o:i:fmh",
                ['step=', 'input=', 'dict1=', 'dict3=' 'output=', 'force', "help", "mute", ])
        except getopt.GetoptError as err:
            print str(err)
            usage()
            sys.exit(2)

        for opt, arg in opts:
            if opt in ['-m', '--mute']:
                P['log'] = False
            elif opt in ['-s', '--step']:
                P['step'] = int(arg)
            elif opt in ['-i', '--input']:
                P['nc_in'] = str(arg).split(',')
            elif opt in ['--dict1']:
                P['dict1'] = str(arg)
            elif opt in ['--dict3']:
                P['dict3'] = str(arg)
            elif opt in ['-o', '--output']:
                P['nc_out'] = str(arg)
            elif opt in ['-f', '--force']:
                P['overwrite'] = True
            
            elif opt in ['-h', '--help']:
                usage()
                sys.exit(2)

    return P

