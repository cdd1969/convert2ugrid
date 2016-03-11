from __future__ import print_function

try:
    import click
    CLICK_IS_HERE = True
    import ui2 as ui
except:
    CLICK_IS_HERE = False
    import ui


def sprint(*args, **kwargs):
    ''' Smart print with:
        * indenting
        * 2 print or not 2 print...
        * mode (warning, fail, False)
    Args:
        *args (str):
            strings to print
        **kwargs:
            log
            mode
            indent
    '''
    log = kwargs.get('log', True)
    mode = kwargs.get('mode', False)
    indent = kwargs.get('indent', False)
    
    if not log:
        return


    message = ''.join([str(arg) for arg in args])


    if mode in ['Warning', 'WARNING', 'warning']:
        message = '\033[93m' + message + '\033[0m'
    elif mode in ['FAIL', 'Fail', 'fail']:
        message = '\033[91m' + message + '\033[0m'
    elif mode in ['OK', 'Ok', 'ok']:
        message = '\033[92m' + message + '\033[0m'
    elif mode in ['BOLD', 'Bold', 'bold']:
        message = '\033[1m' + message + '\033[0m'
    else:
        pass

    if indent is not False:
        message = indent + message

    print (message)
