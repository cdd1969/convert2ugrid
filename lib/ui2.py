import click

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
        click.pause(args[0])
        return
    if 'color' in kwargs:
        color = kwargs.pop('color')
        return click.prompt(click.style(args[0], fg=color), args[1::], **kwargs)
    else:
        return click.prompt(*args, **kwargs)

