try:
    import click
    CLICK_IS_HERE = True
    import ui2 as ui
except:
    CLICK_IS_HERE = False
    import ui
