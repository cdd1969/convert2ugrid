--------------------------
|     convert2ugrid      |
|         v 0.3          |
|    mossco >>> davit    |
--------------------------

The convertor-script has been tested with following configuration:
    Python related:
    -	Python (2.7.6)
    -	netcdf4-python (1.0.8)
    -	numpy (1.8.0)
    General:
    -	HDF5 (1.8.11)
    -	NetCDF C (4.3.0)



A) Documentation can be found under:
        ./info/convert2ugrid_documentation.pdf


B) A tutorial is described in HOWTO section of documentation


C) Basic usage:
   in console type:
        $ python convert2ugrid -t topo
   for help:
        $ python convert2ugrid -h


D) This tool is developed to be used from console, since it requeres user-interaction
   and it is known that some UI (i.e. SublimeText) does not work properly with
   raw_input() command.
    

----------------------------------------------------------------------------------
TODO: update documentation from v0.2 to v0.3

TODO: add support of bathymetry at X_points
      currently only bathymetry located at T_points is accepted

TODO: add more flexible support to choose <depth_variable>, since it my vary and the name
      is not only limited to <GETM_grid_3D> or <levels>

TODO: include numpy.squeeze in fetching ncdf data

TODO: check more usage of new classes cdlVariable() and cdlVariableExt()