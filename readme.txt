--------------------------
|     convert2ugrid      |
|         v 0.2          |
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


C) To use the convertor do the following steps:
        1) open file ./convert2ugrid.py
        2) scroll down to <if __name__ == '__main__':> part
        3) specify inputs (see documentation)
        4) specify parameters in function <create_davit_friendly_netcdf()>
        5) run the script

----------------------------------------------------------------------------------
TODO: add support of bathymetry at X_points
      currently bathymetry should be located at T_points