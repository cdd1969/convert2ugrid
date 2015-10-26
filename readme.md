# convert2ugrid
## Version
0.3

## Description
This script is able to process data stored in netcdf files. It can convert structured to unstructured grid and map data into the new one. Originally it has been developed to visualize output results of [MOSSCO](http://www.mossco.de/) software with [DAVIT](http://www.smileconsult.de/index.php?article_id=26&clang=1) tool. In addition files generated by *convert2ugrid* were succesfully opened with [delft3D QUICKPLOT](http://oss.deltares.nl/web/delft3d) and [NCPLOT](http://wiki.baw.de/methoden/index.php5/NCPLOT).

## Requirements

* OS
    * Linux SUSE/Ubuntu (has not been tested on other systems yet)
* Python related:
    * [python](https://www.python.org/) >= v2.7.6
    * [netcdf4-python](https://github.com/Unidata/netcdf4-python) >= v1.0.8
    * [numpy](http://www.numpy.org/) >= v1.8.0
* General:
    * [HDF5](https://www.hdfgroup.org/HDF5/) >= v1.8.11
    * [NetCDF C](https://github.com/Unidata/netcdf-c) >= v4.3.0


## Basic usage
It is recommended to use console interface in Linux (since script requeres user-interaction through ```raw_input()``` function, which is known not to work properly with some UI (i.e. *SublimeText*)):
```sh
$ python convert2ugrid.py -t topo.nc model_results.nc
```
This will produce three files in your current working directory:
* dictionary_2
* dictionary_4
* ugrid.nc

To get description of possible command-line options:
```sh
$ python convert2ugrid.py -h
```    

## Examples (for BAW internal)
Before you continue set up the relevant path for the ***convert2ugrid*** script:
```sh
$ alias convert2ugrid='python /net/themis/system/akprog/python/qad/convert2ugrid/convert2ugrid.py'
```
### example 1 - *NSBS*
This example covers North and Baltic sea with a coarse grid. File */topo.nc* contains bathymetry defined at T points and lats/lons of grid cell centers. File */netcdf_reference_3d.nc* contains *MOSSCO* simulation output. Vertical layers are defined through sigma coordinates.

key info:
- bathymetry at T points
- lats / lons at T points
- rectangular grid
- vertical sigma layers

Now let's convert it! First define paths to data-files:
```sh
$ export TOPO=/net/widar/home/ak2stud/Nick/python_scripts/convert2ugrid/examples/1_nsbs/data/topo.nc
$ export NC1=/net/widar/home/ak2stud/Nick/python_scripts/convert2ugrid/examples/1_nsbs/data/netcdf_reference_3d.nc
```
Jump to desired directory and run the script:
```sh
$ cd <desired directory>
$ convert2ugrid -t $TOPO $NC1
```

### example 2 - *SNS*
This example covers south part of the North Sea with much finer grid. File */topo.nc* contains bathymetry defined at T points and lats/lons as well as y/x of grid nodes.

key info:
* bathymetry at T points
* lats / lons at X points
* curvilinear grid

First define paths to data-files:
```sh
$ export TOPO=/net/widar/home/ak2stud/Nick/python_scripts/convert2ugrid/examples/2_sns/data/topo_extended.nc
```
Jump to desired directory and run the script:
```sh
$ cd <desired directory>
$ convert2ugrid -t $TOPO
```
**NOTE:** now the file *ugrid.nc* has been produced, and is readable with *QUICKPLOT*, but not with *DAVIT* (due to some limitations). In order to overcome them we will include time-dependent dummy variable.

To do that, we will use predefined *dictionary_2* (you may compare it to the *dictionary_2* you have just generated), and start the script from step 2 (see documentation).
```sh
$ export DICT=/net/widar/home/ak2stud/Nick/python_scripts/convert2ugrid/examples/2_sns/data/dictionary_2
```
```sh
$ convert2ugrid -s 2 --dict2=$DICT -t $TOPO
```

Now *ugrid.nc* has additionaly one dummy-variable and can be viewed in *DAVIT*. 
#
---
## TODO LIST

- [+]: include numpy.squeeze in fetching ncdf data
- [+]: add support of bathymetry at X\_points currently only bathymetry located at T\_points is accepted
- [ ]: test bathymetry at X_points

- [ ]: update documentation from v0.2 to v0.3

- [ ]: add more flexible support to choose **depth\_variable**, since it my vary and the name is not only limited to **GETM\_grid\_3D** or **levels**


- [ ]: check more usage of new classes ```cdlVariable()``` and ```cdlVariableExt()```