#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# This file is part of "convert2ugrid" tool
#
# Author: Nikolai Chernikov, nikolai.chernikov.ru@gmail.com


'''
This module contains functions relevant for processing NETCDF file, which is being createDimension
for further evaluation with DAVIT
'''
from __future__ import division
from netCDF4 import Dataset
import numpy as np
import time
import process_mossco_netcdf
import process_mixed_data


def create_uGrid_ncdf(filename,
                        nMesh2_node, nMesh2_edge, nMesh2_face, nMaxMesh2_face_nodes,
                        Mesh2_edge_nodes, Mesh2_edge_faces, Mesh2_face_nodes, Mesh2_face_edges,
                        Mesh2_node_x, Mesh2_node_y, Mesh2_edge_x, Mesh2_edge_y,
                        Mesh2_face_x, Mesh2_face_y, Mesh2_face_center_x, Mesh2_face_center_y,
                        Mesh2_face_area=None,
                        Mesh2_edge_x_bnd=None, Mesh2_edge_y_bnd=None,
                        Mesh2_face_x_bnd=None, Mesh2_face_y_bnd=None,
                        coord_type='geographic',
                        dim_nMesh2_layer2d=1, dim_nMesh2_layer3d=1, dim_nMesh2_class_names_strlen=20, dim_nMesh2_suspension_classes=1,
                        log=True):
    '''
    Function creates a NETCDF4 file. Data is stored in accordance with
    BAW convention for 2D Unstructured Grid (http://www.baw.de/methoden/index.php5/NetCDF_Unstrukturiertes_Gitter)
    
    input:
        filename - string, containing filename of netcdf file to be created.
        nMesh2_node, nMesh2_edge, nMesh2_face - integers, indicating number of nodes/edges/faces in a grid
        nMaxMesh2_face_nodes - integer, showing maximum number of nodes/edges in a face (could be 3 or 4)

        Mesh2_edge_nodes, Mesh2_edge_faces, Mesh2_face_nodes, Mesh2_face_edges,     |
        Mesh2_node_x, Mesh2_node_y, Mesh2_edge_x, Mesh2_edge_y,                     | => 1D numpy arrays with data
        Mesh2_face_x, Mesh2_face_y, Mesh2_face_center_x, Mesh2_face_center_y,       |
        Mesh2_face_area                                                             |

        coord_type  - string indicating in which variable to store passed data, in x,y or in lon,lat
                by default - 'geographic'
                'cartesian' <> 'geographic'
    '''


    # --------------------------------------------------
    #                   Creating ncdf
    # --------------------------------------------------
    root_grp = Dataset(filename, mode='w', format='NETCDF4')

    root_grp.title = 'mossco >>> uGrid conversion'
    root_grp.history = 'Createded on ' + time.ctime(time.time())
    #root_grp.comment = 'WARNING! Latitude and longitude valueas are stored as X,Y coordinates. Therefore any calculations that involve length or area in X or Y dimension are wrong'
    root_grp.Conventions = 'CF-1.6'
    #root_grp.institution = 'Bundesanstalt fuer Wasserbau - Federal Waterways Engineering and Research Institute'
    root_grp.references = 'http://www.baw.de/ und http://www.baw.de/methoden/index.php5/NetCDF'
    #root_grp.source = ''



    # --------------------------------------------------
    #                   Creating dimensions
    # --------------------------------------------------
    root_grp.createDimension('nMesh2_node', nMesh2_node)
    root_grp.createDimension('nMesh2_edge', nMesh2_edge)
    root_grp.createDimension('nMesh2_face', nMesh2_face)
    root_grp.createDimension('nMaxMesh2_face_nodes', nMaxMesh2_face_nodes)
    root_grp.createDimension('two', 2)
    root_grp.createDimension('three', 3)
    root_grp.createDimension('nMesh2_time', 1)
    root_grp.createDimension('nMesh2_data_time', None)  # None stands for UNLIMITED
    if dim_nMesh2_layer2d is not None:
        root_grp.createDimension('nMesh2_layer_2d', dim_nMesh2_layer2d)
    if dim_nMesh2_layer3d is not None:
        root_grp.createDimension('nMesh2_layer_3d', dim_nMesh2_layer3d)
    root_grp.createDimension('nMesh2_class_names_strlen', dim_nMesh2_class_names_strlen)
    root_grp.createDimension('nMesh2_suspension_classes', dim_nMesh2_suspension_classes)

    if coord_type in ['cartesian', 'both'] :
        # **********************************************************************************************************************************************
        #
        #                  1) Local coordinates
        #
        # **********************************************************************************************************************************************
        # --------------------------------------------------
        #                                   1.1) NODES
        # --------------------------------------------------

        ncVar_Mesh2_node_x = root_grp.createVariable('Mesh2_node_x', 'f8', ('nMesh2_node'), fill_value=False)
        ncVar_Mesh2_node_x.long_name = 'x-Koordinate der Knoten eines 2D-Gitters'
        ncVar_Mesh2_node_x.units = 'm'
        ncVar_Mesh2_node_x.name_id = 1650
        ncVar_Mesh2_node_x.standard_name = 'projection_x_coordinate'
        ncVar_Mesh2_node_x[:] = Mesh2_node_x[:]

        ncVar_Mesh2_node_y = root_grp.createVariable('Mesh2_node_y', 'f8', ('nMesh2_node'), fill_value=False)
        ncVar_Mesh2_node_y.long_name = 'y-Koordinate der Knoten eines 2D-Gitters'
        ncVar_Mesh2_node_y.units = 'm'
        ncVar_Mesh2_node_y.name_id = 1651
        ncVar_Mesh2_node_y.standard_name = 'projection_y_coordinate'
        ncVar_Mesh2_node_y[:] = Mesh2_node_y[:]

        # --------------------------------------------------
        #                                   1.2) EDGES
        # --------------------------------------------------

        ncVar_Mesh2_edge_x = root_grp.createVariable('Mesh2_edge_x', 'f8', ('nMesh2_edge'), fill_value=False)
        ncVar_Mesh2_edge_x.long_name = 'x-Koordinate der Kanten eines 2D-Gitters, Kantenmitte'
        ncVar_Mesh2_edge_x.units = 'm'
        ncVar_Mesh2_edge_x.name_id = 1650
        ncVar_Mesh2_edge_x.bounds = 'Mesh2_edge_x_bnd'
        ncVar_Mesh2_edge_x.standard_name = 'projection_x_coordinate'
        ncVar_Mesh2_edge_x[:] = Mesh2_edge_x[:]

        ncVar_Mesh2_edge_y = root_grp.createVariable('Mesh2_edge_y', 'f8', ('nMesh2_edge'), fill_value=False)
        ncVar_Mesh2_edge_y.long_name = 'y-Koordinate der Kanten eines 2D-Gitters, Kantenmitte'
        ncVar_Mesh2_edge_y.units = 'm'
        ncVar_Mesh2_edge_y.name_id = 1651
        ncVar_Mesh2_edge_y.bounds = 'Mesh2_edge_y_bnd'
        ncVar_Mesh2_edge_y.standard_name = 'projection_y_coordinate'
        ncVar_Mesh2_edge_y[:] = Mesh2_edge_y[:]

        # --------------------------------------------------
        #                                   1.3) FACES
        # --------------------------------------------------

        ncVar_Mesh2_face_x = root_grp.createVariable('Mesh2_face_x', 'f8', ('nMesh2_face'), fill_value=False)
        ncVar_Mesh2_face_x.long_name = 'x-Koordinate der Faces (Polygone) eines 2D-Gitters, Schwerpunkt'
        ncVar_Mesh2_face_x.units = 'm'
        ncVar_Mesh2_face_x.name_id = 1650
        ncVar_Mesh2_face_x.bounds = 'Mesh2_face_x_bnd'
        ncVar_Mesh2_face_x.standard_name = 'projection_x_coordinate'
        ncVar_Mesh2_face_x[:] = Mesh2_face_x[:]

        ncVar_Mesh2_face_y = root_grp.createVariable('Mesh2_face_y', 'f8', ('nMesh2_face'), fill_value=False)
        ncVar_Mesh2_face_y.long_name = 'y-Koordinate der Faces (Polygone) eines 2D-Gitters, Schwerpunkt'
        ncVar_Mesh2_face_y.units = 'm'
        ncVar_Mesh2_face_y.name_id = 1651
        ncVar_Mesh2_face_y.bounds = 'Mesh2_face_y_bnd'
        ncVar_Mesh2_face_y.standard_name = 'projection_y_coordinate'
        ncVar_Mesh2_face_y[:] = Mesh2_face_y[:]

        ncVar_Mesh2_face_center_x = root_grp.createVariable('Mesh2_face_center_x', 'f8', ('nMesh2_face'), fill_value=False)
        ncVar_Mesh2_face_center_x.long_name = 'x-Koordinate der Faces (Polygone) eines 2D-Gitters, Umkreismittelpunkt'
        ncVar_Mesh2_face_center_x.units = 'm'
        ncVar_Mesh2_face_center_x.name_id = 1650
        ncVar_Mesh2_face_center_x.bounds = 'Mesh2_face_x_bnd'
        ncVar_Mesh2_face_center_x.standard_name = 'projection_x_coordinate'
        ncVar_Mesh2_face_center_x[:] = Mesh2_face_center_x[:]

        ncVar_Mesh2_face_center_y = root_grp.createVariable('Mesh2_face_center_y', 'f8', ('nMesh2_face'), fill_value=False)
        ncVar_Mesh2_face_center_y.long_name = 'y-Koordinate der Faces (Polygone) eines 2D-Gitters, Umkreismittelpunkt'
        ncVar_Mesh2_face_center_y.units = 'm'
        ncVar_Mesh2_face_center_y.name_id = 1651
        ncVar_Mesh2_face_center_y.bounds = 'Mesh2_face_y_bnd'
        ncVar_Mesh2_face_center_y.standard_name = 'projection_y_coordinate'
        ncVar_Mesh2_face_center_y[:] = Mesh2_face_center_y[:]

        # --------------------------------------------------
        #                                   1.4) OPTIONAL / BORDERS
        # --------------------------------------------------
        if Mesh2_edge_x_bnd is not None and Mesh2_edge_y_bnd is not None:
            ncVar_Mesh2_edge_x_bnd = root_grp.createVariable('Mesh2_edge_x_bnd', 'f8', ('nMesh2_edge', 'two'), fill_value=False)
            ncVar_Mesh2_edge_y_bnd = root_grp.createVariable('Mesh2_edge_y_bnd', 'f8', ('nMesh2_edge', 'two'), fill_value=False)
            ncVar_Mesh2_edge_x_bnd[...] = Mesh2_edge_x_bnd
            ncVar_Mesh2_edge_y_bnd[...] = Mesh2_edge_y_bnd
        if Mesh2_face_x_bnd is not None and Mesh2_face_y_bnd is not None:
            ncVar_Mesh2_face_x_bnd = root_grp.createVariable('Mesh2_face_x_bnd', 'f8', ('nMesh2_face', 'nMaxMesh2_face_nodes'), fill_value=-999)
            ncVar_Mesh2_face_y_bnd = root_grp.createVariable('Mesh2_face_y_bnd', 'f8', ('nMesh2_face', 'nMaxMesh2_face_nodes'), fill_value=-999)
            ncVar_Mesh2_face_x_bnd[...] = Mesh2_face_x_bnd
            ncVar_Mesh2_face_y_bnd[...] = Mesh2_face_y_bnd

        if coord_type == 'cartesian':
            # ------------------------------------------------------------------------------------
            #                                   3.5) Topology
            # ------------------------------------------------------------------------------------

            ncVar_Mesh2 = root_grp.createVariable('Mesh2', 'i', fill_value=False)
            ncVar_Mesh2.long_name = 'UnTRIM Gitternetz, Drei- und Vierecke gemischt, kein SubGrid'
            ncVar_Mesh2.cf_role = 'mesh_topology'
            ncVar_Mesh2.dimensionality = 2
            ncVar_Mesh2.node_coordinates = 'Mesh2_node_x Mesh2_node_y Mesh2_node_lon Mesh2_node_lat'
            ncVar_Mesh2.edge_coordinates = 'Mesh2_edge_x Mesh2_edge_y Mesh2_edge_lon Mesh2_edge_lat'
            ncVar_Mesh2.face_coordinates = 'Mesh2_face_x Mesh2_face_y Mesh2_face_lon Mesh2_face_lat Mesh2_face_center_x Mesh2_face_center_y Mesh2_face_center_lon Mesh2_face_center_lat'
            ncVar_Mesh2.edge_node_connectivity = 'Mesh2_edge_nodes'
            ncVar_Mesh2.edge_face_connectivity = 'Mesh2_edge_faces'
            ncVar_Mesh2.face_node_connectivity = 'Mesh2_face_nodes'
            ncVar_Mesh2.face_edge_connectivity = 'Mesh2_face_edges'

            # ------------------------------------------------------------------------------------
            #                                   3.6) Area
            # ------------------------------------------------------------------------------------
            if Mesh2_face_area is not None:
                ncVar_Mesh2_face_area = root_grp.createVariable('Mesh2_face_area', 'f8', ('nMesh2_face'), fill_value=1.e+31)
                ncVar_Mesh2_face_area.long_name = 'Zellenflaeche'
                ncVar_Mesh2_face_area.units = 'm2'
                ncVar_Mesh2_face_area.name_id = 1656
                ncVar_Mesh2_face_area.valid_range = [Mesh2_face_area[:].min(), Mesh2_face_area[:].max()]
                ncVar_Mesh2_face_area.coordinates = 'Mesh2_face_x Mesh2_face_y Mesh2_face_lon Mesh2_face_lat'
                ncVar_Mesh2_face_area.grid_mapping = 'Mesh2_crs'
                ncVar_Mesh2_face_area.standard_name = 'cell_area'
                ncVar_Mesh2_face_area.mesh = 'Mesh2'
                ncVar_Mesh2_face_area.location = 'face'
                ncVar_Mesh2_face_area [:] = Mesh2_face_area[:]
    
    elif coord_type in ['geographic', 'both']:
        # *********************************************************************************************************************************************
        #
        #                  2) Georaphical coordinates
        #
        # *********************************************************************************************************************************************
        ncVar_Mesh2_crs = root_grp.createVariable('Mesh2_crs', 'int', (), fill_value=False)
        ncVar_Mesh2_crs[:] = 4326
        ncVar_Mesh2_crs.epsg_code = "EPSG:4326"
        ncVar_Mesh2_crs.comment = "LON, LAT : WGS84, EPSG:4326"
        ncVar_Mesh2_crs.grid_mapping_name = "latitude_longitude"
        #ncVar_Mesh2_crs.grid_mapping_name = "rotated_latitude_longitude"
        ncVar_Mesh2_crs.longitude_of_prime_meridian = 0.0
        ncVar_Mesh2_crs.semi_major_axis = 6378137.0
        ncVar_Mesh2_crs.inverse_flattening = 298.257223563
        # --------------------------------------------------
        #                                   2.1) NODES
        # --------------------------------------------------

        ncVar_Mesh2_node_lon = root_grp.createVariable('Mesh2_node_lon', 'f8', ('nMesh2_node'), fill_value=False)
        ncVar_Mesh2_node_lon.long_name = 'geografische Laenge der 2D-Gitter-Knoten'
        ncVar_Mesh2_node_lon.units = 'degrees_east'
        ncVar_Mesh2_node_lon.name_id = 1653
        ncVar_Mesh2_node_lon.standard_name = 'longitude'
        ncVar_Mesh2_node_lon[:] = Mesh2_node_x[:]

        ncVar_Mesh2_node_lat = root_grp.createVariable('Mesh2_node_lat', 'f8', ('nMesh2_node'), fill_value=False)
        ncVar_Mesh2_node_lat.long_name = 'geografische Breite der 2D-Gitter-Knoten'
        ncVar_Mesh2_node_lat.units = 'degrees_north'
        ncVar_Mesh2_node_lat.name_id = 1652
        ncVar_Mesh2_node_lat.standard_name = 'latitude'
        ncVar_Mesh2_node_lat[:] = Mesh2_node_y[:]

        # --------------------------------------------------
        #                                   2.2) EDGES
        # --------------------------------------------------

        ncVar_Mesh2_edge_lon = root_grp.createVariable('Mesh2_edge_lon', 'f8', ('nMesh2_edge'), fill_value=False)
        ncVar_Mesh2_edge_lon.long_name = 'geografische Laenge der 2D-Gitter-Kanten, Kantenmitte'
        ncVar_Mesh2_edge_lon.units = 'degrees_east'
        ncVar_Mesh2_edge_lon.name_id = 1653
        ncVar_Mesh2_edge_lon.bounds = 'Mesh2_edge_lon_bnd'
        ncVar_Mesh2_edge_lon.standard_name = 'longitude'
        ncVar_Mesh2_edge_lon[:] = Mesh2_edge_x[:]

        ncVar_Mesh2_edge_lat = root_grp.createVariable('Mesh2_edge_lat', 'f8', ('nMesh2_edge'), fill_value=False)
        ncVar_Mesh2_edge_lat.long_name = 'geografische Breite der 2D-Gitter-Kanten, Kantenmitte'
        ncVar_Mesh2_edge_lat.units = 'degrees_north'
        ncVar_Mesh2_edge_lat.name_id = 1652
        ncVar_Mesh2_edge_lat.bounds = 'Mesh2_edge_lat_bnd'
        ncVar_Mesh2_edge_lat.standard_name = 'latitude'
        ncVar_Mesh2_edge_lat[:] = Mesh2_edge_y[:]

        # --------------------------------------------------
        #                                   2.3) FACES
        # --------------------------------------------------

        ncVar_Mesh2_face_lon = root_grp.createVariable('Mesh2_face_lon', 'f8', ('nMesh2_face'), fill_value=False)
        ncVar_Mesh2_face_lon.long_name = 'geografische Laenge der 2D-Gitter-Faces (Polygone), Schwerpunkt'
        ncVar_Mesh2_face_lon.units = 'degrees_east'
        ncVar_Mesh2_face_lon.name_id = 1653
        ncVar_Mesh2_face_lon.bounds = 'Mesh2_face_lon_bnd'
        ncVar_Mesh2_face_lon.standard_name = 'longitude'
        ncVar_Mesh2_face_lon[:] = Mesh2_face_x[:]

        ncVar_Mesh2_face_lat = root_grp.createVariable('Mesh2_face_lat', 'f8', ('nMesh2_face'), fill_value=False)
        ncVar_Mesh2_face_lat.long_name = 'geografische Breite der 2D-Gitter-Faces (Polygone), Schwerpunkt'
        ncVar_Mesh2_face_lat.units = 'degrees_north'
        ncVar_Mesh2_face_lat.name_id = 1652
        ncVar_Mesh2_face_lat.bounds = 'Mesh2_face_lat_bnd'
        ncVar_Mesh2_face_lat.standard_name = 'latitude'
        ncVar_Mesh2_face_lat[:] = Mesh2_face_y[:]

        ncVar_Mesh2_face_center_lon = root_grp.createVariable('Mesh2_face_center_lon', 'f8', ('nMesh2_face'), fill_value=False)
        ncVar_Mesh2_face_center_lon.long_name = 'geografische Laenge der 2D-Gitter-Faces (Polygone), Umkreismittelpunkt'
        ncVar_Mesh2_face_center_lon.units = 'degrees_east'
        ncVar_Mesh2_face_center_lon.name_id = 1653
        ncVar_Mesh2_face_center_lon.standard_name = 'longitude'
        ncVar_Mesh2_face_center_lon[:] = Mesh2_face_center_x[:]

        ncVar_Mesh2_face_center_lat = root_grp.createVariable('Mesh2_face_center_lat', 'f8', ('nMesh2_face'), fill_value=False)
        ncVar_Mesh2_face_center_lat.long_name = 'geografische Breite der 2D-Gitter-faces (Polygone), Umkreismittelpunkt'
        ncVar_Mesh2_face_center_lat.units = 'degrees_north'
        ncVar_Mesh2_face_center_lat.name_id = 1652
        ncVar_Mesh2_face_center_lat.standard_name = 'latitude'
        ncVar_Mesh2_face_center_lat[:] = Mesh2_face_center_y[:]

        # --------------------------------------------------
        #                                   2.4) OPTIONAL / BORDERS
        # --------------------------------------------------
        if Mesh2_edge_x_bnd is not None and Mesh2_edge_y_bnd is not None:
            ncVar_Mesh2_edge_lon_bnd = root_grp.createVariable('Mesh2_edge_lon_bnd', 'f8', ('nMesh2_edge', 'two'), fill_value=False)
            ncVar_Mesh2_edge_lat_bnd = root_grp.createVariable('Mesh2_edge_lat_bnd', 'f8', ('nMesh2_edge', 'two'), fill_value=False)
            ncVar_Mesh2_edge_lon_bnd[...] = Mesh2_edge_x_bnd
            ncVar_Mesh2_edge_lat_bnd[...] = Mesh2_edge_y_bnd
        if Mesh2_face_x_bnd is not None and Mesh2_face_y_bnd is not None:
            ncVar_Mesh2_face_lon_bnd = root_grp.createVariable('Mesh2_face_lon_bnd', 'f8', ('nMesh2_face', 'nMaxMesh2_face_nodes'), fill_value=-999)
            ncVar_Mesh2_face_lat_bnd = root_grp.createVariable('Mesh2_face_lat_bnd', 'f8', ('nMesh2_face', 'nMaxMesh2_face_nodes'), fill_value=-999)
            ncVar_Mesh2_face_lon_bnd[...] = Mesh2_face_x_bnd
            ncVar_Mesh2_face_lat_bnd[...] = Mesh2_face_y_bnd

        # ------------------------------------------------------------------------------------
        #                                   3.5) Topology
        # ------------------------------------------------------------------------------------

        ncVar_Mesh2 = root_grp.createVariable('Mesh2', 'i', fill_value=False)
        #ncVar_Mesh2.long_name = 'UnTRIM Gitternetz, Drei- und Vierecke gemischt, kein SubGrid'
        ncVar_Mesh2.long_name = 'UnTRIM Gitternetz, Vierecke, kein SubGrid'
        ncVar_Mesh2.cf_role = 'mesh_topology'
        ncVar_Mesh2.dimensionality = 2
        ncVar_Mesh2.node_coordinates = 'Mesh2_node_lon Mesh2_node_lat Mesh2_node_x Mesh2_node_y'
        ncVar_Mesh2.edge_coordinates = 'Mesh2_edge_lon Mesh2_edge_lat Mesh2_edge_x Mesh2_edge_y'
        ncVar_Mesh2.face_coordinates = 'Mesh2_face_lon Mesh2_face_lat Mesh2_face_x Mesh2_face_y Mesh2_face_center_lon Mesh2_face_center_lat Mesh2_face_center_x Mesh2_face_center_y '
        ncVar_Mesh2.edge_node_connectivity = 'Mesh2_edge_nodes'
        ncVar_Mesh2.edge_face_connectivity = 'Mesh2_edge_faces'
        ncVar_Mesh2.face_node_connectivity = 'Mesh2_face_nodes'
        ncVar_Mesh2.face_edge_connectivity = 'Mesh2_face_edges'


    
    else:
        err_msg = 'passed <coord_type = {0}> is invalid. Choose "cartesian", "geographic" or "both"'.format(coord_type)
        raise TypeError(err_msg)
    # **************************************************
    #
    #                  3) Topological coordinates
    #
    # **************************************************
    # --------------------------------------------------
    #                                   3.1) EDGE >>> NODES
    # --------------------------------------------------

    ncVar_Mesh2_edge_nodes = root_grp.createVariable('Mesh2_edge_nodes', 'i', ('nMesh2_edge', 'two'))
    ncVar_Mesh2_edge_nodes.long_name = 'Knotenverzeichnis der Kanten, Anfangs- und Endpunkt'
    ncVar_Mesh2_edge_nodes.cf_role = 'edge_node_connectivity'
    ncVar_Mesh2_edge_nodes.start_index = 0
    ncVar_Mesh2_edge_nodes[:] = Mesh2_edge_nodes[:]

    # --------------------------------------------------
    #                                   3.2) EDGE >>> FACES
    # --------------------------------------------------

    ncVar_Mesh2_edge_faces = root_grp.createVariable('Mesh2_edge_faces', 'i', ('nMesh2_edge', 'two'), fill_value=-999)
    ncVar_Mesh2_edge_faces.long_name = 'Face- (Polygon-) Verzeichnis der Kanten, linker und rechter Nachbar'
    ncVar_Mesh2_edge_faces.cf_role = 'edge_face_connectivity'
    ncVar_Mesh2_edge_faces.start_index = 0
    ncVar_Mesh2_edge_faces[:] = Mesh2_edge_faces[:]

    # --------------------------------------------------
    #                                   3.3) FACE >>> NODES
    # --------------------------------------------------

    ncVar_Mesh2_face_nodes = root_grp.createVariable('Mesh2_face_nodes', 'i', ('nMesh2_face', 'nMaxMesh2_face_nodes'), fill_value=-999)
    ncVar_Mesh2_face_nodes.long_name = 'Knotenverzeichnis der Faces (Polygone),entgegen dem Uhrzeigersinn'
    ncVar_Mesh2_face_nodes.cf_role = 'face_node_connectivity'
    ncVar_Mesh2_face_nodes.start_index = 0
    ncVar_Mesh2_face_nodes[:] = Mesh2_face_nodes[:]

    # ------------------------------------------------------------------------------------
    #                                   3.4) FACE >>> EDGES
    # ------------------------------------------------------------------------------------

    ncVar_Mesh2_face_edges = root_grp.createVariable('Mesh2_face_edges', 'i', ('nMesh2_face', 'nMaxMesh2_face_nodes'), fill_value=-999)
    ncVar_Mesh2_face_edges.long_name = 'Kantenverzeichnis der Faces (Polygone),entgegen dem Uhrzeigersinn'
    ncVar_Mesh2_face_edges.cf_role = 'face_edge_connectivity'
    ncVar_Mesh2_face_edges.start_index = 0
    ncVar_Mesh2_face_edges[:] = Mesh2_face_edges[:]


    # --------------------------------------------------
    #                   Closing ncdf
    # --------------------------------------------------
    root_grp.close()

    if log: print 'File created succesfully: %s' % (filename)
    




























def append_test_Mesh2_face_z_3d_and_Mesh2_face_z_3d_bnd(fname_davit, nx=0, ny=0, mask=None, log=False):
    '''
                    DUMMY VALUES!!!
    Function appends to DAVIT netcdf following variables:
        Mesh2_face_z_face_3d
        Mesh2_face_z_face_bnd_3d
        Mesh2_edge_z_edge_3d
        Mesh2_edge_z_edge_bnd_3d

    The data stored is artificial (!!!):
        - layers are in down-positiv order (first layer is at surface)
        - each layer is 1m thick
        - layer-borders are (0,1), (1,2), (2,3), ..., (NLAYERS-1, NLAYERS) in meters
    '''
    '''
    INPUT:
        fname_davit, fname_mossco  - strings, pathes to files
        mask - an boolean mask 2D array. See function process_mossco_netcdf.make_mask_array_from_mossco_bathymetry()
        log  - boolean Flag to print output in console


    CDL EXAMPLE:
    float Mesh2_face_z_face_3d(nMesh2_data_time, nMesh2_layer_3d, nMesh2_face) ;
        Mesh2_face_z_face_3d:long_name = "z_face [ face ]" ;
        Mesh2_face_z_face_3d:units = "m" ;
        Mesh2_face_z_face_3d:name_id = 1702 ;
        Mesh2_face_z_face_3d:positive = "down" ;
        Mesh2_face_z_face_3d:bounds = "Mesh2_face_z_face_bnd_3d" ;
        Mesh2_face_z_face_3d:standard_name = "depth" ;
    '''

    if log: print 'running: append_test_Mesh2_face_z_3d_and_Mesh2_face_z_3d_bnd()'
    
    # get dimensions....
    # ----------------------------------------------------------------------------   
    __nc = Dataset(fname_davit, mode='r')
    try:
        nt = __nc.dimensions['nMesh2_data_time'].__len__()
    except:
        nt = __nc.dimensions['nMesh2_time'].__len__()
    nz     = __nc.dimensions['nMesh2_layer_3d'].__len__()
    nFaces = __nc.dimensions['nMesh2_face'].__len__()
    nEdges = __nc.dimensions['nMesh2_edge'].__len__()
    __nc.close()

    if log:
        print 'nZ, nY, nX =', nz, ny, nx
        print 'nTimesteps =', nt
        print 'nFaces =', nFaces
        print 'nEdges =', nEdges


    # FACE middle values.....
    # ----------------------------------------------------------------------------
    variable = dict()
    variable['vname'] = 'Mesh2_face_z_face_3d'
    variable['dtype'] = 'f4'
    variable['dims'] = ('nMesh2_data_time', 'nMesh2_layer_3d', 'nMesh2_face')
    variable['_FillValue'] = None
    
    ATTRS = dict()
    ATTRS['long_name']     = 'z_face [ face ]'
    ATTRS['units']         = 'm'
    ATTRS['positive']      = 'down'
    ATTRS['name_id']       = 1702
    ATTRS['bounds']        = 'Mesh2_face_z_face_bnd_3d'
    ATTRS['standard_name'] = 'depth'
    ATTRS['comment']       = "warning: dummy values are used"
    variable['attributes'] = ATTRS

    davit_dummy = np.zeros((nt, nz, nFaces))
    for z in xrange(nz):
        davit_dummy[:, z, :] = nz-z-0.5  # this equation will produce values like 0.5, 1.5, 2.5, ..., nz-1.0-0.5
    
    variable['data'] = davit_dummy

    append_VariableData_to_netcdf(fname_davit, variable, log=log)
    del variable

    # FACE bounds values...
    # ----------------------------------------------------------------------------
    if log:
        print '+'*50
        print '+'*50
        print 'Bounds'
        print '+'*50
        print '+'*50
    variable = dict()
    variable['vname'] = 'Mesh2_face_z_face_bnd_3d'
    variable['dtype'] = 'f4'
    variable['dims'] = ('nMesh2_data_time', 'nMesh2_layer_3d', 'nMesh2_face', 'two')
    variable['_FillValue'] = None

    ATTRS = dict()
    ATTRS['name_id']       = 1703
    #ATTRS['units']         = 'm'
    ATTRS['comment']       = "warning: dummy values are used"
    variable['attributes'] = ATTRS

    davit_dummy = np.zeros((nt, nz, nFaces, 2))
    for z in xrange(nz):
        davit_dummy[:, z, :, 0] = nz-z-1.  # this equation will produce values like 0, 1, 2, ... nz-z
        davit_dummy[:, z, :, 1] = nz-z+0.  # this equation will produce values like 1, 2, 3, ... nz-z+1

    variable['data'] = davit_dummy

    append_VariableData_to_netcdf(fname_davit, variable, log=log)
    del variable


    if log:
        print '+'*50
        print '+'*50
        print 'EDGE'
        print '+'*50
        print '+'*50
    # EDGE middle values.....
    # ----------------------------------------------------------------------------
    variable = dict()
    variable['vname'] = 'Mesh2_edge_z_edge_3d'
    variable['dtype'] = 'f4'
    variable['dims'] = ('nMesh2_data_time', 'nMesh2_layer_3d', 'nMesh2_edge')
    variable['_FillValue'] = None
    
    ATTRS = dict()
    ATTRS['long_name']     = 'z_edge [ edge ]'
    ATTRS['units']         = 'm'
    ATTRS['positive']      = 'down'
    ATTRS['name_id']       = 1707
    ATTRS['bounds']        = 'Mesh2_edge_z_edge_bnd_3d'
    ATTRS['standard_name'] = 'depth'
    ATTRS['comment']       = "warning: dummy values are used"
    variable['attributes'] = ATTRS


    davit_dummy = np.zeros((nt, nz, nEdges))
    for z in xrange(nz):
        davit_dummy[:, z, :] = nz-z-0.5  # this equation will produce values like 0.5, 1.5, 2.5, ..., nz-1.0-0.5
    
    variable['data'] = davit_dummy

    append_VariableData_to_netcdf(fname_davit, variable, log=log)
    del variable

    # EDGE bounds values...
    # ----------------------------------------------------------------------------
    if log:
        print '+'*50
        print '+'*50
        print 'EDGE Bounds'
        print '+'*50
        print '+'*50
    variable = dict()
    variable['vname'] = 'Mesh2_edge_z_edge_bnd_3d'
    variable['dtype'] = 'f4'
    variable['dims'] = ('nMesh2_data_time', 'nMesh2_layer_3d', 'nMesh2_edge', 'two')
    variable['_FillValue'] = None

    ATTRS = dict()
    ATTRS['name_id']       = 1708
    #ATTRS['units']         = 'm'
    ATTRS['comment']       = "warning: dummy values are used"
    variable['attributes'] = ATTRS
    
    davit_dummy = np.zeros((nt, nz, nEdges, 2))

    for z in xrange(nz):
        davit_dummy[:, z, :, 0] = nz-z-1.  # this equation will produce values like 0, 1, 2, ... nz-z
        davit_dummy[:, z, :, 1] = nz-z+0.  # this equation will produce values like 1, 2, 3, ... nz-z+1

    variable['data'] = davit_dummy

    append_VariableData_to_netcdf(fname_davit, variable, log=log)
    if log: print 'finished: append_test_Mesh2_face_z_3d_and_Mesh2_face_z_3d_bnd()'



def append_Time_andDatetime_to_netcdf(fname_davitnc, fname_mossconc=None, time_var_name='time', dummy_values=False, log=True):
        # --------------------------------------------------
        # 3.1) fill netcdf with time variables (TIME and DATATIME)
        # --------------------------------------------------
        if dummy_values is False:

            __nc = Dataset(fname_mossconc, mode='r')
            time_var_units = __nc.variables[time_var_name].units
            max_t          = __nc.variables[time_var_name][:].max()
            min_t          = __nc.variables[time_var_name][:].min()
            __nc.close()
        else:
            time_var_units = 'seconds since 2000-01-01 00:00:00'
            max_t = 60*60*24
            min_t = 0



        var_t = dict()
        var_t['vname'] = 'nMesh2_time'
        var_t['dtype'] = 'f8'
        var_t['dims'] = ('nMesh2_time', )
        var_t['_FillValue'] = False

        ATTRS = dict()
        # because currently MOSSCO writes a time output as "hours since 2009-01-02 00:00:00",
        # it is nessesary to modify it. Will be removed, after corrections in MOSSCO.
        ATTRS['long_name'] = "time"
        units = time_var_units+' 01:00'  # '01:00 is indicating timezone, see CF conventions for details'
        ATTRS['units'] = units
        ATTRS['name_id'] = 1640
        ATTRS['axis'] = "T"
        ATTRS['bounds'] = "nMesh2_time_bnd"
        ATTRS['calendar'] = "gregorian"
        ATTRS['standard_name'] = "time"
        if dummy_values:
            ATTRS['comment'] = "warning: dummy values are used"
        var_t['attributes'] = ATTRS
        var_t['data'] = np.array([min_t])
        append_VariableData_to_netcdf(fname_davitnc, var_t, log=log)
        
        var_t_bnd = dict()
        var_t_bnd['vname'] = 'nMesh2_time_bnd'
        var_t_bnd['dtype'] = 'f8'
        var_t_bnd['dims'] = ('nMesh2_time', 'two')
        if dummy_values:
            ATTRS = dict()
            ATTRS['comment'] = "warning: dummy values are used"
            var_t['attributes'] = ATTRS
        var_t_bnd['_FillValue'] = False

        var_t_bnd['data'] = np.array([min_t, max_t])
        #var_t_bnd['data'] = np.array([0, 62985600])
        append_VariableData_to_netcdf(fname_davitnc, var_t_bnd, log=log)


        if dummy_values is False:
            var_dt = dict()
            var_dt['vname'] = 'nMesh2_data_time'
            var_dt['dtype'] = 'f8'
            var_dt['dims'] = ('nMesh2_data_time', )
            var_dt['_FillValue'] = False

            ATTRS = dict()
            ATTRS['long_name'] = "time"
            
            ATTRS['units'] = time_var_units + ' 01:00'
            ATTRS['name_id'] = 1640
            ATTRS['axis'] = "T"
            ATTRS['calendar'] = "gregorian"
            ATTRS['standard_name'] = "time"
            var_t['attributes'] = ATTRS
            
            var_dt['attributes'] = ATTRS
            var_dt['data'] = process_mossco_netcdf.read_mossco_nc_1d(fname_mossconc, 'time')
            append_VariableData_to_netcdf(fname_davitnc, var_dt, log=log)





def append_VariableData_to_netcdf(filename, variable, log=True):
    '''
    Function appends data to a NETCDF4 file created by function "create_uGrid_ncdf()". Check description there.
    Data is stored in accordance with BAW convention for 2D Unstructured Grid
    (http://www.baw.de/methoden/index.php5/NetCDF_Unstrukturiertes_Gitter)
    
    input:
        filename - string, containing filename of netcdf file to be created.

        variables - dictionary with data...
            variable['vname'] = string, 'Mesh2_node_depth'
            variable['dtype'] = string, 'f8' or 'f4' or 'i'
            variable['dims'] = tuple of strings, ('time', 'nMesh2_node')
            variable['_FillValue'] = False
            variable['attributes'] = dictionary with attributes

            variable['data'] = numpy array of dimensions specified in variable['dims'] or None.
    '''


    _n = 'append_VariableData_to_netcdf():'
    if log: print _n, 'appending variable: ', variable['vname']
    if log: print _n, 'input     datashape:', variable['data'].shape
    #if log: print _n, 'input    data-array:', variable['data']

    # --------------------------------------------------
    #                   Appending ncdf
    # --------------------------------------------------
    root_grp = Dataset(filename, mode='a')

    # check if this variable already exists!
    if variable['vname'] in root_grp.variables.keys():
        root_grp.close()
        if log:
            print _n, 'Variable <%s> skipped. Already exising' % (variable['vname'])
        return



    # now create Variable
    ncVar_data = root_grp.createVariable(variable['vname'], variable['dtype'], dimensions=variable['dims'], fill_value=variable['_FillValue'])
    
    # add attributes
    if 'attributes' in variable.keys():
        for attr_name, attr_value in variable['attributes'].iteritems():
            if not attr_name.startswith('_'):  # reserved for special attributes
                # remove one element lists...
                if isinstance(attr_value, list):
                    if len(attr_value) == 1:
                        attr_value = attr_value[0]

                # change long to int
                # since even if declared int() on 64bit python - it will become long()
                if isinstance(attr_value, np.int64) or isinstance(attr_value, long) or isinstance(attr_value, int):
                    attr_value = np.int32(attr_value)


                # treat real range values
                if attr_name in ['valid_range', 'range']:
                    attr_value = [np.float32(variable['data'].min()), np.float32(variable['data'].max())]
                elif attr_name in ['valid_max', 'max']:
                    attr_value = np.float32(variable['data'].max())
                elif attr_name in ['valid_min', 'min']:
                    attr_value = np.float32(variable['data'].min())

                if log: print _n, '\t adding attribute <{0} = {1}> of type {2}'.format(attr_name, attr_value, type(attr_value))

                # actually append
                ncVar_data.setncattr(attr_name, attr_value)
    
    # fill data
    if variable['data'] is not None:

        #print 'NC  before', ncVar_data.shape
        #print 'DAT before', variable['data'].shape

        # HARDCODE FOR UNLIMITED dim
        # will fail if first dim is nMesh2_layer_3d. Need to rework this !
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        if len (variable['dims']) == 1:  # if we have a 1d variable (time)
            ncVar_data[:] = variable['data'][:]

        
        elif len (variable['dims']) == 2:
            if variable['dims'][0] in ['nMesh2_data_time']:  # we have a 2D variable here (time, face)
                nrec = variable['data'].shape[0]  # number of records in unlimited dimension
                for t in xrange(nrec):
                    ncVar_data[t, ...] = variable['data'][t, ...]
            else:
                ncVar_data[:] = variable['data'][:]
        

        elif len (variable['dims']) == 3:
            if variable['dims'][0] in ['nMesh2_data_time']:  # todo: implement for "nMesh2_time"
                if variable['dims'][1] in ['nMesh2_layer_3d']:  # we have a 3D variable here (time, layer, face)
                    nrec = variable['data'].shape[0]  # number of records in unlimited dimension
                    nlay = variable['data'].shape[1]  # number of layers
                    for t in xrange(nrec):
                        for layer in xrange(nlay):
                            ncVar_data[t, layer, ...] = variable['data'][t, layer, ...]
                
                elif variable['dims'][1] in ['nMesh2_layer_2d']:  # if we have fake 3D variable here (time, 1, face)
                    nrec = variable['data'].shape[0]  # number of records in unlimited dimension
                    nlay = ncVar_data.shape[1]  # number of layers. Should be 1
                    if nlay != 1:  # since it is 2D, nlay should be always nlay=1
                        print '\n\n'+'ERROR\n'+'*'*50+'\n'+'declared 2D variable appeared to be 3D. Skipping...'+'\n'+'*'*50+'\n'
                        return
                    for t in xrange(nrec):
                        ncVar_data[t, 0, ...] = variable['data'][t,  ...]
        

        elif len (variable['dims']) == 4:
            if variable['dims'][0] in ['nMesh2_data_time']:
                if variable['dims'][1] in ['nMesh2_layer_3d']:
                    if variable['dims'][3] in ['two']:    # we have a 4D-bound variable (time, layer, face, 2)
                        nrec = variable['data'].shape[0]  # number of records in unlimited dimension
                        nlay = variable['data'].shape[1]  # number of layers
                        two  = variable['data'].shape[3]  # should be 2
                        for t in xrange(nrec):
                            for layer in xrange(nlay):
                                for zero_or_one in xrange(two):
                                    ncVar_data[t, layer, ..., zero_or_one] = variable['data'][t, layer, ..., zero_or_one]
                    
                    else:  #NOT IMPLEMENTED YET
                        pass

        elif len(variable['dims']) > 4:
            a = raw_input('Unsupported number of dimensions nDim>4.   Not yet implemented')
            pass


    if log: print _n, 'output    datashape:', ncVar_data.shape
    if log: print _n, 'Variable appended succesfully: %s' % (variable['vname'])
    root_grp.close()













def fill_bound_variable(Mesh2_var_x_bnd, Mesh2_var_y_bnd, Mesh2_face_nodes=None, Mesh2_node_lat=None, Mesh2_node_lon=None, log=False):
    """
    WARNING; THIS FUNCTION IS DEPRECATED! ALREADY BUILT IN MESH2 !!!

    Function calculates bounds for the elements. Based on indxing
    
    !!! WARNING: this routine will work only if starting index is 0 !!!
    !!! these lines >>> Mesh2_node_lon[dataIndexes_i[j]]            !!!
    !!! these lines >>> Mesh2_node_lat[dataIndexes_i[j]]            !!!


    int Mesh2_var_x_bnd (nMesh2_data, number_of_elements_in_data)
    int Mesh2_var_y_bnd (nMesh2_data, number_of_elements_in_data)
        example:
            int Mesh2_face_x_bnd (nMesh2_face, nMaxMesh2_face_nodes)
            int Mesh2_face_y_bnd (nMesh2_face, nMaxMesh2_face_nodes)
            int Mesh2_edge_x_bnd (nMesh2_edge, two)
            int Mesh2_edge_y_bnd (nMesh2_edge, two)
    """
    print 'filling bound variable'
    nData = Mesh2_var_x_bnd.shape[0]
    nElem = Mesh2_var_x_bnd.shape[1]



    if Mesh2_face_nodes is not None:
        print 'Filling  face bounds...'
        pass
    if Mesh2_node_lat is not None and Mesh2_node_lon is not None:
        print 'Filling  Mesh2_face_lon_bnd, Mesh2_face_lat_bnd'
        for i in xrange(nData):  #cycle over data(faces/edges)....
            # i is the index of face/edge...
            dataIndexes_i  = Mesh2_face_nodes[i, :]
            lon_elements_i = np.zeros(nElem)  # arrays with coordinates
            lat_elements_i = np.zeros(nElem)  # arrays with coordinates

            for j in xrange(nElem):
                lon_elements_i[j] = Mesh2_node_lon[dataIndexes_i[j]]
                lat_elements_i[j] = Mesh2_node_lat[dataIndexes_i[j]]
            if log:
                print '\t face =', i
                print '\t face_nodes[{0}] ='.format(i), dataIndexes_i
                print '\t face_x_bnd[{0}] ='.format(i), lon_elements_i
                print '\t face_y_bnd[{0}] ='.format(i), lat_elements_i
            Mesh2_var_x_bnd[i, :] = lon_elements_i
            Mesh2_var_y_bnd[i, :] = lat_elements_i
    print 'bounds have been filled'













def append_sigma_vertical_coord_vars(list_with_filenames, nLayers, filename, add_eta=False, add_depth=False, mask=None, log=False):
    '''
        function will append following variables to passed file...
            - Mesh2_sigma_layers
            - Mesh2_face_Wasserstand_2d (optional)  - added if <add_eta>=True
            - Mesh2_face_depth_2d (optional)  - added if <add_depth>=True
        

        add_eta, add_depth = False , it means that these vars won't be created
        why?
        because, most likely they are already appended or will be appended to file based on DICTIONARY4
    '''
    if log:
        print '-'*25+'\n appending sigma_vertical coordinates ...'
    # ----------------------------------------
    # ----------------------------------------
    # ---------------  SIGMA   ---------------
    # ----------------------------------------
    # ----------------------------------------
    var = dict()
    var['vname'] = 'nMesh2_layer_3d'
    var['dtype'] = 'f8'
    var['dims'] = ('nMesh2_layer_3d', )
    var['_FillValue'] = False

    ATTRS = dict()
    ATTRS['standard_name'] = 'ocean_sigma_coordinate'
    ATTRS['long_name'] = "sigma at layer midpoints"
    ATTRS['positive'] = 'up'
    ATTRS['formula_terms'] = 'sigma: nMesh2_layer_3d eta: Mesh2_face_Wasserstand_2d depth: Mesh2_face_depth_2d'

    var['attributes'] = ATTRS
    sigma, sigma_type = process_mossco_netcdf.get_sigma_coordinates(list_with_filenames, nLayers, varname='level')
    if sigma_type == 'center':
        pass
    elif sigma_type == 'border':
        sigma = process_mixed_data.create_sigma_coords_of_layer_center(sigma)

    var['data'] = sigma
    append_VariableData_to_netcdf(filename, var, log=log)
    del var

    # ----------------------------------------
    # ----------------------------------------
    # ---------------  ETA   -----------------
    # ----------------------------------------
    # ----------------------------------------
    if add_eta:  # actually append

        var = dict()
        var['vname'] = 'Mesh2_face_Wasserstand_2d'
        var['dtype'] = 'f4'
        var['dims'] = ('nMesh2_data_time', 'nMesh2_face')
        var['_FillValue'] = False

        ATTRS = dict()
        ATTRS['long_name']        = "Wasserstand, Face (Polygon)"
        ATTRS['standard_name']    = 'sea_surface_height'
        ATTRS['units']            = 'm'
        ATTRS['name_id']          = 3
        ATTRS['cell_measures']    = 'area: Mesh2_face_wet_area'
        ATTRS['cell_measures']    = 'nMesh2_data_time: point area: mean'
        ATTRS['coordinates']      = 'Mesh2_face_x Mesh2_face_y Mesh2_face_lon Mesh2_face_lat'
        ATTRS['grid_mapping']     = 'Mesh2_crs'
        ATTRS['mesh']             = 'Mesh2'
        ATTRS['location']         = 'face'
        
        var['attributes'] = ATTRS


        water_depth_at_soil_surface = None
        water_depth_at_soil_surface_vname = 'water_depth_at_soil_surface'
        bathymetry = None
        bathymetry_vname = 'bathymetry'

        for nc_file in list_with_filenames:
            root_grp = Dataset(nc_file, mode='r')
            if water_depth_at_soil_surface_vname in root_grp.variables.keys() and not water_depth_at_soil_surface:
                water_depth_at_soil_surface = process_mossco_netcdf.read_mossco_nc_rawvar(nc_file, water_depth_at_soil_surface_vname)
            
            if bathymetry_vname in root_grp.variables.keys() and not bathymetry:
                bathymetry = process_mossco_netcdf.read_mossco_nc_rawvar(nc_file, bathymetry_vname)

            root_grp.close()
            
            if water_depth_at_soil_surface is not None and bathymetry is not None:
                break


        water_level = process_mossco_netcdf.get_water_level(list_with_filenames, varname='water_level',
                            water_depth_at_soil_surface=water_depth_at_soil_surface, bathymetry=bathymetry,
                            log=log)
        var['data'] = process_mixed_data.flatten_xy_data(water_level, mask=mask)
        append_VariableData_to_netcdf(filename, var, log=log)
        del var

    # ----------------------------------------
    # ----------------------------------------
    # --------------  DEPTH   ----------------
    # ----------------------------------------
    # ----------------------------------------
    if add_depth:

        var = dict()
        var['vname'] = 'Mesh2_face_depth_2d'
        var['dtype'] = 'f8'
        var['dims'] = ('nMesh2_time', 'nMesh2_face')
        var['_FillValue'] = False

        ATTRS = dict()
        ATTRS['long_name']        = "Topographie"
        ATTRS['standard_name']    = 'sea_floor_depth_below_geoid'
        ATTRS['units']            = 'm'
        ATTRS['name_id']          = 17
        ATTRS['cell_measures']    = 'area: Mesh2_face_area'
        ATTRS['cell_measures']    = 'nMesh2_time: mean area: mean'
        ATTRS['coordinates']      = 'Mesh2_face_x Mesh2_face_y Mesh2_face_lon Mesh2_face_lat'
        ATTRS['grid_mapping']     = 'Mesh2_crs'
        ATTRS['mesh']             = 'Mesh2'
        ATTRS['location']         = 'face'

        var['attributes'] = ATTRS

        var['data'] = process_mixed_data.flatten_xy_data(bathymetry, mask=mask)
        append_VariableData_to_netcdf(filename, var, log=log)
        del var




    # ----------------------------------------
    # ----------------------------------------
    # ------------  ELEVATION   --------------
    # ----------------------------------------
    
    # ----------------------------------------
    # FACE middle values.....
    # ----------------------------------------
    var = dict()
    var['vname'] = 'Mesh2_face_z_face_3d'
    var['dtype'] = 'f4'
    var['dims']  = ('nMesh2_data_time', 'nMesh2_layer_3d', 'nMesh2_face')
    var['_FillValue'] = None
    
    ATTRS = dict()
    ATTRS['long_name']     = 'z_face [ face ]'
    ATTRS['units']         = 'm'
    ATTRS['positive']      = 'up'
    ATTRS['name_id']       = 1702
    ATTRS['bounds']        = 'Mesh2_face_z_face_bnd_3d'
    ATTRS['standard_name'] = 'depth'
    var['attributes']      = ATTRS
   

    var1 = dict()
    var1['vname'] = 'Mesh2_face_z_face_bnd_3d'
    var1['dtype'] = 'f4'
    var1['dims']  = ('nMesh2_data_time', 'nMesh2_layer_3d', 'nMesh2_face', 'two')
    var1['_FillValue'] = None
    
    ATTRS = dict()
    ATTRS['long_name'] = 'elevations of lower and upper layer-boundaries'
    ATTRS['units']     = 'm'
    var1['attributes'] = ATTRS





    elev, elev_borders = process_mixed_data.create_layer_elevation_from_sigma_coords(water_level, sigma, bathymetry, log=log)
    

    var['data']  = process_mixed_data.flatten_xy_data(elev, mask=mask)
    var1['data'] = process_mixed_data.flatten_xy_data(elev_borders, mask=mask)
    
    append_VariableData_to_netcdf(filename, var,  log=log)
    append_VariableData_to_netcdf(filename, var1, log=log)
    del var, var1
    



    if log:
        print '-'*25+'\n'
