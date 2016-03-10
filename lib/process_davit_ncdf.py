#!/usr/bin/python
#-*- coding: utf-8 -*-
# Copyright (C) 2015, Nikolai Chernikov <nikolai.chernikov.ru@gmail.com>
#
# "convert2ugrid" is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License v3+. "convert2ugrid" is distributed in the
# hope that it will be useful, but WITHOUT ANY WARRANTY.  Consult the file
# LICENSE.GPL or www.gnu.org/licenses/gpl-3.0.txt for the full license terms.


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
import process_cdl
from . import ui


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
                        start_index=0,
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
        ncVar_Mesh2_crs = root_grp.createVariable('Mesh2_crs', 'i', (), fill_value=False)
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
    ncVar_Mesh2_edge_nodes.start_index = start_index
    ncVar_Mesh2_edge_nodes[:] = Mesh2_edge_nodes[:]

    # --------------------------------------------------
    #                                   3.2) EDGE >>> FACES
    # --------------------------------------------------

    ncVar_Mesh2_edge_faces = root_grp.createVariable('Mesh2_edge_faces', 'i', ('nMesh2_edge', 'two'), fill_value=-999)
    ncVar_Mesh2_edge_faces.long_name = 'Face- (Polygon-) Verzeichnis der Kanten, linker und rechter Nachbar'
    ncVar_Mesh2_edge_faces.cf_role = 'edge_face_connectivity'
    ncVar_Mesh2_edge_faces.start_index = start_index
    ncVar_Mesh2_edge_faces[:] = Mesh2_edge_faces[:]

    # --------------------------------------------------
    #                                   3.3) FACE >>> NODES
    # --------------------------------------------------

    ncVar_Mesh2_face_nodes = root_grp.createVariable('Mesh2_face_nodes', 'i', ('nMesh2_face', 'nMaxMesh2_face_nodes'), fill_value=-999)
    ncVar_Mesh2_face_nodes.long_name = 'Knotenverzeichnis der Faces (Polygone),entgegen dem Uhrzeigersinn'
    ncVar_Mesh2_face_nodes.cf_role = 'face_node_connectivity'
    ncVar_Mesh2_face_nodes.start_index = start_index
    ncVar_Mesh2_face_nodes[:] = Mesh2_face_nodes[:]

    # ------------------------------------------------------------------------------------
    #                                   3.4) FACE >>> EDGES
    # ------------------------------------------------------------------------------------

    ncVar_Mesh2_face_edges = root_grp.createVariable('Mesh2_face_edges', 'i', ('nMesh2_face', 'nMaxMesh2_face_nodes'), fill_value=-999)
    ncVar_Mesh2_face_edges.long_name = 'Kantenverzeichnis der Faces (Polygone),entgegen dem Uhrzeigersinn'
    ncVar_Mesh2_face_edges.cf_role = 'face_edge_connectivity'
    ncVar_Mesh2_face_edges.start_index = start_index
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
    _nc = Dataset(fname_davit, mode='r')
    try:
        nt = _nc.dimensions['nMesh2_data_time'].__len__()
    except:
        nt = _nc.dimensions['nMesh2_time'].__len__()
    
    nz     = _nc.dimensions['nMesh2_layer_3d'].__len__()
    nFaces = _nc.dimensions['nMesh2_face'].__len__()
    nEdges = _nc.dimensions['nMesh2_edge'].__len__()
    _nc.close()

    if log:
        print 'nZ, nY, nX =', nz, ny, nx
        print 'nTimesteps =', nt
        print 'nFaces =', nFaces
        print 'nEdges =', nEdges


    # FACE middle values.....
    # ----------------------------------------------------------------------------
    var = process_cdl.cdlVariable()
    var.set_name('Mesh2_face_z_face_3d')
    var.set_dtype('float')
    var.set_dims(('nMesh2_data_time', 'nMesh2_layer_3d', 'nMesh2_face'))
    var.set_fillvalue(None)
    var.set_attr('long_name', 'z_face [ face ]')
    var.set_attr('units', 'm')
    var.set_attr('positive', 'down')
    var.set_attr('name_id', 1702)
    var.set_attr('bounds', 'Mesh2_face_z_face_bnd_3d')
    var.set_attr('standard_name', 'depth')
    var.set_attr('comment', "warning: dummy values are used")

    davit_dummy = np.zeros((nt, nz, nFaces))
    for z in xrange(nz):
        davit_dummy[:, z, :] = nz-z-0.5  # this equation will produce values like 0.5, 1.5, 2.5, ..., nz-1.0-0.5
    

    var_data = davit_dummy
    append_VariableData_to_netcdf(fname_davit, var, var_data, fv=var.get_fillvalue(), log=log)
    del var

    # FACE bounds values...
    # ----------------------------------------------------------------------------
    if log: print '{0}{0}{1}\n{0}{0}'.format('+'*50+'\n', 'Face Bounds')

    var = process_cdl.cdlVariable()
    var.set_name('Mesh2_face_z_face_bnd_3d')
    var.set_dtype('float')
    var.set_dims(('nMesh2_data_time', 'nMesh2_layer_3d', 'nMesh2_face', 'two'))
    var.set_fillvalue(None)
    var.set_attr('name_id', 1703)
    #var.set_attr('units', 'm')
    var.set_attr('comment', "warning: dummy values are used")

    davit_dummy = np.zeros((nt, nz, nFaces, 2))
    for z in xrange(nz):
        davit_dummy[:, z, :, 0] = nz-z-1.  # this equation will produce values like 0, 1, 2, ... nz-z
        davit_dummy[:, z, :, 1] = nz-z+0.  # this equation will produce values like 1, 2, 3, ... nz-z+1


    var_data = davit_dummy
    append_VariableData_to_netcdf(fname_davit, var, var_data, fv=var.get_fillvalue(), log=log)
    del var


    if log: print '{0}{0}{1}\n{0}{0}'.format('+'*50+'\n', 'EDGE')

    # EDGE middle values.....
    # ----------------------------------------------------------------------------
    var = process_cdl.cdlVariable()
    var.set_name('Mesh2_edge_z_edge_3d')
    var.set_dtype('float')
    var.set_dims(('nMesh2_data_time', 'nMesh2_layer_3d', 'nMesh2_edge'))
    var.set_fillvalue(None)
    var.set_attr('long_name', 'z_edge [ edge ]')
    var.set_attr('units', 'm')
    var.set_attr('positive', 'down')
    var.set_attr('name_id', 1707)
    var.set_attr('bounds', 'Mesh2_edge_z_edge_bnd_3d')
    var.set_attr('standard_name', 'depth')
    var.set_attr('comment', "warning: dummy values are used")

    davit_dummy = np.zeros((nt, nz, nEdges))
    for z in xrange(nz):
        davit_dummy[:, z, :] = nz-z-0.5  # this equation will produce values like 0.5, 1.5, 2.5, ..., nz-1.0-0.5
    

    var_data = davit_dummy
    append_VariableData_to_netcdf(fname_davit, var, var_data, fv=var.get_fillvalue(), log=log)
    del var

    # EDGE bounds values...
    # ----------------------------------------------------------------------------
    if log: print '{0}{0}{1}\n{0}{0}'.format('+'*50+'\n', 'EDGE Bounds')

    var = process_cdl.cdlVariable()
    var.set_name('Mesh2_edge_z_edge_bnd_3d')
    var.set_dtype('float')
    var.set_dims(('nMesh2_data_time', 'nMesh2_layer_3d', 'nMesh2_edge', 'two'))
    var.set_fillvalue(None)
    var.set_attr('name_id', 1708)
    #var.set_attr('units', 'm')
    var.set_attr('comment', "warning: dummy values are used")
    
    davit_dummy = np.zeros((nt, nz, nEdges, 2))

    for z in xrange(nz):
        davit_dummy[:, z, :, 0] = nz-z-1.  # this equation will produce values like 0, 1, 2, ... nz-z
        davit_dummy[:, z, :, 1] = nz-z+0.  # this equation will produce values like 1, 2, 3, ... nz-z+1


    var_data = davit_dummy
    append_VariableData_to_netcdf(fname_davit, var, var_data, fv=var.get_fillvalue(), log=log)
    del var
    if log: print 'finished: append_test_Mesh2_face_z_3d_and_Mesh2_face_z_3d_bnd()'



def append_Time_andDatetime_to_netcdf(fname_davitnc, fname_mossconc=None, time_var_name='time', dummy_values=False, log=True):
        # --------------------------------------------------
        # 3.1) fill netcdf with time variables (TIME and DATATIME)
        # --------------------------------------------------
        if dummy_values is False:

            _nc = Dataset(fname_mossconc, mode='r')
            time_var_units = _nc.variables[time_var_name].units
            max_t          = _nc.variables[time_var_name][:].max()
            min_t          = _nc.variables[time_var_name][:].min()
            _nc.close()
        else:
            time_var_units = 'seconds since 2000-01-01 00:00:00'
            max_t = 60*60*24
            min_t = 0



        var = process_cdl.cdlVariable()
        var.set_name('nMesh2_time')
        var.set_dtype('double')
        var.set_dims(('nMesh2_time', ))
        var.set_fillvalue(False)
        var.set_attr('long_name', "time")
        # because currently MOSSCO writes a time output as "hours since 2009-01-02 00:00:00",
        # it is nessesary to modify it. Will be removed, after corrections in MOSSCO.
        var.set_attr('units', time_var_units+' 01:00')  # '01:00 is indicating timezone, see CF conventions for details'
        var.set_attr('name_id', 1640)
        var.set_attr('axis', "T")
        var.set_attr('bounds', "nMesh2_time_bnd")
        var.set_attr('calendar', "gregorian")
        var.set_attr('standard_name', "time")
        if dummy_values:
            var.set_attr('comment', "warning: dummy values are used")


        var_data = np.array([min_t])
        append_VariableData_to_netcdf(fname_davitnc, var, var_data, fv=var.get_fillvalue(), log=log)
        del var
        

        var = process_cdl.cdlVariable()
        var.set_name('nMesh2_time_bnd')
        var.set_dtype('double')
        var.set_dims(('nMesh2_time', 'two'))
        var.set_fillvalue(False)
        if dummy_values:
            var.set_attr('comment', "warning: dummy values are used")

        var_data = np.array([min_t, max_t])
        append_VariableData_to_netcdf(fname_davitnc, var, var_data, fv=var.get_fillvalue(), log=log)
        del var


        if dummy_values is False:
            var = process_cdl.cdlVariable()
            var.set_name('nMesh2_data_time')
            var.set_dtype('double')
            var.set_dims(('nMesh2_data_time', ))
            var.set_fillvalue(False)
            var.set_attr('long_name', "time")
            var.set_attr('units', time_var_units + ' 01:00')
            var.set_attr('name_id', 1640)
            var.set_attr('axis', "T")
            var.set_attr('calendar', "gregorian")
            var.set_attr('standard_name', "time")
            var_data = process_mossco_netcdf.read_mossco_nc_1d(fname_mossconc, 'time')
            append_VariableData_to_netcdf(fname_davitnc, var, var_data, fv=var.get_fillvalue(), log=log)
            del var





def append_VariableData_to_netcdf(nc_out_fname, var, var_data, fv=None, log=True):
    ''' Procedure appends variable (meta+data) to a `nc_out_fname` netCDF file.
    Args:
    -----
        nc_out_fname (str):
            filename of a netcdf file to which variable `var` will be appended

        var (process_cdl.cdlVariable|process_mixed_data.cdlVariableExt):
            meta information of the variable which will be appended

        var_data (nd-array):
            n-dimensional numpy array with variable data. It should has shape,
            that corresponds with dimensions of netcdf file `nc_out_fname` declared
            with `var.get_dims()`

        fv (None or float|int):
            fill value. Default: None

    '''
    _n = 'append_VariableData_to_netcdf():'
    if isinstance(var, process_cdl.cdlVariable):
        meta = var
    elif isinstance(var, process_mixed_data.cdlVariableExt):
        meta = var.get_parent()
    else:
        raise TypeError(_n+'Invalid datatype. Passed argument `var` is of type {0}. Should be <cdlVariable> or <cdlVariableExt>'.format(type(var)))

    if log: print _n, 'Now try to append variable <{0}>'.format(meta.get_name())
    # --------------------------------------------------
    #                   Appending ncdf
    # --------------------------------------------------
    nc = Dataset(nc_out_fname, mode='a')

    # check if this variable already exists!
    if meta.get_name() in nc.variables.keys():
        if log: print _n, 'Variable <{0}> skipped. Already exists in file <{1}>'.format(meta.get_name(), nc_out_fname)
        nc.close()
        return

    # now compare shape of the array allocated within netcdf and shape of passed data
    nc_var_shape = tuple([nc.dimensions[d].__len__() for d in meta.get_dims()])
    # first check data shape
    if var_data is not None:
        if nc_var_shape != var_data.shape:
            # what could have gone wrong?
            # 1) unlimited dimension has not been yet initialized!
            #    (t, y, x) where <t> is UNLIMITED and NO data has been specified yet: (0, ny, nx) where <ny> and <nx> - integers.
            #    And lets imagine we try to pass array of shape (3, ny, nx) with 3 timesteps written... This error will trigger
            #    Let's get through!
            if len(nc_var_shape) != len(var_data.shape):
                msg = _n+' Invalid dimensions. Passed metadata {0} declares following dimension names {1} for variable {2}. NetCDF file {3} has following lengths of previously mentioned dimensions {4}. Therefore exacly this shape {4} will be allocated for an array within NetCDF file. However passed data-array has shape {5} and cannot be fit into allocated shape'.format(type(var), meta.get_dims(), meta.get_name(), nc_out_fname, nc_var_shape, var_data.shape)
                msg2 = _n+' Now I will squeeze (see numpy.squeeze()) passed data-array and try to fit it into allocated array within NetCDF file once again'
                print msg
                print msg2
                # here can be something like.... (1, 2) and (2, )
                squeezed_ncvarshape = tuple([d for d in nc_var_shape if d != 1])
                squeezed_vdatashape = tuple([d for d in var_data.shape if d != 1])
                if len(squeezed_ncvarshape) != len(squeezed_vdatashape):
                    print _n, _n+' Squeezed arrays do not match each other (allocated squeezed shape {0} != data-array squeezed shape {1}). Proceeding...'.format(squeezed_ncvarshape, squeezed_vdatashape)
                    ui.promt(_n+' Press Enter to skip appending variable <{0}>'.format(meta.get_name()), color='yellow', pause=True)
                    nc.close()
                    return
                else:
                    print _n+' Squeezed arrays match each other (allocated squeezed shape {0} = data-array squeezed shape {1}). Proceeding...'.format(squeezed_ncvarshape, squeezed_vdatashape)
                    #if not ui.promtYesNo(_n+' Continue appending variable `{0}` ("yes") or skip it ("no")?'.format(meta.get_name())):
                    #    nc.close()
                    #    return

            else:  # if shape is equal
                for i, dim_length_nc, dim_length_data in zip(xrange(len(nc_var_shape)), nc_var_shape, var_data.shape):
                    if (dim_length_nc != dim_length_data and dim_length_nc == 0 and nc.dimensions[meta.get_dims()[i]].isunlimited()):
                        # no error, since unlimited dimension of undefined length
                        pass
                    else:
                        print msg
                        ui.promt(_n+' Press Enter to skip appending variable <{0}>'.format(meta.get_name()), color='yellow', pause=True)
                        nc.close()
                        return

    # at this point all tests have been passed.....

    # now create Variable
    nc_var = nc.createVariable(meta.get_name(), meta.get_dtype(syntax='python-netcdf'), dimensions=meta.get_dims(), fill_value=fv)
    
    # add attributes
    if meta.get_attrs():
        for attr_name, attr_value in meta.get_attrs().iteritems():
            if not attr_name.startswith('_'):  # reserved for special attributes
                # treat real range values
                if attr_name in ['valid_range', 'range']:
                    attr_value = [np.float32(var_data.min()), np.float32(var_data.max())]
                elif attr_name in ['valid_max', 'max']:
                    attr_value = np.float32(var_data.max())
                elif attr_name in ['valid_min', 'min']:
                    attr_value = np.float32(var_data.min())

                if log: print '\t adding attribute <{0} = {1}> of type {2}'.format(attr_name, attr_value, type(attr_value))

                # actually append
                nc_var.setncattr(attr_name, attr_value)
    
    # fill data
    nc_var[:] = var_data

    if log: print _n, 'output datashape:', nc_var.shape
    if log: print _n, 'Variable appended succesfully: %s' % (meta.get_name())
    nc.close()




def append_sigma_vertical_coord_vars(list_with_filenames, nLayers, nc_out_fname, add_eta=False, add_depth=False, mask=None, sigma_varname='level', log=False):
    ''' Appends variables that define vertical position of layers...
            - nMesh2_layer_3d                       >>> 1d sigma coords of cell center
            - Mesh2_face_z_face_3d                  >>> elevation cell center values
            - Mesh2_face_z_face_3d_bnd              >>> elevation cell border values
            - (optionaly) Mesh2_face_Wasserstand_2d >>> water level
            - (optionaly) Mesh2_face_depth_2d       >>> bottom depth
        ... to passed `nc_out_fname` netcdf file
        
    Args:
    -----
        list_with_filenames (list of str):
            list with names of netcdf files. The var `sigma_varname` will be searched within this files
        nLayers (int):
            number of vertical layers
        nc_out_fname (str):
            filename of the output netcdf file
        add_eta (bool):
            flag to add variable "Mesh2_face_Wasserstand_2d" to file `nc_out_fname`.
            This is useful because, most likely this var is already appended or will be
            appended to file based on DICTIONARY4
        add_depth (bool):
            flag to add variable "Mesh2_face_depth_2d" to file `nc_out_fname`
            This is useful because, most likely this var is already appended or will be
            appended to file based on DICTIONARY4
        mask (2D array of bool):
            2d array of (y, x) dimensions with boolean mask (to treat NaN cells)
        sigma_varname (str):
            name of the varable to get sigma-layer info from
        log (bool):
            flag to print output
        
    '''
    if log:
        print '-'*25+'\n appending sigma_vertical coordinates ...'
    # ----------------------------------------
    # ----------------------------------------
    # ---------------  SIGMA   ---------------
    # ----------------------------------------
    # ----------------------------------------
    var = process_cdl.cdlVariable()
    var.set_name('nMesh2_layer_3d')
    var.set_dtype('double')
    var.set_dims(('nMesh2_layer_3d', ))
    var.set_fillvalue(False)
    var.set_attr('standard_name', 'ocean_sigma_coordinate')
    var.set_attr('long_name', "sigma at layer midpoints")
    var.set_attr('positive', 'up')
    var.set_attr('formula_terms', 'sigma: nMesh2_layer_3d eta: Mesh2_face_Wasserstand_2d depth: Mesh2_face_depth_2d')

    sigma, sigma_type = process_mossco_netcdf.get_sigma_coordinates(list_with_filenames, nLayers, sigma_varname=sigma_varname, waterdepth_varname='water_depth_at_soil_surface', layerdepth_varname='getmGrid3D_getm_layer',)
    if sigma_type == 'center':
        pass
    elif sigma_type == 'border':
        sigma = process_mixed_data.create_sigma_coords_of_layer_center(sigma)

    append_VariableData_to_netcdf(nc_out_fname, var, sigma, fv=var.get_fillvalue(), log=log)
    del var

    # ----------------------------------------
    # ----------------------------------------
    # ---------------  ETA   -----------------
    # ----------------------------------------
    # ----------------------------------------
    if add_eta:  # actually append

        var = process_cdl.cdlVariable()
        var.set_name('Mesh2_face_Wasserstand_2d')
        var.set_dtype('float')
        var.set_dims(('nMesh2_data_time', 'nMesh2_face'))
        var.set_fillvalue(False)
        var.set_attr('long_name', "Wasserstand, Face (Polygon)")
        var.set_attr('standard_name', 'sea_surface_height')
        var.set_attr('units', 'm')
        var.set_attr('name_id', 3)
        var.set_attr('cell_measures', 'area: Mesh2_face_wet_area')
        var.set_attr('cell_measures', 'nMesh2_data_time: point area: mean')
        var.set_attr('coordinates', 'Mesh2_face_x Mesh2_face_y Mesh2_face_lon Mesh2_face_lat')
        var.set_attr('grid_mapping', 'Mesh2_crs')
        var.set_attr('mesh', 'Mesh2')
        var.set_attr('location', 'face')


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
        var_data = process_mixed_data.flatten_xy_data(water_level, mask=mask)
        append_VariableData_to_netcdf(nc_out_fname, var, var_data, fv=var.get_fillvalue(), log=log)
        del var

    # ----------------------------------------
    # ----------------------------------------
    # --------------  DEPTH   ----------------
    # ----------------------------------------
    # ----------------------------------------
    if add_depth:

        var = process_cdl.cdlVariable()
        var.set_name('Mesh2_face_depth_2d')
        var.set_dtype('double')
        var.set_dims(('nMesh2_time', 'nMesh2_face'))
        var.set_fillvalue(False)
        var.set_attr('long_name', "Topographie")
        var.set_attr('standard_name', 'sea_floor_depth_below_geoid')
        var.set_attr('units', 'm')
        var.set_attr('name_id', 17)
        var.set_attr('cell_measures', 'area: Mesh2_face_area')
        var.set_attr('cell_measures', 'nMesh2_time: mean area: mean')
        var.set_attr('coordinates', 'Mesh2_face_x Mesh2_face_y Mesh2_face_lon Mesh2_face_lat')
        var.set_attr('grid_mapping', 'Mesh2_crs')
        var.set_attr('mesh', 'Mesh2')
        var.set_attr('location', 'face')

        var_data = process_mixed_data.flatten_xy_data(bathymetry, mask=mask)
        append_VariableData_to_netcdf(nc_out_fname, var, var_data, fv=var.get_fillvalue(), log=log)
        del var




    # ----------------------------------------
    # ----------------------------------------
    # ------------  ELEVATION   --------------
    # ----------------------------------------
    
    # ----------------------------------------
    # FACE cell center values.....
    # ----------------------------------------
    var1 = process_cdl.cdlVariable()
    var1.set_name('Mesh2_face_z_face_3d')
    var1.set_dtype('float')
    var1.set_dims(('nMesh2_data_time', 'nMesh2_layer_3d', 'nMesh2_face'))
    var1.set_fillvalue(None)
    var1.set_attr('long_name', 'z_face [ face ]')
    var1.set_attr('units', 'm')
    var1.set_attr('positive', 'up')
    var1.set_attr('name_id', 1702)
    var1.set_attr('bounds', 'Mesh2_face_z_face_bnd_3d')
    var1.set_attr('standard_name', 'depth')
   
    # ----------------------------------------
    # FACE cell border values.....
    # ----------------------------------------
    var2 = process_cdl.cdlVariable()
    var2.set_name('Mesh2_face_z_face_bnd_3d')
    var2.set_dtype('float')
    var2.set_dims(('nMesh2_data_time', 'nMesh2_layer_3d', 'nMesh2_face', 'two'))
    var2.set_fillvalue(None)
    var2.set_attr('long_name', 'elevations of lower and upper layer-boundaries')
    var2.set_attr('units', 'm')

    elev, elev_borders = process_mixed_data.create_layer_elevation_from_sigma_coords(water_level, sigma, bathymetry, log=log)
    

    var_data1  = process_mixed_data.flatten_xy_data(elev, mask=mask)
    var_data2 = process_mixed_data.flatten_xy_data(elev_borders, mask=mask)
    
    append_VariableData_to_netcdf(nc_out_fname, var1, var_data1, fv=var1.get_fillvalue(), log=log)
    append_VariableData_to_netcdf(nc_out_fname, var2, var_data2, fv=var2.get_fillvalue(), log=log)
    del var1, var2

    if log:
        print '-'*25+'\n'
