#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# This file is part of "convert2ugrid" tool
#
# Author: Nikolai Chernikov, nikolai.chernikov.ru@gmail.com


import os
import sys
import inspect
import getopt
import process_davit_ncdf
import process_mossco_netcdf
import Mesh2
from process_mossco_netcdf import make_mask_array_from_mossco_bathymetry, read_mossco_lon_vector, read_mossco_lat_vector
# use this if you want to include modules from a subfolder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0], "lib")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)
    


def make_2d_rectangular_grid(x_vector, y_vector, mask=None, log=False, startingindex=1):
    '''
    function creates grid of size len(x_vector)*len(y_vector)

    x_vector - numpy array, vector of x-coordinates of cell centers
    y_vector - numpy array, vector of y-coordinates of cell centers

    mask - is an 2D array of size (len(x_vector), len(y_vector)), filled with True, False
            (True - for masked cell, False - for non-masked)
    startingindex - an integer [0 , 1], for starting indes of elements.
                    startingindex=0 >>> 0, 1, 2, ... N-1
                    startingindex=1 >>> 1, 2, 3, ... N

    log - flag for showing console log
    '''
    nx, ny = len(x_vector), len(y_vector)
    
    grid = Mesh2.Grid2D(x_vector, y_vector, mask=mask)


    #dims, topo, nodes, edges, faces
    nMesh2_node = grid.get_nMesh2_node()
    nMesh2_edge = grid.get_nMesh2_edge()
    nMesh2_face = grid.get_nMesh2_face()
    nMaxMesh2_face_nodes = grid.get_nMaxMesh2_face_nodes()

    Mesh2_edge_nodes = grid.get_Mesh2_edge_nodes()
    Mesh2_edge_faces = grid.get_Mesh2_edge_faces()
    Mesh2_face_nodes = grid.get_Mesh2_face_nodes()
    Mesh2_face_edges = grid.get_Mesh2_face_edges()

    Mesh2_node_x = grid.get_Mesh2_node_x()
    Mesh2_node_y = grid.get_Mesh2_node_y()
    
    Mesh2_edge_x = grid.get_Mesh2_edge_x()
    Mesh2_edge_y = grid.get_Mesh2_edge_y()
    
    Mesh2_face_x = grid.get_Mesh2_face_x()
    Mesh2_face_y = grid.get_Mesh2_face_y()
    Mesh2_face_center_x = grid.get_Mesh2_face_center_x()
    Mesh2_face_center_y = grid.get_Mesh2_face_center_y()
    Mesh2_face_area = grid.get_Mesh2_face_area()

    Mesh2_edge_x_bnd, Mesh2_edge_y_bnd = grid.get_Mesh2_edge_bnd()
    Mesh2_face_x_bnd, Mesh2_face_y_bnd = grid.get_Mesh2_face_bnd()
    
    if log:
        print '-'*80
        print 'Grid: {0}x{1} nodes'.format(nx, ny)
        print grid
        print '-'*80

    # ------------------------------------------------------------
    # ------------------------------------------------------------
    # ------------------------------------------------------------
    if startingindex == 1:
        # Do nothing, cause by default starting index is 1
        # WARNING!!! change 'start_index' attribute to = 1 in process_davit_ncdf.create_uGrid_ncdf() manually
        # for all variables like node_x, node_y, node_lat, etc...
        
        if log: print "WARNING!!! change 'start_index' attribute to = 1 in process_davit_ncdf.create_uGrid_ncdf() manually"
        if log: print "for all variables like node_x, node_y, node_lat, etc..."
        pass
    # ------------------------------------------------------------
    # ------------------------------------------------------------
    # ------------------------------------------------------------
    elif startingindex == 0:
        # use python-style indexing. Convert from [1..N] to [0..N-1]
        # WARNING!!! change 'start_index' attribute to = 0 in process_davit_ncdf.create_uGrid_ncdf() manually
        # for all variables like node_x, node_y, node_lat, etc...

        if log: print "WARNING!!! change 'start_index' attribute to = 0 in process_davit_ncdf.create_uGrid_ncdf() manually"
        if log: print "for all variables like node_x, node_y, node_lat, etc..."


        def find_max_indixes():
            max_i_n = 0
            max_i_e = 0
            max_i_f = 0
            for ne in xrange(nMesh2_edge):
                for i in xrange(2):
                    if Mesh2_edge_nodes[ne, i] > max_i_n:
                        max_i_n = Mesh2_edge_nodes[ne, i]
                    if Mesh2_edge_faces[ne, i] > max_i_f:
                        max_i_f = Mesh2_edge_faces[ne, i]
            for nf in xrange(nMesh2_face):
                for nfn in xrange(nMaxMesh2_face_nodes):
                    if Mesh2_face_nodes[nf, nfn] > max_i_n:
                        max_i_n = Mesh2_face_nodes[nf, nfn]
                    if Mesh2_face_edges[nf, nfn] > max_i_e:
                        max_i_e = Mesh2_face_edges[nf, nfn]
            print 'maximum indexes of Nodes, Edges, Faces:', max_i_n, max_i_e, max_i_f

        if log:
            print '-'*50
            # find maximum indexes used
            find_max_indixes()

        for ne in xrange(nMesh2_edge):
            for i in xrange(2):
                if Mesh2_edge_nodes[ne, i] != -999:
                    Mesh2_edge_nodes[ne, i] += -1
                if Mesh2_edge_faces[ne, i] != -999:
                    Mesh2_edge_faces[ne, i] += -1

        for nf in xrange(nMesh2_face):
            for nfn in xrange(nMaxMesh2_face_nodes):
                if Mesh2_face_nodes[nf, nfn] != -999:
                    Mesh2_face_nodes[nf, nfn] += -1
                if Mesh2_face_edges[nf, nfn] != -999:
                    Mesh2_face_edges[nf, nfn] += -1

        if log:
            # find maximum indexes used
            print '\t indexes have been reworked'
            find_max_indixes()
            print '-'*50

    
    
    # ------------------------------------------------------------
    # ------------------------------------------------------------
    # ------------------------------------------------------------
    else:
        raise ValueError('{1}Passed starting index for Mesh2 elements is {0}. Choose from [0, 1]{1}'.format(startingindex, '\n'+'-'*50+'\n'))
    
    # ------------------------------------------------------------
    # ------------------------------------------------------------
    # ------------------------------------------------------------
    return [nMesh2_node, nMesh2_edge, nMesh2_face, nMaxMesh2_face_nodes],\
           [Mesh2_edge_nodes, Mesh2_edge_faces, Mesh2_face_nodes, Mesh2_face_edges],\
           [Mesh2_node_x, Mesh2_node_y],\
           [Mesh2_edge_x, Mesh2_edge_y],\
           [Mesh2_face_x, Mesh2_face_y, Mesh2_face_center_x, Mesh2_face_center_y, Mesh2_face_area],\
           [Mesh2_edge_x_bnd, Mesh2_edge_y_bnd, Mesh2_face_x_bnd, Mesh2_face_y_bnd]


def test():
    nc_name1 = 'deep_lake-1x1-reference_3d.2d.0000.nc'
    m = make_mask_array_from_mossco_bathymetry(nc_name1)
    x = read_mossco_lon_vector(nc_name1)
    y = read_mossco_lat_vector(nc_name1)
    make_2d_rectangular_grid(x, y, mask=m, log=True)
    

def test2():
        #filenameINPUT1 = 'data/NSBS/netcdf_reference_3d.nc'
        filenameINPUT1 = '../data/square/netcdf_reference_3d.nc'
        #filenameINPUT2 = 'data/NSBS/topo.nc'
        filenameINPUT2 = '../data/square/deep_lake-1x1-reference_3d.2d.0000.nc'
        # --------------------------------------------------
        # 0) Read, x any y vectors from the netcdf
        # --------------------------------------------------
        #x_vector, _ = process_mossco_netcdf.read_mossco_nc_1d(filenameINPUT2, 'lon')
        #y_vector, _ = process_mossco_netcdf.read_mossco_nc_1d(filenameINPUT2, 'lat')
        x_vector, _ = process_mossco_netcdf.read_mossco_nc_1d(filenameINPUT2, 'lonc')
        y_vector, _ = process_mossco_netcdf.read_mossco_nc_1d(filenameINPUT2, 'latc')
        print '-'*50
        for i in x_vector: print i
        print '-'*50
        for i in y_vector: print i
        print '-'*50
        m = process_mossco_netcdf.make_mask_array_from_mossco_bathymetry(filenameINPUT2, varname='bathymetry', fillvalue=-10., transpose=True)

        # --------------------------------------------------
        # 1) Create grid, and unpack values
        # --------------------------------------------------
        g = Mesh2.Grid2D(x_vector, y_vector, mask=None)
        n = g.get_nodes()
        for i in range(36):
            print n[i]
        """
        dims, topo, local_nodes, local_edges, local_faces = make_2d_rectangular_grid(x_vector, y_vector, mask=None, log=True)

        nMesh2_node = dims[0]
        nMesh2_edge = dims[1]
        nMesh2_face = dims[2]
        nMaxMesh2_face_nodes = dims[3]

        Mesh2_edge_nodes = topo[0]
        Mesh2_edge_faces = topo[1]
        Mesh2_face_nodes = topo[2]
        Mesh2_face_edges = topo[3]

        Mesh2_node_x = local_nodes[0]
        Mesh2_node_y = local_nodes[1]

        Mesh2_edge_x = local_edges[0]
        Mesh2_edge_y = local_edges[1]

        Mesh2_face_x = local_faces[0]
        Mesh2_face_y = local_faces[1]
        Mesh2_face_center_x = local_faces[2]
        Mesh2_face_center_y = local_faces[3]
        """

if __name__ == '__main__':
    test2()
