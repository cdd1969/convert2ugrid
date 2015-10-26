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

# use this if you want to include modules from a subfolder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0], "lib")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)
    


def make_2d_qudratic_grid_or_curvilinear(x, y, mask=None, log=True, startingindex=1, data_location='T_points'):
    '''
    -----------------------------------------------
    Function creates unstructured grid of size.
    -----------------------------------------------
    See DOC string at Mesh2.Grid2D()
    -----------------------------------------------
    Inputs:
    -----------------------------------------------
    x - numpy array, x-coordinates of T or X points. Can be 1d or 2d
    y - numpy array, y-coordinates of T or X points. Can be 1d or 2d


    mask - is an 2D array of shape
            - (len(x), len(y)), if x and y are 1D
            - shape(x)=shape(y), if x and y are 2D
        filled with True, False (True - for masked cell, False - for non-masked)

    startingindex - an integer [0 , 1], for starting indes of elements.
                    startingindex=0 >>> 0, 1, 2, ... N-1
                    startingindex=1 >>> 1, 2, 3, ... N

    log - flag for showing console log
    '''
    #nx, ny = len(x_vector), len(y_vector)
    
    grid = Mesh2.Grid2D(x, y, mask=mask, data_location=data_location)


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
        #print 'Grid: {0}x{1} nodes'.format(nx, ny)
        print grid
        print '-'*80

    # ------------------------------------------------------------
    # ------------------------------------------------------------
    # ------------------------------------------------------------
    if startingindex == 1:
        # Do nothing, cause by default starting index is 1
        # WARNING!!! change 'start_index' attribute to = 1 in process_davit_ncdf.create_uGrid_ncdf() manually
        # for all variables like node_x, node_y, node_lat, etc...
        
        if log: print "WARNING! make sure 'start_index' attribute is set to <1> in process_davit_ncdf.create_uGrid_ncdf()"
        pass
    # ------------------------------------------------------------
    # ------------------------------------------------------------
    # ------------------------------------------------------------
    elif startingindex == 0:
        # use python-style indexing. Convert from [1..N] to [0..N-1]
        # WARNING!!! change 'start_index' attribute to = 0 in process_davit_ncdf.create_uGrid_ncdf() manually
        # for all variables like node_x, node_y, node_lat, etc...

        if log: print "WARNING! make sure 'start_index' attribute is set to <0> in process_davit_ncdf.create_uGrid_ncdf()"


        def find_max_indixes():
            max_i_n = 0
            max_i_e = 0
            max_i_f = 0
            for ne in xrange(nMesh2_edge):
                for i in xrange(2):
                    if Mesh2_edge_nodes[ne, i] > max_i_n:
                        max_i_n = int(Mesh2_edge_nodes[ne, i])
                    if Mesh2_edge_faces[ne, i] > max_i_f:
                        max_i_f = int(Mesh2_edge_faces[ne, i])
            for nf in xrange(nMesh2_face):
                for nfn in xrange(nMaxMesh2_face_nodes):
                    if Mesh2_face_nodes[nf, nfn] > max_i_n:
                        max_i_n = int(Mesh2_face_nodes[nf, nfn])
                    if Mesh2_face_edges[nf, nfn] > max_i_e:
                        max_i_e = int(Mesh2_face_edges[nf, nfn])
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
