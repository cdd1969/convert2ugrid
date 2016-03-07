#!/usr/bin/python
#-*- coding: utf-8 -*-
# Copyright (C) 2015, Nikolai Chernikov <nikolai.chernikov.ru@gmail.com>
#
# "convert2ugrid" is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License v3+. "convert2ugrid" is distributed in the
# hope that it will be useful, but WITHOUT ANY WARRANTY.  Consult the file
# LICENSE.GPL or www.gnu.org/licenses/gpl-3.0.txt for the full license terms.

from __future__ import division
import numpy as np
import numpy.ma as ma


def gridhelp():
    errmsg = 'Currently limited grid-types are implemented:\n\
              \t a) <x_coord> and <y_coord> are 1D, at T-points, 2D <bathymetry> at T_points => Rectengular uniform\n\
              \t b) <x_coord> and <y_coord> are 1D, at X-points, 2D <bathymetry> at T_points => Rectengular (uniform/non-uniform)\n\
              \t c) <x_coord> and <y_coord> are 2D, at X-points, 2D <bathymetry> at T_points => Curvilinear'
    return errmsg


class Mesh2_node(object):
    """ Class describes a Point (Node) on a 2D plane"""
    def __init__(self, x, y, index):
        # ------------------
        # <x, y>  - floats, describing cartesian coordinates of a point
        # <index> - integer indicating id of current element
        self.__node_x = float(x)
        self.__node_y = float(y)

        self.__index = int(index)

    def set(self, x, y):
        self.__node_x = float(x)
        self.__node_y = float(y)
    
    # >>> getters >>>
    
    def get_node_x(self):
        return self.__node_x
    
    def get_node_y(self):
        return self.__node_y

    def get_vector(self):
        return np.array([self.__node_x, self.__node_y])

    def get_index(self):
        return self.__index

    def __eq__(self, node_to_compare):
        if not isinstance(node_to_compare, Mesh2_node):
            err_msg = '<Mesh2_edge> id={2}\n Type mismatch. Received object <{0}> of type <{1}>. Must be type <Mesh2_node>'.format(node_to_compare, type(node_to_compare), self.__index)
            raise TypeError(err_msg)
        if np.array_equal(self.get_vector(), node_to_compare.get_vector()):
            return True
        else:
            return False
    
    def __ne__(self, node_to_compare):
        if not isinstance(node_to_compare, Mesh2_node):
            err_msg = '<Mesh2_edge> id={2}\n Type mismatch. Received object <{0}> of type <{1}>. Must be type <Mesh2_node>'.format(node_to_compare, type(node_to_compare), self.__index)
            raise TypeError(err_msg)
        if not np.array_equal(self.get_vector(), node_to_compare.get_vector()):
            return True
        else:
            return False
    
    def __repr__(self):
        return ("Node: id={2} ( x: {0} , y: {1} )".format(self.__node_x, self.__node_y, self.__index))


class Mesh2_edge(object):
    """ Class describing an Edge connecting two Nodes (vertexes) on a 2D plane """
    def __init__(self, node1, node2, index):
        # ------------------
        # <node1, node2> - ire objects of type <Mesh2_node>
        # <index>        - integer indicating id of current element

        self.__node1 = node1
        self.__node2 = node2
        self.__index = int(index)
        self.update()
    
    def set(self, node1, node2):
        self.__node1 = node1
        self.__node2 = node2
        self.update()

    def set_index(self, index):
        self.__index = int(index)

    def update(self):
        self.validateInput()

    def validateInput(self):
        for item in [self.__node1, self.__node2]:
            if not isinstance(item, Mesh2_node):
                err_msg = '<Mesh2_edge> id={2}\nobject <{0}> of type <{1}> passed as input to <Mesh2_edge> must be of type <Mesh2_node>'.format(item, type(item), self.__index)
                raise TypeError(err_msg)
    
    def calculateGeometry(self):
        self.__vector = np.array([self.__node2.get_node_x() - self.__node1.get_node_x(), self.__node2.get_node_y() - self.__node1.get_node_y()])

    # >>> getters >>>

    def get_edge_x(self):
        return (self.__node1.get_node_x() + self.__node2.get_node_x()) / 2.
    
    def get_edge_y(self):
        return(self.__node1.get_node_y() + self.__node2.get_node_y()) / 2.

    def get_edge_nodes(self):
        return np.array([self.__node1.get_index(), self.__node2.get_index()])

    def get_vector(self):
        return np.array([self.__node2.get_node_x() - self.__node1.get_node_x(), self.__node2.get_node_y() - self.__node1.get_node_y()])

    def get_length(self):
        vector = self.get_vector()
        return np.linalg.norm(vector)
    
    def get_index(self):
        return self.__index
    
    def get_node1(self):
        return self.__node1
    
    def get_node2(self):
        return self.__node2

    def __repr__(self):
        self.update()
        return "Edge: id={1} ( nodes_id: {0}, Mesh2_edge_x: {2}, Mesh2_edge_y: {3})".format(
                self.get_edge_nodes(), self.__index, self.get_edge_x(), self.get_edge_y())


class Mesh2_face(object):
    ''' Class describing a Face (Polygon) - Triangle or Rectangular (!!!) shape for 2D unstructured grid '''
    def __init__(self, nodes, index):
        # ------------------
        # <nodes>  - is a list of objects of type <Mesh2_node>
        #            nodes should be ordered and go counter clockwise
        # <index>  - integer indicating id of current element
        self.__nodes = nodes
        self.__index = index
        
        #attrs = vars(self)
        #print '\n'.join("%s: %s" % item for item in attrs.items())
        self.update()
    
    def update(self):
        self.__nNodes = len(self.__nodes)
        self.create_attributes()
        self.validateInput()
        self.calculateGeometry()
 
    def create_attributes(self):
        # check number of nodes 3 or 4
        if self.__nNodes == 4:
            self.__node1 = self.__nodes[0]
            self.__node2 = self.__nodes[1]
            self.__node3 = self.__nodes[2]
            self.__node4 = self.__nodes[3]
            self.__nodelist = [self.__node1, self.__node2, self.__node3, self.__node4]
            
            self.__edge1 = Mesh2_edge(self.__node1, self.__node2, -999)
            self.__edge2 = Mesh2_edge(self.__node2, self.__node3, -999)
            self.__edge3 = Mesh2_edge(self.__node3, self.__node4, -999)
            self.__edge4 = Mesh2_edge(self.__node4, self.__node1, -999)
            self.__edgelist = [self.__edge1, self.__edge2, self.__edge3, self.__edge4]
       
        elif self.__nNodes == 3:
            self.__node1 = self.__nodes[0]
            self.__node2 = self.__nodes[1]
            self.__node3 = self.__nodes[2]
            self.__nodelist = [self.__node1, self.__node2, self.__node3]
            
            self.__edge1 = Mesh2_edge(self.__node1, self.__node2, -999)
            self.__edge2 = Mesh2_edge(self.__node2, self.__node3, -999)
            self.__edge3 = Mesh2_edge(self.__node3, self.__node1, -999)
            self.__edgelist = [self.__edge1, self.__edge2, self.__edge3]
        else:
            err_msg = '<Mesh2_face> id={1}\npassed <{0}> nodes to <Mesh2_face>. Polygon must have 3 or 4 nodes'.format(
                        self.__nNodes, self.__index())
            raise ValueError(err_msg)
    
    def validateInput(self):
        # VALIDATE DATA
        # --------------------------------
        
        for item in self.__nodes:
            if not isinstance(item, Mesh2_node):
                err_msg = '<Mesh2_face> id={2}\nobject <{0}> of type <{1}> passed as input to <Mesh2_face> object must be of type <Mesh2_node>'.format(
                            item, type(item), self.get_index())
                raise TypeError(err_msg)
        
        # check if points are colinear
        matrix = np.array([edge.get_vector() for edge in self.__edgelist])
        rank_matrix = np.linalg.matrix_rank(matrix)
        if rank_matrix == 1:  # if vectors are not colinear >>> rank > 1
            err_msg = '<Mesh2_face> id={2}\nNodes ...\n\t{0}\n...passed as input to <Mesh2_face> object are laying on one line (forming colinear vectors)'.format(
                            '\n\t'.join([str(n) for n in self.__nodelist]), type(item), self.get_index())
            raise ValueError(err_msg)
        
        # check if 4nodes makes a rectangular shape
        if self.__nNodes == 4:
            v1 = self.__edge1
            v2 = self.__edge2
            v3 = self.__edge3
            v4 = self.__edge4
            for edge1, edge2 in zip([v1, v2, v3, v4], [v4, v1, v2, v3]):
                product = np.dot(edge1.get_vector(), edge2.get_vector())  # if vectors are orthogonal dot product == 0
                error_margin = 1.e-5  # error margin for determining orthogonallity
                #if product != 0:
                if abs(product)-error_margin > 0:
                    err_msg = ('Mesh2_face: index={2}\n\tdot product of Edge1,Edge2 = {5} (if=0 => orthogonal)\n\t{0} >>> {3}\n\t{1} >>> {4}\nVectors (calculated in <Mesh2_face> object) formed from element nodes are not orthogonal.\nPolygon is not a rectangle'.format(edge1, edge2, self.get_index(), [edge1.get_node1(), edge1.get_node2()], [edge2.get_node1(), edge2.get_node2()], product ))
                    #raise ValueError(err_msg)
                    #print err_msg
    
    def polygonArea(self, X=None, Y=None):
        '''
        X, Y        Arrays of the x and y coordinates of the vertices, traced in a clockwise direction, starting at any vertex.
                    If you trace them counterclockwise, the result will be correct but have a negative sign.
        Returns     the area of the polygon
        '''
        # if user has not forced his own X,Y array work with those, we have...
        if X is None or Y is None:
            X = [n.get_node_x() for n in self.__nodelist]
            Y = [n.get_node_y() for n in self.__nodelist]


        area = 0.         # Accumulates area in the loop
        ox, oy = X[0], Y[0]
        for x, y in zip(X[1:], Y[1:]):
            area  += (x*oy-y*ox)
            ox, oy = x, y
        return abs(area*0.5)
    
    def polygonCentroid(self, X=None, Y=None):
        '''
        X, Y        Arrays of the x and y coordinates of the vertices, traced in a clockwise direction, starting at any vertex.
                    If you trace them counterclockwise, the result will be correct but have a negative sign.
        Returns     the centroid coords of the polygon [x, y]
        '''
        # if user has not forced his own X,Y array work with those, we have...
        if X is None or Y is None:
            X = [n.get_node_x() for n in self.__nodelist]
            Y = [n.get_node_y() for n in self.__nodelist]
        sum_x = 0.
        sum_y = 0.
        for x, y in zip(X, Y):
            sum_x += x
            sum_y += y
        return [sum_x/float(len(X)), sum_y/float(len(Y))]

    def calculateGeometry(self):
        # for rectangular select nodes 1, 2 and 4 , which are create orthogonal vectors
        if self.__nNodes == 3:
            n1, n2, n3 = self.__node1, self.__node2, self.__node3
        elif self.__nNodes == 4:
            n1, n2, n3 = self.__node1, self.__node2, self.__node4
        
        # calculating polygon area
        # ------------------------
        self.__area = self.polygonArea()
        
        # calculating polygon centroid
        # ----------------------------
        centroid_center = self.polygonCentroid()

        self.__face_x = centroid_center[0]
        self.__face_y = centroid_center[1]

        # calculating polygon circumscribed circle center
        # -----------------------------------------------
        if self.__nNodes == 4:
            self.__face_center_x = centroid_center[0]
            self.__face_center_y = centroid_center[1]
        
        elif self.__nNodes == 3:  # if triangle
            # equation from   http://en.wikipedia.org/wiki/Circumscribed_circle
            
            D = 2*np.linalg.norm(np.cross((n1.get_vector()-n2.get_vector()), (n2.get_vector()-n3.get_vector())))**2
            #alpha
            #np.linalg.norm(vector2)  - return modulus of vector
            alpha = np.linalg.norm(n2.get_vector()-n3.get_vector())**2*np.dot((n1.get_vector()-n2.get_vector()), (n1.get_vector()-n3.get_vector()))/D
            #betta
            betta = np.linalg.norm(n1.get_vector()-n3.get_vector())**2*np.dot((n2.get_vector()-n1.get_vector()), (n2.get_vector()-n3.get_vector()))/D
            #gamma
            gamma = np.linalg.norm(n1.get_vector()-n2.get_vector())**2*np.dot((n3.get_vector()-n1.get_vector()), (n3.get_vector()-n2.get_vector()))/D

            circumscribed_center = alpha*n1.get_vector()+betta*n2.get_vector()+gamma*n3.get_vector()
            self.__face_center_x = circumscribed_center[0]
            self.__face_center_y = circumscribed_center[1]
    
    def set_index(self, index):
        self.__index = int(index)

    def set_face_x(self, x):
        self.__face_x = x
    
    def set_face_y(self, y):
        self.__face_y = y

    def set_face_center_x(self, x):
        self.__face_center_x = x
    
    def set_face_center_y(self, y):
        self.__face_center_y = y

    def set_edges_indexes(self, index_list):
        if len(index_list) != self.__nNodes:
            print self
            raise ValueError('number of indexes does not match number of edges')
        i = 0
        for edge, index in zip(self.__edgelist, index_list):
            edge.set_index(index)
            self.__edgelist[i] = edge
            i += 1
    

    # >>> getters >>>

    def get_face_x(self):
        #self.update()
        return self.__face_x
    
    def get_face_y(self):
        #self.update()
        return self.__face_y

    def get_face_center_x(self):
        #self.update()
        return self.__face_center_x

    def get_face_center_y(self):
        #self.update()
        return self.__face_center_y

    def get_area(self):
        #self.update()
        return self.__area
    
    def get_index(self):
        return self.__index
    
    def get_face_edges(self):
        #self.update()
        return np.array([e.get_index() for e in self.__edgelist])
    
    def get_face_nodes(self):
        #self.update()
        return np.array([n.get_index() for n in self.__nodelist])

    def get_nNodes(self):
        #self.update()
        return self.__nNodes
    
    def get_nodes(self):
        #self.update()
        return self.__nodelist
    
    def get_edges(self):
        #self.update()
        return self.__edgelist
    
    def get_edges_indexes(self):
        #self.update()
        l = []
        for edge in self.__edgelist:
            l.append(edge.get_index())
        return l


    def __repr__(self):
        #self.update()
        return "Polygon: id={5} ( nodes_id: {0}, edges_id: {6} Mesh2_face_x: {1}, Mesh2_face_y: {2}, Mesh2_face_center_x: {3}, Mesh2_face_center_y: {4} )".format(
                [n.get_index() for n in self.__nodelist], self.__face_x, self.__face_y, self.__face_center_x, self.__face_center_y, self.__index, self.get_edges_indexes())

#
#
#
# --------------------------------------------------------------
# -------------------------- GRID ------------------------------
# --------------------------------------------------------------
#
#


class Grid2D(object):  #Mesh2_node, Mesh2_edge, Mesh2_face, object):

    def __init__(self, x_coord, y_coord, mask=None, data_location='T_points'):
        '''
        -----------------------------------------------
        Function creates unstructured grid of size.
        -----------------------------------------------
        Currently 4 input possibilities are accepted:

                1) Rectengular uniform grid. (all cells have equal size), data_location='T_points'
                    inputs:
                        x_coord - 1d array of length nx
                        y_coord - 1d array of length ny
                    Bathymetry (at T_points) is specified on grid (ny, nx)

                2) Rectengular uniform grid. (all cells have equal size), data_location='X_points'
                    inputs:
                        x_coord - 1d array of length nx
                        y_coord - 1d array of length ny
                    Bathymetry (at T_points) is specified on grid (nx-1, ny-1)
         
                3) Rectengular grid. (cells may have non-qual size), data_location='X_points'
                    Bathymetry (at T_points) is specified on grid (ny-1, nx-1)
                
                4) Curvilinear grid. data_location='X_points'
                    inputs:
                        x_coord - 2d array of shape (ny, nx)
                        y_coord - 2d array of shape (ny, nx)
                    Bathymetry (at T_points) is specified on grid (ny-1, nx-1)

        -----------------------------------------------
        Inputs:
        -----------------------------------------------
        x_coord - numpy array, x-coordinates of T or X points. Can be 1d or 2d
        y_coord - numpy array, y-coordinates of T or X points. Can be 1d or 2d


        mask - is an 2D array of shape
                - (len(x_coord), len(y_coord)), if x and y are 1D
                - shape(x_coord)=shape(y_coord), if x and y are 2D
            filled with True, False (True - for masked cell, False - for non-masked)
        
        data_location - string, 'T_points' or 'X_points', defines where do the arrays x_coord and y_coord show the values
        '''
        self._x_coord = x_coord
        self._y_coord = y_coord
        

        # ------------------------------------
        # --             CASE 1             --
        # --    rectangular uniform grid    --
        # --        data at T_points        --
        # ------------------------------------
        if (len(x_coord.shape) == 1) and (len(y_coord.shape) == 1) and (data_location == 'T_points'):
            self._nx = len(x_coord)
            self._ny = len(y_coord)
            
            # we cannot calculate the node coordinates in this case if the cells are not uniform!
            # therefore here WE ASSUME that all the cells are uniform!
            
            # now let us check the dictance between cells
            tolerance = 0.01  # tolerance in [units] to check grid uniformity
            
            for i, x in enumerate(x_coord[2:]):
                if abs(abs(x_coord[i]-x_coord[i+1])-abs(x_coord[i+1]-x)) > tolerance:
                    errmsg = 'Coordinates to generate grid are passed as two 1D arrays at T_points. Grid must be uniform in this case. Therefore the distance between grid cells should be equal:\n\tabs( abs(x[i]-x[i+1])-abs(x[i+1]-x[i+2]) ) < tolerance\nBut the passed x_coordinate array does not comply:\n\tabs( abs(x[{0}]-x[{1}])-abs(x[{1}]-x[{2}]) ) > tolerance\n\tabs( abs({3} - {4}) - abs({4} - {5}) ) > {6}\nYou can try following to fix the problem:\n\t a) increase tolerance\n\t b) pass coordinates of X_points'.format(i, i+1, i+2, x_coord[i], x_coord[i+1], x, tolerance)
                    raise ValueError(errmsg)
            
            for i, y in enumerate(y_coord[2:]):
                if abs(abs(y_coord[i]-y_coord[i+1])-abs(y_coord[i+1]-y)) > tolerance:
                    errmsg = 'Coordinates to generate grid are passed as two 1D arrays at T_points. Grid must be uniform in this case. Therefore the distance between grid cells should be equal:\n\tabs( abs(x[i]-x[i+1])-abs(x[i+1]-x[i+2]) ) < tolerance\nBut the passed x_coordinate array does not comply:\n\tabs( abs(x[{0}]-x[{1}])-abs(x[{1}]-x[{2}]) ) > tolerance\n\tabs( abs({3} - {4}) - abs({4} - {5}) ) > {6}\nYou can try following to fix the problem:\n\t a) increase tolerance\n\t b) pass coordinates of X_points'.format(i, i+1, i+2, x_coord[i], x_coord[i+1], x, tolerance)
                    raise ValueError(errmsg)
            
            # Since cells are uniform it doesnot make sence to find <dx> and <dy> for each cell within the loop,
            # and we can simply do it one time
            self._dx = abs(x_coord[0]-x_coord[1])/2.  # x- projection of distance from cell center to its left/right border
            self._dy = abs(y_coord[0]-y_coord[1])/2.  # y- projection of distance from cell center to its top/bottom border
        
        # ------------------------------------
        # --             CASE 2             --
        # --        rectangular grid        --
        # --        data at X_points        --
        # ------------------------------------
        elif (len(x_coord.shape) == 1) and (len(y_coord.shape) == 1) and (data_location == 'X_points'):
            self._nx = len(x_coord)-1
            self._ny = len(y_coord)-1


        # ------------------------------------
        # --             CASE 3             --
        # --        curvilinear grid        --
        # --        data at X_points        --
        # ------------------------------------
        elif (len(x_coord.shape) == 2) and (len(y_coord.shape) == 2) and data_location == 'X_points':
            self._nx = x_coord.shape[1]-1
            self._ny = x_coord.shape[0]-1
        
        else:
            errmsg = 'INVALID GRID TYPE. '+gridhelp()+'\nPassed arrays are:\n\t x_coord shape: ({0})\n\t y_coord shape: ({1})\n\t data location: {2}'.format(
                     x_coord.shape, y_coord.shape, data_location)
            raise ValueError(errmsg)
        

        self._faceMap = np.zeros((self._ny, self._nx))
        self._faceMap[:] = -999  # -999 is a missingValue indicator
    
        if mask is not None:
            if self._faceMap.shape != np.array(mask).shape:
                print 'faceMap shape :{0}\nmask shape: {1}'.format(self._faceMap.shape, np.array(mask).shape)
                raise ValueError('Masked array has invalid shape')
            else:
                self._faceMap = ma.array(self._faceMap, mask=mask)
        self._nNodes  = 0
        self._nEdges  = 0
        self._nFaces  = 0
        self._nodes   = []
        self._edges   = []
        self._faces   = []
        self._nMaxMesh2_face_nodes = 4  # HARDCODE for mossco ( rectangular elements)
        self.update_mossco_grid(data_location=data_location)


    

    def update_mossco_grid(self, data_location='T_points'):
        # coordinates (0,0) is located at top left
        # y axis from top left to bottom left
        # x axis from top left to top right
        nNodes = 0
        nEdges = 0
        nFaces = 0
        # cycle through all valid cells
        for j in xrange(self._ny):
            for i in xrange(self._nx):
                if self._faceMap[j, i]:
                    #a valid cell center is found
                    self._nFaces += 1
                    self._faceMap[j, i] = int(self._nFaces)

                    # now get coordinates of the face center and nodes
                    # case 1
                    if (len(self._x_coord.shape) == 1) and (len(self._y_coord.shape) == 1) and (data_location == 'T_points'):
                        x, y = self._x_coord[i], self._y_coord[j]  # x,y center coord
                        lt_x, lt_y = x-self._dx, y-self._dy  # left top node
                        lb_x, lb_y = x-self._dx, y+self._dy  # left bottom node
                        rb_x, rb_y = x+self._dx, y+self._dy  # right bottom node
                        rt_x, rt_y = x+self._dx, y-self._dy  # right top node
                    
                    # case 2
                    elif (len(self._x_coord.shape) == 1) and (len(self._y_coord.shape) == 1) and (data_location == 'X_points'):
                        lt_x, lt_y = self._x_coord[i  ], self._y_coord[j  ]  # left top node
                        lb_x, lb_y = self._x_coord[i  ], self._y_coord[j+1]  # left bottom node
                        rb_x, rb_y = self._x_coord[i+1], self._y_coord[j+1]  # right bottom node
                        rt_x, rt_y = self._x_coord[i+1], self._y_coord[j  ]  # right top node


                    # case 3
                    elif (len(self._x_coord.shape) == 2) and (len(self._y_coord.shape) == 2) and data_location == 'X_points':
                        #x, y = self._x_coord[i], self._y_coord[j]  # x,y center coord
                        lt_x, lt_y = self._x_coord[j  , i  ], self._y_coord[j,   i  ]  # left top node
                        lb_x, lb_y = self._x_coord[j+1, i  ], self._y_coord[j+1, i  ]  # left bottom node
                        rb_x, rb_y = self._x_coord[j+1, i+1], self._y_coord[j+1, i+1]  # right bottom node
                        rt_x, rt_y = self._x_coord[j  , i+1], self._y_coord[j  , i+1]  # right top node


                    tr = False  # flag indicating, that there is a neighbour cell on top-right (diagonally)
                    t  = False  # flag indicating, that there is a neighbour cell on top
                    tl = False  # flag indicating, that there is a neighbour cell on top-left (diagonally)
                    l  = False  # flag indicating, that there is a neighbour cell on left

                    # now check if there are aleready existing polygons 'left', or 'top'
                    #('bottom' and 'right' are not possible due to i,j cycling mechanism)
                    #
                    # also take care of first vertical and horizontal lines i=0 and j=0
                    # UPDATE: need more control. should be remade ito zones

                    if 0 < i < self._nx-1 and 0 < j:
                        # BODY
                        faceIndex = int(self._nFaces)
                        nNodes = 4
                        nEdges = 4
                        #
                        if self._faceMap[j-1, i+1]:
                            tr = True
                            tr_face = self.get_face_by_index(self._faceMap[j-1, i+1])
                        if self._faceMap[j-1, i]:
                            t = True
                            t_face = self.get_face_by_index(self._faceMap[j-1, i])
                        if self._faceMap[j-1, i-1]:
                            tl = True
                            tl_face = self.get_face_by_index(self._faceMap[j-1, i-1])
                        if self._faceMap[j, i-1]:
                            l = True
                            l_face = self.get_face_by_index(self._faceMap[j, i-1])

                        if tr:
                            nNodes += -1
                            node_rt = tr_face.get_nodes()[1]  # right top node
                        if tl:
                            nNodes += -1
                            node_lt = tl_face.get_nodes()[2]  # left top node
                        
                        if t and not tr and not tl:
                            nNodes += -2
                            nEdges += -1
                            node_rt = t_face.get_nodes()[2]  # right top node
                            node_lt = t_face.get_nodes()[1]  # left top node
                            edge_t = t_face.get_edges()[1]  #bottom edge of the top cell is the top edge of current cell
                        elif t and tr and not tl:
                            nNodes += -1
                            nEdges += -1
                            node_lt = t_face.get_nodes()[1]  # left top node
                            edge_t = t_face.get_edges()[1]  #bottom edge of the top cell is the top edge of current cell
                        elif t and not tr and tl:
                            nNodes += -1
                            nEdges += -1
                            node_rt = t_face.get_nodes()[2]  # right top node
                            edge_t = t_face.get_edges()[1]  #bottom edge of the top cell is the top edge of current cell
                        elif t and tr and tl:
                            nEdges += -1
                            edge_t = t_face.get_edges()[1]  #bottom edge of the top cell is the top edge of current cell

                        if l and not (tl or t):
                            nNodes += -2
                            nEdges += -1
                            edge_l = l_face.get_edges()[2]  #right edge of the left cell is the left edge of current cell
                            node_lt = l_face.get_nodes()[3]  # left top node
                            node_lb = l_face.get_nodes()[2]  # left-bottom node
                        elif l and (tl or t):
                            nNodes += -1
                            nEdges += -1
                            edge_l = l_face.get_edges()[2]  #right edge of the left cell is the left edge of current cell
                            node_lb = l_face.get_nodes()[2]  # left-bottom node

                        self._nNodes += nNodes
                        self._nEdges += nEdges


                        #create nodes and append them
                        if not (t or tl or l):
                            node_lt = Mesh2_node(lt_x, lt_y, self._nNodes-nNodes+1)  # left top node
                            self._nodes.append(node_lt)
                            nNodes += -1
                        if not l:
                            node_lb = Mesh2_node(lb_x, lb_y, self._nNodes-nNodes+1)  # left bottom node
                            self._nodes.append(node_lb)
                            nNodes += -1
                        node_rb = Mesh2_node(rb_x, rb_y, self._nNodes-nNodes+1)  # right bottom node
                        self._nodes.append(node_rb)
                        nNodes += -1
                        if not (tr or t):
                            node_rt = Mesh2_node(rt_x, rt_y, self._nNodes-nNodes+1)  # right top node
                            self._nodes.append(node_rt)
                            nNodes += -1
                        
                        # create FACE
                        face = Mesh2_face([node_lt, node_lb, node_rb, node_rt], faceIndex)
                        # reIndex its edges and append them
                        if t and not l:
                            face.set_edges_indexes([self._nEdges-2, self._nEdges-1, self._nEdges, edge_t.get_index()])
                            edge_l, edge_b, edge_r, edge_t = face.get_edges()
                            for edge in [edge_l, edge_b, edge_r]:
                                self._edges.append(edge)
                        elif t and l:
                            face.set_edges_indexes([edge_l.get_index(), self._nEdges-1, self._nEdges, edge_t.get_index()])
                            edge_l, edge_b, edge_r, edge_t = face.get_edges()
                            for edge in [edge_b, edge_r]:
                                self._edges.append(edge)
                        elif not t and l:
                            face.set_edges_indexes([edge_l.get_index(), self._nEdges-2, self._nEdges-1, self._nEdges])
                            edge_l, edge_b, edge_r, edge_t = face.get_edges()
                            for edge in [edge_b, edge_r, edge_t]:
                                self._edges.append(edge)
                        elif not t and not l:
                            face.set_edges_indexes([self._nEdges-3, self._nEdges-2, self._nEdges-1, self._nEdges])
                            edge_l, edge_b, edge_r, edge_t = face.get_edges()
                            for edge in [edge_l, edge_b, edge_r, edge_t]:
                                self._edges.append(edge)
                        # now, after edges were changed in face, we can append it
                        self._faces.append(face)
                        #
                        #
                        #
                        #
                        #
                        #
                        #
                    elif i == 0 and j == 0:
                        #print "found cell on the border: i=0, j=0. Not implemented yet. Border cells should be empty"
                        faceIndex = int(self._nFaces)
                        nNodes = 4
                        nEdges = 4
                        self._nNodes += nNodes
                        self._nEdges += nEdges
                        node_lt = Mesh2_node(lt_x, lt_y, self._nNodes-nNodes+1)  # left top node
                        nNodes += -1
                        node_lb = Mesh2_node(lb_x, lb_y, self._nNodes-nNodes+1)  # left bottom node
                        nNodes += -1
                        node_rb = Mesh2_node(rb_x, rb_y, self._nNodes-nNodes+1)  # right bottom node
                        nNodes += -1
                        node_rt = Mesh2_node(rt_x, rt_y, self._nNodes-nNodes+1)  # right top node
                        for node in [node_lt, node_lb, node_rb, node_rt]:
                            self._nodes.append(node)
                        face = Mesh2_face([node_lt, node_lb, node_rb, node_rt], faceIndex)
                        face.set_edges_indexes([self._nEdges-3, self._nEdges-2, self._nEdges-1, self._nEdges])
                        edge_l, edge_b, edge_r, edge_t = face.get_edges()
                        for edge in [edge_l, edge_b, edge_r, edge_t]:
                            self._edges.append(edge)
                        self._faces.append(face)
                        #
                        #
                        #
                        #
                        #
                        #
                        #
                    elif i == 0 and 0 < j :
                        #print "found cell on the border: i=0, 0<j. Not implemented yet. Border cells should be empty"

                        faceIndex = int(self._nFaces)
                        nNodes = 4  # to be placed... at this moment
                        nEdges = 4  # to be placed... at this moment
                        #
                        if self._faceMap[j-1, i+1]:
                            tr = True
                            tr_face = self.get_face_by_index(self._faceMap[j-1, i+1])
                        if self._faceMap[j-1, i]:
                            t = True
                            t_face = self.get_face_by_index(self._faceMap[j-1, i])
    

                        if tr:
                            nNodes += -1
                            node_rt = tr_face.get_nodes()[1]  # right top node
                        
                        if t and not tr:
                            nNodes += -2
                            nEdges += -1
                            node_rt = t_face.get_nodes()[2]  # right top node
                            node_lt = t_face.get_nodes()[1]  # left top node
                            edge_t = t_face.get_edges()[1]  #bottom edge of the top cell is the top edge of current cell
                        elif t and tr:
                            nNodes += -1
                            nEdges += -1
                            node_lt = t_face.get_nodes()[1]  # left top node
                            edge_t = t_face.get_edges()[1]  #bottom edge of the top cell is the top edge of current cell
                        
                        self._nNodes += nNodes
                        self._nEdges += nEdges

                        #create nodes and append them
                        if not t:
                            node_lt = Mesh2_node(lt_x, lt_y, self._nNodes-nNodes+1)  # left top node
                            self._nodes.append(node_lt)
                            nNodes += -1
                        node_lb = Mesh2_node(lb_x, lb_y, self._nNodes-nNodes+1)  # left bottom node
                        nNodes += -1
                        self._nodes.append(node_lb)
                        node_rb = Mesh2_node(rb_x, rb_y, self._nNodes-nNodes+1)  # right bottom node
                        self._nodes.append(node_rb)
                        nNodes += -1
                        if not tr:
                            node_rt = Mesh2_node(rt_x, rt_y, self._nNodes-nNodes+1)  # right top node
                            self._nodes.append(node_rt)
                            nNodes += -1
                        
                        # create FACE
                        face = Mesh2_face([node_lt, node_lb, node_rb, node_rt], faceIndex)
                        # reIndex its edges and append them
                        if t:
                            face.set_edges_indexes([self._nEdges-2, self._nEdges-1, self._nEdges, edge_t.get_index()])
                            edge_l, edge_b, edge_r, edge_t = face.get_edges()
                            for edge in [edge_l, edge_b, edge_r]:
                                self._edges.append(edge)
                        elif not t:
                            face.set_edges_indexes([self._nEdges-3, self._nEdges-2, self._nEdges-1, self._nEdges])
                            edge_l, edge_b, edge_r, edge_t = face.get_edges()
                            for edge in [edge_l, edge_b, edge_r, edge_t]:
                                self._edges.append(edge)
                        # now, after edges were changed in face, we can append it
                        self._faces.append(face)
                        #
                        #
                        #
                        #
                        #
                        #
                        #
                        #
                    elif 0 < i and j == 0:
                        #print "found cell on the border:  0 < i, j = 0. Not implemented yet. Border cells should be empty"
                        
                        faceIndex = int(self._nFaces)
                        nNodes = 4
                        nEdges = 4
                        self._nNodes += nNodes
                        self._nEdges += nEdges

                        if self._faceMap[j, i-1]:
                            l = True
                            l_face = self.get_face_by_index(self._faceMap[j, i-1])
                        if l:
                            nNodes += -2
                            nEdges += -1
                            edge_l = l_face.get_edges()[2]  #right edge of the left cell is the left edge of current cell
                            node_lt = l_face.get_nodes()[3]  # left top node
                            node_lb = l_face.get_nodes()[2]  # left-bottom node

                        self._nNodes += nNodes
                        self._nEdges += nEdges


                        #create nodes and append them
                        if not l:
                            node_lt = Mesh2_node(lt_x, lt_y, self._nNodes-nNodes+1)  # left top node
                            nNodes += -1
                            self._nodes.append(node_lt)
                            node_lb = Mesh2_node(lb_x, lb_y, self._nNodes-nNodes+1)  # left bottom node
                            nNodes += -1
                            self._nodes.append(node_lb)

                        node_rb = Mesh2_node(rb_x, rb_y, self._nNodes-nNodes+1)  # right bottom node
                        nNodes += -1
                        node_rt = Mesh2_node(rt_x, rt_y, self._nNodes-nNodes+1)  # right top node
                        
                        for node in [node_rb, node_rt]:
                            self._nodes.append(node)
                        face = Mesh2_face([node_lt, node_lb, node_rb, node_rt], faceIndex)
                        if not l:
                            face.set_edges_indexes([self._nEdges-3, self._nEdges-2, self._nEdges-1, self._nEdges])
                            edge_l, edge_b, edge_r, edge_t = face.get_edges()
                            for edge in [edge_l, edge_b, edge_r, edge_t]:
                                self._edges.append(edge)
                        elif l:
                            face.set_edges_indexes([edge_l.get_index(), self._nEdges-2, self._nEdges-1, self._nEdges])
                            edge_l, edge_b, edge_r, edge_t = face.get_edges()
                            for edge in [edge_b, edge_r, edge_t]:
                                self._edges.append(edge)
                        self._faces.append(face)
                        #
                        #
                        #
                        #
                        #
                        #
                        #
                        #
                    elif i == self._nx-1 and 0 < j :
                        #print "found cell on the border:  i =nx and 0 < j. Not implemented yet. Border cells should be empty"
                        # BODY
                        faceIndex = int(self._nFaces)
                        nNodes = 4
                        nEdges = 4

                        if self._faceMap[j-1, i]:
                            t = True
                            t_face = self.get_face_by_index(self._faceMap[j-1, i])
                        if self._faceMap[j-1, i-1]:
                            tl = True
                            tl_face = self.get_face_by_index(self._faceMap[j-1, i-1])
                        if self._faceMap[j, i-1]:
                            l = True
                            l_face = self.get_face_by_index(self._faceMap[j, i-1])

                        if tl:
                            nNodes += -1
                            node_lt = tl_face.get_nodes()[2]  # left top node
                        
                        if t and not tl:
                            nNodes += -2
                            nEdges += -1
                            node_rt = t_face.get_nodes()[2]  # right top node
                            node_lt = t_face.get_nodes()[1]  # left top node
                            edge_t = t_face.get_edges()[1]  #bottom edge of the top cell is the top edge of current cell

                        elif t and tl:
                            nNodes += -1
                            nEdges += -1
                            node_rt = t_face.get_nodes()[2]  # right top node
                            edge_t = t_face.get_edges()[1]  #bottom edge of the top cell is the top edge of current cell

                        if l and not (tl or t):
                            nNodes += -2
                            nEdges += -1
                            edge_l = l_face.get_edges()[2]  #right edge of the left cell is the left edge of current cell
                            node_lt = l_face.get_nodes()[3]  # left top node
                            node_lb = l_face.get_nodes()[2]  # left-bottom node
                        elif l and (tl or t):
                            nNodes += -1
                            nEdges += -1
                            edge_l = l_face.get_edges()[2]  #right edge of the left cell is the left edge of current cell
                            node_lb = l_face.get_nodes()[2]  # left-bottom node

                        self._nNodes += nNodes
                        self._nEdges += nEdges

                        #create nodes and append them
                        if not (t or tl or l):
                            node_lt = Mesh2_node(lt_x, lt_y, self._nNodes-nNodes+1)  # left top node
                            self._nodes.append(node_lt)
                            nNodes += -1
                        if not l:
                            node_lb = Mesh2_node(lb_x, lb_y, self._nNodes-nNodes+1)  # left bottom node
                            self._nodes.append(node_lb)
                            nNodes += -1
                        
                        node_rb = Mesh2_node(rb_x, rb_y, self._nNodes-nNodes+1)  # right bottom node
                        self._nodes.append(node_rb)
                        nNodes += -1
                        if not t:
                            node_rt = Mesh2_node(rt_x, rt_y, self._nNodes-nNodes+1)  # right top node
                            self._nodes.append(node_rt)
                            nNodes += -1
                        
                        # create FACE
                        face = Mesh2_face([node_lt, node_lb, node_rb, node_rt], faceIndex)
                        # reIndex its edges and append them
                        if t and not l:
                            face.set_edges_indexes([self._nEdges-2, self._nEdges-1, self._nEdges, edge_t.get_index()])
                            edge_l, edge_b, edge_r, edge_t = face.get_edges()
                            for edge in [edge_l, edge_b, edge_r]:
                                self._edges.append(edge)
                        elif t and l:
                            face.set_edges_indexes([edge_l.get_index(), self._nEdges-1, self._nEdges, edge_t.get_index()])
                            edge_l, edge_b, edge_r, edge_t = face.get_edges()
                            for edge in [edge_b, edge_r]:
                                self._edges.append(edge)
                        elif not t and l:
                            face.set_edges_indexes([edge_l.get_index(), self._nEdges-2, self._nEdges-1, self._nEdges])
                            edge_l, edge_b, edge_r, edge_t = face.get_edges()
                            for edge in [edge_b, edge_r, edge_t]:
                                self._edges.append(edge)
                        elif not t and not l:
                            face.set_edges_indexes([self._nEdges-3, self._nEdges-2, self._nEdges-1, self._nEdges])
                            edge_l, edge_b, edge_r, edge_t = face.get_edges()
                            for edge in [edge_l, edge_b, edge_r, edge_t]:
                                self._edges.append(edge)
                        # now, after edges were changed in face, we can append it
                        self._faces.append(face)



    def get_face_by_index(self, index):
        """
        # list self._faces contains faces in ascending order:
        # NOTE: object indixes start here from <1> and not from <0>
        # self._faces[0].get_index() >>> 1
        # self._faces[1].get_index() >>> 2
        # self._faces[2].get_index() >>> 3 etc.
        """
        #print "\t\t\tgetting face by index:"
        #print '\t\t\tcurrently ', len(self._faces), ' faces'
        return self._faces[int(index)-1]


    def update(self):
        pass
    
    def get_faces(self):
        return self._faces

    def get_nodes(self):
        return self._nodes
    

    def get_edges(self):
        return self._edges
    #
    #
    # ------- > What we need....
    #
    #

    def get_Mesh2_edge_nodes(self):

        a = np.zeros((self._nEdges, 2))
        for i, e in enumerate(self._edges):
            a[i, 0], a[i, 1] = e.get_node1().get_index(), e.get_node2().get_index()
        return a


    def get_Mesh2_edge_bnd(self):
        x = np.zeros((self._nEdges, 2))
        y = np.zeros((self._nEdges, 2))
        for i, e in enumerate(self._edges):
            x[i, 0], x[i, 1] = e.get_node1().get_node_x(), e.get_node2().get_node_x()
            y[i, 0], y[i, 1] = e.get_node1().get_node_y(), e.get_node2().get_node_y()
        return x, y


    def get_Mesh2_edge_faces(self):
        # initializing
        self._Mesh2_edge_faces = np.zeros((self._nEdges, 2))
        self._Mesh2_edge_faces[:] = -1

        for j in xrange(self._ny):
            for i in xrange(self._nx):
                if self._faceMap[j, i]:
                    #print 'found_valid polygon'
                    # if at j,i there is a valid cell
                    #
                    #
                    #
                    # because we iterate over the cells directly the same way, that initialization was run,
                    # we know that if cell has a neighbour at TOP or LEFT it uses edges, generated with those neighbour cells,
                    # this fact should be taken into account when calculating connectivity (left/right side of the edges)
                    # ... and it is used :)
                    #
                    # here the index of an enge correspondes to its position in list, thus we can use
                    #    >>> self._Mesh2_edge_faces[ind_e_t-1, 0]

                    face = self.get_face_by_index(self._faceMap[j, i])
                    ind_e_l, ind_e_b, ind_e_r, ind_e_t = face.get_edges_indexes()

                    if 0 < i < self._nx-1 and 0 < j < self._ny-1:
                        #print '__ in the middle'
                        #print '\t it has edges:', ind_e_l, ind_e_b, ind_e_r, ind_e_t
                        # MAIN BODY

                        # if top has neighbour
                        if self._faceMap[ j-1, i]:
                            #print "\t polygon has top neigbour"
                            pass  # already included

                        else:
                            #print "\t polygon has no top neigbour"
                            self._Mesh2_edge_faces[ind_e_t-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_t-1, 1] = -999
                        
                        # if left has neighbour
                        if self._faceMap[ j, i-1]:
                            #print "\t polygon has left neigbour"
                            pass  # already included
                        else:
                            #print "\t polygon has no left neigbour"
                            self._Mesh2_edge_faces[ind_e_l-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_l-1, 1] = -999

                        # if bottom has neighbour
                        if self._faceMap[ j+1, i]:
                            #print "\t polygon has bottom neigbour"
                            self._Mesh2_edge_faces[ind_e_b-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_b-1, 1] = self.get_face_by_index(self._faceMap[ j+1, i]).get_index()
                        else:
                            #print "\t polygon has no bottom neigbour"
                            self._Mesh2_edge_faces[ind_e_b-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_b-1, 1] = -999

                        # if right has neighbour
                        if self._faceMap[ j, i+1]:
                            #print "\t polygon has right neigbour"
                            self._Mesh2_edge_faces[ind_e_r-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_r-1, 1] = self.get_face_by_index(self._faceMap[ j, i+1]).get_index()
                        else:
                            #print "\t polygon has no right neigbour"
                            self._Mesh2_edge_faces[ind_e_r-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_r-1, 1] = -999

                    elif i == 0 and j == 0:
                        self._Mesh2_edge_faces[ind_e_t-1, 0] = face.get_index()
                        self._Mesh2_edge_faces[ind_e_t-1, 1] = -999
                        
                        self._Mesh2_edge_faces[ind_e_l-1, 0] = face.get_index()
                        self._Mesh2_edge_faces[ind_e_l-1, 1] = -999
                        
                        # if bottom has neighbour
                        if self._faceMap[ j+1, i]:
                            self._Mesh2_edge_faces[ind_e_b-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_b-1, 1] = self.get_face_by_index(self._faceMap[ j+1, i]).get_index()
                        else:
                            self._Mesh2_edge_faces[ind_e_b-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_b-1, 1] = -999

                        # if right has neighbour
                        if self._faceMap[ j, i+1]:
                            self._Mesh2_edge_faces[ind_e_r-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_r-1, 1] = self.get_face_by_index(self._faceMap[ j, i+1]).get_index()
                        else:
                            self._Mesh2_edge_faces[ind_e_r-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_r-1, 1] = -999
                    
                    elif i == 0 and 0 < j < self._ny-1:
                        self._Mesh2_edge_faces[ind_e_l-1, 0] = face.get_index()
                        self._Mesh2_edge_faces[ind_e_l-1, 1] = -999

                        # if bottom has neighbour
                        if self._faceMap[ j+1, i]:
                            self._Mesh2_edge_faces[ind_e_b-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_b-1, 1] = self.get_face_by_index(self._faceMap[ j+1, i]).get_index()
                        else:
                            self._Mesh2_edge_faces[ind_e_b-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_b-1, 1] = -999

                        # if right has neighbour
                        if self._faceMap[ j, i+1]:
                            self._Mesh2_edge_faces[ind_e_r-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_r-1, 1] = self.get_face_by_index(self._faceMap[ j, i+1]).get_index()
                        else:
                            self._Mesh2_edge_faces[ind_e_r-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_r-1, 1] = -999
                        # if top has neighbour
                        if self._faceMap[ j-1, i]:
                            pass  # already included
                        else:
                            self._Mesh2_edge_faces[ind_e_t-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_t-1, 1] = -999

                    elif i == 0 and j == self._ny-1:
                        self._Mesh2_edge_faces[ind_e_l-1, 0] = face.get_index()
                        self._Mesh2_edge_faces[ind_e_l-1, 1] = -999
                        
                        self._Mesh2_edge_faces[ind_e_b-1, 0] = face.get_index()
                        self._Mesh2_edge_faces[ind_e_b-1, 1] = -999
                        # if right has neighbour
                        if self._faceMap[ j, i+1]:
                            self._Mesh2_edge_faces[ind_e_r-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_r-1, 1] = self.get_face_by_index(self._faceMap[ j, i+1]).get_index()
                        else:
                            self._Mesh2_edge_faces[ind_e_r-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_r-1, 1] = -999
                        # if top has neighbour
                        if self._faceMap[ j-1, i]:
                            pass  # already included
                        else:
                            self._Mesh2_edge_faces[ind_e_t-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_t-1, 1] = -999
                    
                    elif 0 < i < self._nx-1 and j == 0:
                        self._Mesh2_edge_faces[ind_e_t-1, 0] = face.get_index()
                        self._Mesh2_edge_faces[ind_e_t-1, 1] = -999

                        # if left has neighbour
                        if self._faceMap[ j, i-1]:
                            pass  # already included
                        else:
                            self._Mesh2_edge_faces[ind_e_l-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_l-1, 1] = -999

                        # if bottom has neighbour
                        if self._faceMap[ j+1, i]:
                            self._Mesh2_edge_faces[ind_e_b-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_b-1, 1] = self.get_face_by_index(self._faceMap[ j+1, i]).get_index()
                        else:
                            self._Mesh2_edge_faces[ind_e_b-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_b-1, 1] = -999

                        # if right has neighbour
                        if self._faceMap[ j, i+1]:
                            self._Mesh2_edge_faces[ind_e_r-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_r-1, 1] = self.get_face_by_index(self._faceMap[ j, i+1]).get_index()
                        else:
                            self._Mesh2_edge_faces[ind_e_r-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_r-1, 1] = -999

                    elif 0 < i < self._nx-1 and j == self._ny-1:
                        self._Mesh2_edge_faces[ind_e_b-1, 0] = face.get_index()
                        self._Mesh2_edge_faces[ind_e_b-1, 1] = -999

                        # if top has neighbour
                        if self._faceMap[ j-1, i]:
                            pass  # already included
                        else:
                            self._Mesh2_edge_faces[ind_e_t-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_t-1, 1] = -999
                        
                        # if left has neighbour
                        if self._faceMap[ j, i-1]:
                            pass  # already included
                        else:
                            self._Mesh2_edge_faces[ind_e_l-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_l-1, 1] = -999

                        # if right has neighbour
                        if self._faceMap[ j, i+1]:
                            self._Mesh2_edge_faces[ind_e_r-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_r-1, 1] = self.get_face_by_index(self._faceMap[ j, i+1]).get_index()
                        else:
                            self._Mesh2_edge_faces[ind_e_r-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_r-1, 1] = -999
                    
                    elif i == self._nx-1 and j == 0:
                        self._Mesh2_edge_faces[ind_e_r-1, 0] = face.get_index()
                        self._Mesh2_edge_faces[ind_e_r-1, 1] = -999

                        self._Mesh2_edge_faces[ind_e_t-1, 0] = face.get_index()
                        self._Mesh2_edge_faces[ind_e_t-1, 1] = -999

                        # if left has neighbour
                        if self._faceMap[ j, i-1]:
                            pass  # already included
                        else:
                            self._Mesh2_edge_faces[ind_e_l-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_l-1, 1] = -999

                        # if bottom has neighbour
                        if self._faceMap[ j+1, i]:
                            self._Mesh2_edge_faces[ind_e_b-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_b-1, 1] = self.get_face_by_index(self._faceMap[ j+1, i]).get_index()
                        else:
                            self._Mesh2_edge_faces[ind_e_b-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_b-1, 1] = -999

                    elif i == self._nx-1 and 0 < j < self._ny-1:
                        self._Mesh2_edge_faces[ind_e_r-1, 0] = face.get_index()
                        self._Mesh2_edge_faces[ind_e_r-1, 1] = -999

                        # if top has neighbour
                        if self._faceMap[ j-1, i]:
                            pass  # already included
                        else:
                            self._Mesh2_edge_faces[ind_e_t-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_t-1, 1] = -999
                        
                        # if left has neighbour
                        if self._faceMap[ j, i-1]:
                            pass  # already included
                        else:
                            self._Mesh2_edge_faces[ind_e_l-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_l-1, 1] = -999

                        # if bottom has neighbour
                        if self._faceMap[ j+1, i]:
                            self._Mesh2_edge_faces[ind_e_b-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_b-1, 1] = self.get_face_by_index(self._faceMap[ j+1, i]).get_index()
                        else:
                            self._Mesh2_edge_faces[ind_e_b-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_b-1, 1] = -999

                    elif i == self._nx-1 and j == self._ny-1:
                        self._Mesh2_edge_faces[ind_e_r-1, 0] = face.get_index()
                        self._Mesh2_edge_faces[ind_e_r-1, 1] = -999
                    
                        self._Mesh2_edge_faces[ind_e_b-1, 0] = face.get_index()
                        self._Mesh2_edge_faces[ind_e_b-1, 1] = -999

                        # if top has neighbour
                        if self._faceMap[ j-1, i]:
                            pass  # already included
                        else:
                            self._Mesh2_edge_faces[ind_e_t-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_t-1, 1] = -999
                        
                        # if left has neighbour
                        if self._faceMap[ j, i-1]:
                            pass  # already included
                        else:
                            self._Mesh2_edge_faces[ind_e_l-1, 0] = face.get_index()
                            self._Mesh2_edge_faces[ind_e_l-1, 1] = -999

        if self._nEdges == len(self._Mesh2_edge_faces):
            return self._Mesh2_edge_faces
        else:
            err_msg = 'self._nEdges = {0}\nlen(self._Mesh2_edge_faces) = {1}'.format(self._nEdges, len(self._Mesh2_edge_faces))
            raise ValueError(err_msg)

    def get_Mesh2_face_nodes(self):
        a = np.zeros((self._nFaces, self._nMaxMesh2_face_nodes))
        a[:] = -999
        for i, f in enumerate(self._faces):
            if f.get_nNodes() == 3:
                a[i, 0], a[i, 1], a[i, 2] = f.get_nodes()[0].get_index(), f.get_nodes()[1].get_index(), f.get_nodes()[2].get_index()
            elif f.get_nNodes() == 4:
                a[i, 0], a[i, 1], a[i, 2], a[i, 3] = f.get_nodes()[0].get_index(), f.get_nodes()[1].get_index(), f.get_nodes()[2].get_index(), f.get_nodes()[3].get_index()
        return a
    
    
    def get_Mesh2_face_bnd(self):
        x = np.zeros((self._nFaces, self._nMaxMesh2_face_nodes))
        y = np.zeros((self._nFaces, self._nMaxMesh2_face_nodes))
        x[:] = -999
        y[:] = -999
        for i, f in enumerate(self._faces):
            if f.get_nNodes() == 3:
                x[i, 0], x[i, 1], x[i, 2] = f.get_nodes()[0].get_node_x(), f.get_nodes()[1].get_node_x(), f.get_nodes()[2].get_node_x()
                y[i, 0], y[i, 1], y[i, 2] = f.get_nodes()[0].get_node_y(), f.get_nodes()[1].get_node_y(), f.get_nodes()[2].get_node_y()
            elif f.get_nNodes() == 4:
                x[i, 0], x[i, 1], x[i, 2], x[i, 3] = f.get_nodes()[0].get_node_x(), f.get_nodes()[1].get_node_x(), f.get_nodes()[2].get_node_x(), f.get_nodes()[3].get_node_x()
                y[i, 0], y[i, 1], y[i, 2], y[i, 3] = f.get_nodes()[0].get_node_y(), f.get_nodes()[1].get_node_y(), f.get_nodes()[2].get_node_y(), f.get_nodes()[3].get_node_y()
        return x, y
    

    def get_Mesh2_face_edges(self):

        a = np.zeros((self._nFaces, self._nMaxMesh2_face_nodes))
        a[:] = -999
        for i, f in enumerate(self._faces):
            if f.get_nNodes() == 3:
                a[i, 0] = f.get_edges()[0].get_index()
                a[i, 1] = f.get_edges()[1].get_index()
                a[i, 2] = f.get_edges()[2].get_index()
            elif f.get_nNodes() == 4:
                a[i, 0] = f.get_edges()[0].get_index()
                a[i, 1] = f.get_edges()[1].get_index()
                a[i, 2] = f.get_edges()[2].get_index()
                a[i, 3] = f.get_edges()[3].get_index()
        return a

    def print_nodes_info(self):
        for n in self._nodes:
            print n

    def get_nMesh2_node(self):
        return self._nNodes

    def get_nMesh2_edge(self):
        return self._nEdges

    def get_nMesh2_face(self):
        return self._nFaces

    def get_nMaxMesh2_face_nodes(self):
        return self._nMaxMesh2_face_nodes
    
    def get_Mesh2_node_x(self):
        a = np.zeros(self._nNodes)
        for i, node in enumerate(self._nodes):
            a[i] = node.get_node_x()
        return a

    def get_Mesh2_node_y(self):
        a = np.zeros(self._nNodes)
        for i, node in enumerate(self._nodes):
            a[i] = node.get_node_y()
        return a

    def get_Mesh2_edge_x(self):
        a = np.zeros(self._nEdges)
        for i, edge in enumerate(self._edges):
            a[i] = edge.get_edge_x()
        return a
    
    def get_Mesh2_edge_y(self):
        a = np.zeros(self._nEdges)
        for i, edge in enumerate(self._edges):
            a[i] = edge.get_edge_y()
        return a

    def get_Mesh2_face_x(self):
        a = np.zeros(self._nFaces)
        for i, face in enumerate(self._faces):
            a[i] = face.get_face_x()
        return a

    def get_Mesh2_face_y(self):
        a = np.zeros(self._nFaces)
        for i, face in enumerate(self._faces):
            a[i] = face.get_face_y()
        return a
    
    def get_Mesh2_face_center_x(self):
        a = np.zeros(self._nFaces)
        for i, face in enumerate(self._faces):
            a[i] = face.get_face_center_x()
        return a
    
    def get_Mesh2_face_center_y(self):
        a = np.zeros(self._nFaces)
        for i, face in enumerate(self._faces):
            a[i] = face.get_face_center_y()
        return a

    def get_Mesh2_face_area(self):
        a = np.zeros(self._nFaces)
        for i, face in enumerate(self._faces):
            a[i] = face.get_area()
        return a


    def __repr__(self):
        #self.update()
        return "uGrid2D:\n\tnodes: {0}, \n\tedges: {1}, \n\tfaces: {2}".format(
                self._nNodes, self._nEdges, self._nFaces)


if __name__ == '__main__':
    def test1():
        x = np.linspace(5, 6, 6)
        y = np.linspace(3, 4, 5)

        msk = [ [True, True, True, True, True],
                [True, True, False, False, True],
                [True, False, False, True, True],
                [True, True, False, False, True],
                [True, False, True, False, True],
                [True, True, True, True, True] ]
        #print x
        #print y
        #print msk
        #grid = Grid2D(x, y)
        grid = Grid2D(x, y, mask=msk)
        print grid

        n = grid.get_nodes()
        e = grid.get_edges()
        f = grid.get_faces()
        for i in xrange(7):
            print n[i]
            #print e[i]

        a = np.arange(15).reshape(3,5)
        for v1, v2 in zip(a.flatten(1), a.ravel(order='F')):
            if v1 != v2:
                print v1, v2, 'are not equal'





    test1()

