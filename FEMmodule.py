# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 09:47:52 2015

@author: claudio
"""

from __future__ import division
import numpy as np

def Nodes(nx, ny, l=10, w=2):
    """
    function to compute nodes and elements coordinates
    input 
    nx : number of elements in x direction
    ny : number of elements in y direction
    l :  length of plate (default = 10)
    w :  width of plate (default = 2)
    
    output
    - matrix of nodes coordinates
    - matrix of elements (which nodes belongs to each element)
    - matrix of constrained elements
    """
    nodes = np.zeros(( (nx+1)*(ny+1),2 ))
    
    nodes[:,0] = np.tile(np.linspace(0,l,nx+1),ny+1)
    nodes[:,1] = np.tile(np.linspace(0,w,ny+1),(nx+1,1)).reshape((nx+1)*(ny+1), order='F') 
    
    elements = np.zeros((nx*ny,4))
    
    elements[:,0] =  np.reshape(np.tile(np.arange(1,nx+1),(ny,1)).T + \
                                np.tile(np.arange(0,nx*ny+(ny-nx)+1,nx+1),(nx,1)), \
                                        (nx*ny),order='F')
    elements[:,1] = elements[:,0] + 1
    elements[:,2] = elements[:,0] + 1 + nx + 1
    elements[:,3] = elements[:,0] + nx + 1
    
    constr_nodes = np.arange(1,(nx+1)*(ny+1)+1, nx+1 )
    
    return (nodes, elements, constr_nodes)
    

#n, e, cn = Nodes(2,3)
#
#print n
#print e
#print cn


class GaussIntegr2D2x2:
    npoints = 4
    points = 1./np.sqrt(3.)*np.array([[-1, -1],[1,-1],[-1,1],[1,1]])
    weights = np.ones((4,1))

    def Integrate(self, f):
        return np.sum([f(self.points[i,:])*self.weights[i] for i in range(self.npoints)])



class PlaneStressElasticTensor:
    def __init__(self, E, nu):
        self.ElTens = np.array([[1., nu, 0.],[nu, 1., 0.],[0., 0., (1. - nu) / 2.]]) * (E / (1 - nu**2))



class BilinearWeights:
    """
    Bilinear interpolation functions
    """
    ndims = 2
    nnodes = 4
    
    def Weights(self, x):
        """
        bilinear weights
        x must be a 1x2 vector
        weights is a 4x1 vector [N1 N2 N3 N4]
        """
        weights = np.zeros(self.nnodes)
        weights[0] = (x[0]-1)*(x[1]-1)*0.25
        weights[1] = (x[0]+1)*(x[1]-1)*0.25
        weights[2] = (x[0]+1)*(x[1]+1)*0.25
        weights[3] = (x[0]-1)*(x[1]+1)*0.25
        
        return weights
    
    def Ders(self, x):
        """
        Bilinear weight derivatives
        - x must be a (1x2) vector
        - ders is a (2x4) vector: ders(i,j) is the partial derivative of weight j with respect to the coordinate i
        """
        ders = np.zeros((self.ndims, self.nnodes))
        ders[0, 0] = (x[1] - 1.) *  0.25	#-1 -1
        ders[0, 1] = (x[1] - 1.) * -0.25	# 1 -1
        ders[0, 2] = (x[1] + 1.) *  0.25	# 1  1
        ders[0, 3] = (x[1] + 1.) * -0.25    #-1  1
        ders[1, 0] = (x[0] - 1.) *  0.25    #-1 -1
        ders[1, 1] = (x[0] + 1.) * -0.25    # 1 -1
        ders[1, 2] = (x[0] + 1.) *  0.25	# 1  1
        ders[1, 3] = (x[0] - 1.) * -0.25	#-1  1
        return ders
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        