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
    
    return (nodes, elements,constr_nodes)
    

#n, e, cn = Nodes(2,3)
#
#print n
#print e
#print cn