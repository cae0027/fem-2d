
import numpy as np
from oned_gauss import oned_gauss 

def twod_edge_shape(x_local,r,w):
    '''
    #---------------------------------------------------------------------------------#
    # Description: Computes test functions on an edge at Gauss nodes.                 #
    #                                                                                 #
    # Usage: x_g,w_g,phi = twod_edge_shape(x_local,r,w)                               #
    #                                                                                 #
    #                                                                                 #
    # Inputs:                                                                         #
    #                                                                                 #
    #   x_local: double, (n_dof, 2) coordinates of local edge nodes                   #
    #                                                                                 #
    #   r: double, (n_gauss,1) Gauss nodes on reference interval [-1,1]               #
    #                                                                                 #
    #   w: double, (n_gauss,1) Gauss weights on reference interval [-1,1]             #
    #                                                                                 #
    #                                                                                 #
    # Outputs:                                                                        #
    #                                                                                 #
    #   x_g: double, Gauss nodes on edge                                              #
    #                                                                                 #
    #   w_g: double >0, Gauss weights on edge                                         #
    #                                                                                 #
    #   phi: double, (n_gauss,n_dof) basis functions evaluated at physical Gauss      #
    #       nodes.                                                                    #
    #                                                                                 #
    #                                                                                 #
    # Modified: Hans-Werner van Wyk 6/14/2017                                         #
    # --------------------------------------------------------------------------------#
    '''

    x = x_local
    n, _ = x.shape

    # Edge coordinates 
    x1 = x[0,0];  x2 = x[1,0]
    y1 = x[0,1];  y2 = x[1,1]

    # Gauss nodes
    rule = len(r)
    x_g = np.zeros((rule, 2))
    x_g[:,0] = 0.5*(x1+x2) + r*0.5*(x2-x1)
    x_g[:,1] = 0.5*(y1+y2) + r*0.5*(y2-y1)

    # Gauss weights
    w_g = 0.5*w*np.sqrt((x2-x1)**2 + (y2-y1)**2)

    if n == 2:
        # Linear element
        phi = np.zeros((rule, n))
        phi[:,0] = (1 - r)/2
        phi[:,1] = (1 + r)/2
       
    elif n == 3:
        # Quadratic element    
        phi = np.zeros((rule, n)) 
        phi[:,0] = 0.5*r*(r-1)
        phi[:,1] = 0.5*r*(r+1)
        phi[:,2] = -(r+1)*(r-1)
    
    elif n == 4:
        # Cubic element
        phi = np.zeros((rule,n))
        phi[:,0] = -9/2*(r-1/3)*(r-2/3)*(r-1)
        phi[:,1] =   9/2*r*(r-1/3)*(r-2/3)
        phi[:,2] =  27/2*r*(r-2/3)*(r-1)
        phi[:,3] = -27/2*r*(r-1/3)*(r-1)

    else:
        raise Exception('Only linear, quadratic, and cubic elements are supported')
    return x_g, w_g, phi

#
# Test twod_edge_shape
#
n_dof = 2                                        # n_dof for quadratic element is 3
x_local = x[index_b[:n_dof],:]                   # The first boundary element nodes
r, w = oned_gauss(7)

if __name__ == '__main__':
    x_g, w_g, phi = twod_edge_shape(x_local,r,w)
