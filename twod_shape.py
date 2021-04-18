
import numpy as np
import warnings

def twod_shape(x_local,r,s,w):
    '''
    ----------------------------------------------------------------------------------
    #  twod_shape.m - computes test functions and derivatives on an                  #
    #                 element given element coordinates and Gauss points.            #
    #                                                                                #
    #                 ! Note: optimized for straight-sided elements.  Use            #
    #                 ! `twod_shapeiso' for isoparametric elements.                  #
    #                                                                                #
    #  Usage:    [x_g,w_g,phi,p_x,p_y] = twod_shape(x_local,r,s,w)                   #
    #                                                                                #
    #  Variables:     x_local                                                        #
    #                        Coordinates of the element nodes                        #
    #                 (r,s)                                                          #
    #                        Coordinates of Gauss points in unit triangle            #
    #                 w                                                              #
    #                        Gauss weights associated with (r,s)                     #
    #                                                                                #
    #                 x_g                                                            #
    #                        Coordinates of Gauss points in the element              #
    #                 w_g                                                            #
    #                        Gauss weights scaled by the element Jacobian            #
    #                 phi                                                            #
    #                        Value of element shape functions at x_g                 #
    #                 p_x                                                            #
    #                 p_y                                                            #
    #                        First spatial derivatives of phi                        #
    #                                                                                #
    ----------------------------------------------------------------------------------
    '''
    x        = x_local
    n,ncoord = x.shape

    # fprintf('Space Dimension           : %i\n',ncoord)
    # fprintf('Number of shape functions : %i\n',n)
    
    if ncoord != 2:
        warnings.warn("You are messing with a wrong dimension, please stick to 2D!")

    rule = len(r)
   
    # Compute (r,s) -> (x,y) transformation for straight-sided elements
    c0 =  x[0,:]
    c1 = -x[0,:] + x[1,:]
    c2 = -x[0,:]          + x[2,:]

    x_g = np.zeros(( rule, n ))
    x_g[:,0] = c0[0] + c1[0]*r + c2[0]*s
    xr  = c1[0]
    xs  = c2[0]

    x_g[:,1] = c0[1] + c1[1]*r + c2[1]*s
    yr  = c1[1]
    ys  = c2[1]

    # Compute the Jacobian of the (r,s) -> (x,y) transformation
    jac = xr*ys - yr*xs
    w_g = jac*w

    rx  = ys/jac
    sx  =-yr/jac
    ry  =-xs/jac
    sy  = xr/jac

    if n == 3:
        # Compute shape function and derivatives at Gauss points
        phi = np.zeros((rule,n))
        phi[:,0] = 1. - r  - s
        phi[:,1] =      r     
        phi[:,2] =           s

        p_x = np.zeros((rule,n))
        p_x[:,0] =    - rx - sx
        p_x[:,1] =      rx     
        p_x[:,2] =           sx
        
        p_y = np.zeros((rule,n))
        p_y[:,0] =    - ry - sy
        p_y[:,1] =      ry     
        p_y[:,2] =           sy

    elif n == 6:
        # Compute shape function and derivatives at Gauss points
        phi = np.zeros((rule,n))
        phi[:,0] = 1. - 3.*r - 3.*s + 2.*r*r + 4.*r*s + 2.*s*s
        phi[:,1] =    - 1.*r        + 2.*r*r                    
        phi[:,2] =           - 1.*s                     + 2.*s*s
        phi[:,3] =      4.*r        - 4.*r*r - 4.*r*s          
        phi[:,4] =                              4.*r*s          
        phi[:,5] =             4.*s           - 4.*r*s - 4.*s*s

        p_x = np.zeros(( rule,n ));
        p_x[:,0] = ( -3. + 4.*r + 4.*s )*rx + ( -3. + 4.*r + 4.*s )*sx
        p_x[:,1] = ( -1. + 4.*r        )*rx                            
        p_x[:,2] =                             ( -1.        + 4.*s )*sx
        p_x[:,3] = (  4. - 8.*r - 4.*s )*rx + (     - 4.*r        )*sx
        p_x[:,4] = (              4.*s )*rx + (       4.*r        )*sx
        p_x[:,5] = (            - 4.*s )*rx + (  4. - 4.*r - 8.*s )*sx

        p_y = np.zeros(( rule,n ));
        p_y[:,0] = ( -3. + 4.*r + 4.*s )*ry + ( -3. + 4.*r + 4.*s )*sy
        p_y[:,1] = ( -1. + 4.*r        )*ry                            
        p_y[:,2] =                             ( -1.        + 4.*s)*sy
        p_y[:,3] = (  4. - 8.*r - 4.*s )*ry + (     - 4.*r        )*sy
        p_y[:,4] = (              4.*s )*ry + (       4.*r        )*sy
        p_y[:,5] = (            - 4.*s )*ry + (  4. - 4.*r - 8.*s )*sy
    
    elif n == 7:
        # Compute shape function and derivatives at Gauss points
        phi = np.zeros((rule,n))
        phi[:,0] = (1-r-s)*(2.*(1-r-s)-1) +  3.*(1.-r-s)*r*s
        phi[:,1] = r*(2.*r-1)             +  3.*(1.-r-s)*r*s
        phi[:,2] = s*(2.*s-1)             +  3.*(1.-r-s)*r*s
        phi[:,3] = 4.*(1-r-s)*r           - 12.*(1.-r-s)*r*s
        phi[:,4] = 4.*r*s                 - 12.*(1.-r-s)*r*s
        phi[:,5] = 4.*s*(1-r-s)           - 12.*(1.-r-s)*r*s
        phi[:,6] = 27.*(1-r-s)*r*s                   

        p_r = np.zeros((rule,n))
        p_r[:,0] = -3 + 4.*r + 7.*s - 6.*r*s - 3.*(s**2)
        p_r[:,1] = -1 + 4.*r + 3.*s - 6.*r*s - 3.*(s**2)
        p_r[:,2] =             3.*s - 6.*r*s - 3.*(s**2)
        p_r[:,3] =  4 - 8.*r -16.*s +24.*r*s +12.*(s**2)
        p_r[:,4] =           - 8.*s +24.*r*s +12.*(s**2)
        p_r[:,5] =           -16.*s +24.*r*s +12.*(s**2)
        p_r[:,6] =            27.*s -54.*r*s -27.*(s**2)

        p_s = np.zeros((rule,n))
        p_s[:,0] = -3 + 7.*r + 4.*s - 6.*r*s - 3.*(r**2)
        p_s[:,1] =      3.*r        - 6.*r*s - 3.*(r**2)
        p_s[:,2] = -1 + 3.*r + 4.*s - 6.*r*s - 3.*(r**2)
        p_s[:,3] =    -16.*r        +24.*r*s +12.*(r**2)
        p_s[:,4] =    - 8.*r        +24.*r*s +12.*(r**2)
        p_s[:,5] =  4 -16.*r - 8.*s +24.*r*s +12.*(r**2)
        p_s[:,6] =     27.*r        -54.*r*s -27.*(r**2)

        p_x = p_r*rx + p_s*sx
        p_y = p_r*ry + p_s*sy

    else:
        raise Exception('Elements with {} interior nodes are not currently supported'
                       .format(x.shape[0]))
    return x_g,w_g,phi,p_x,p_y



# Test twod_shape
x_local = x[e_conn[0,:],:]  # This is the first element 
#
# x and e_conn are from 
#
# x, e_conn, _ = twod_mesh(x_l,x_r,y_l,y_r,etype,n_nodesx,n_nodesy, *args)

if __name__ == '__main__':
    x_g,w_g,phi,p_x,p_y = twod_shape(x_local, r, s, w)
    
    
