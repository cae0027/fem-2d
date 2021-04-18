
import numpy as np

def twod_bilinear( kernel, phi, test, w_g ):
    
    '''
    #--------------------------------------------------------------------------#
    #  twod_bilinear.py - routine to compute \int{ kernel*phi*test }           #
    #                                                                          #
    #  Usage:    M = twod_bilinear(kernel, phi, test, w_g)                     #
    #                                                                          #
    #  Variables:     kernel                                                   #
    #                        Kernel function in the integral evaluated         #
    #                        at the Gauss points                               #
    #                                                                          #
    #                 phi                                                      #
    #                        matrix of element test functions evaluated        #
    #                        at the Gauss points (dim: n_gauss, n_dof)         #
    #                                                                          #
    #                 test                                                     #
    #                        matrix of test functions evaluated at the         #
    #                        Gauss points (dim: n_gauss, n_dof)                #
    #                                                                          #
    #                 w_g                                                      #
    #                        Row vector of Gauss weights                       #
    #--------------------------------------------------------------------------#
    #   [n_gauss,n_row] = size(test)                                           #
    #   [n_g1   ,n_col] = size(phi )                                           #
    #                                                                          #
    #   M = np.zeros((n_row,n_col))                                            #
    #   for i in range(n_row):                                                 #
    #     for j=1:n_col                                                        #
    #       M[i,j] = np.dot(( w_g.T    * test[:,i]' ),( kernel * phi[:,j] )) #
    #--------------------------------------------------------------------------#
    '''                                                        
    M = np.dot( np.dot(test.T, np.diag( w_g*kernel )), phi)
    return M



# Test twod_bilinear 
x_local = x[e_conn[0,:],:]                   # The first element nodes
#
# x and e_conn are from 
#
# x, e_conn, _ = twod_mesh(x_l,x_r,y_l,y_r,etype,n_nodesx,n_nodesy, *args)

q = lambda x,y:x+y
x_g,w_g,phi,p_x,p_y = twod_shape(x_local, r, s, w)
kernel = q(r,s)
test = phi                                   # For Galerkin method

if __name__ == '__main__':
    M = twod_bilinear( kernel, phi, test, w_g )
