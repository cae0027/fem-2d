
import numpy as np

def twod_mesh(x_l,x_r,y_l,y_r,etype,n_nodesx,n_nodesy, *args):
    '''
    --------------------------------------------------------------------------------
    #  twod_mesh.m - Generate a rectangular mesh with a prescribed density.        #
    #                This routine returns nodal coordinates,                       #
    #                connectivity, and the nodal indices of of an element          #
    #                                                                              #
    #  Usage:    [x,e_conn,index_b] = twod_mesh(x_l,x_r,y_l,y_r,etype,...          #
    #                                           n_nodex,n_nodesy,fig_num)          #
    #                                                                              #
    #  Variables:     (x_l,y_l)                                                    #
    #                        coordinates of lower left corner                      #
    #                 (x_r,y_r)                                                    #
    #                                                                              #
    #                 etype                                                        #
    #                        element type                                          #
    #                          'linear', 'quadratic', 'cubic'                      #
    #                                                                              #
    #                 n_nodesx, n_nodesy                                           #
    #                        number of nodes in x and y directions                 #
    #                        (must be compatible with element type)                #
    #                                                                              #
    #                 args                                                         #
    #                     add optional argument 'plot' to display plot             #
    #                                                                              #
    #  Outputs:                                                                    #
    #                 x                                                            #
    #                        node coordinates of mesh                              #
    #                 e_conn                                                       #
    #                        element connectivity of mesh                          #
    #                 index_b                                                      #
    #                        node numbers of boundary nodes                        #
    --------------------------------------------------------------------------------
    '''
    if etype == 'linear':
        n_nodes = n_nodesx*n_nodesy
        n_elements = 2*(n_nodesx-1)*(n_nodesy-1)
        
        # Generate node coordinates
        dx = (x_r-x_l)/(n_nodesx-1)
        dy = (y_r-y_l)/(n_nodesy-1)
        
        x = np.zeros((n_nodes,2))
        for i in range(n_nodesx):
            for j in range(n_nodesy):
                k = i+j*n_nodesx;
                x[k,0] = x_l + dx*i
                x[k,1] = y_l + dy*j


    # Generate element connectivity
        e_conn = np.zeros((n_elements, 3)).astype(int);
        ie     = -1
        for j in range(n_nodesy-1):
            for i in range(n_nodesx-1):
                ie = ie + 1 
                k = i + j*n_nodesx
                e_conn[ie,:] = [k, k+1+n_nodesx, k+n_nodesx]
                ie = ie + 1
                e_conn[ie,:] = [k, k+1, k+1+n_nodesx]

    elif etype == 'quadratic':
        #  Generate a mesh with quadratic elements
        n_nodes  = n_nodesx*n_nodesy;
        n_elements = int(2*(n_nodesx-1)*(n_nodesy-1)/4)

        dx = (x_r-x_l)/(n_nodesx-1)
        dy = (y_r-y_l)/(n_nodesy-1)
        
        x = np.zeros((n_nodes,2))
        for i in range(n_nodesx):
            for j in range(n_nodesy):
                k = i+j*n_nodesx;
                x[k,0] = x_l + dx*i
                x[k,1] = y_l + dy*j

        #  Generate element connectivity
        e_conn = np.zeros((n_elements, 6)).astype(int)
        ie     = -1
        for j in range(0,n_nodesy-1,2):
            for i in range(0,n_nodesx-1,2):
                ie = ie + 1
                k = i + j*n_nodesx
                e_conn[ie,:] = [k, k+2+2*n_nodesx, k+2*n_nodesx, 
                                k + 1 + n_nodesx, k+1+2*n_nodesx, k+n_nodesx ]
                ie = ie + 1
                e_conn[ie,:] = [k, k+2, k+2+2*n_nodesx, 
                                k+1, k+2+n_nodesx, k+1+n_nodesx ]

    elif etype == 'cubic':
        raise Exception('twod_mesh: cubic elements are not currently implemented')
    else:
        raise Exception('twod_mesh: {} is not a valid element type'.format(etype))
    
    # Get boundary indices
    row_1 = np.arange(1,n_nodesx-1)
    row_2 = np.arange(0,(n_nodesy-1)*n_nodesx, n_nodesx)
    row_3 = np.arange(n_nodesx-1, n_nodesy*n_nodesx, n_nodesx)
    row_4 = np.arange((n_nodesy-1)*n_nodesx, n_nodes-1)

    index_b = np.concatenate((row_1, row_2, row_3, row_4), axis=0)
    
    if len(args) >= 1:
        twod_plotm2(x,e_conn,'ro')
        
    return x, e_conn, index_b


# Test twod_mesh
x_min = 0; x_max = 1; nx = 5 
y_min = 0; y_max = 1; ny = 3
etype = 'quadratic'

if __name__ == '__main__':
    x,e_conn,index_b = twod_mesh(x_min,x_max,y_min,y_max,etype,
                                           nx,ny)
  
