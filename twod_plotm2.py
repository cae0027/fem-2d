
import matplotlib.pyplot as plt

def twod_plotm2(x,e_conn,symbol):
    '''
    #-----------------------------------------------------------------------#
    #  twod_plotm2.m - Plots linear triagular mesh                          #
    #                  (also works for quadratic straight sided elements)   #
    #                                                                       #
    #  Copyright (c) 2001, Jeff Borggaard, Virginia Tech                    #
    #  Version: 1.0                                                         #
    #                                                                       #
    #  Usage:    [] = twod_plotm2(n_fig,x,e_conn,symbol)                    #
    #                                                                       #
    #  Variables:     n_fig                                                 #
    #                        Makes figure n_fig active                      #
    #                 x                                                     #
    #                        Nodal coordinates                              #
    #                 e_conn                                                #
    #                        Element connectivity                           #
    #                 symbol                                                #
    #                        (optional) symbol for the nodes                #
    #                                                                       #
    #  See twod_plotm1 for linear elements                                  #
    #-----------------------------------------------------------------------#
    '''
    n_elems,ndof = e_conn.shape
  
    if ndof==6:
        for i in range(n_elems):
            plt.plot([x[e_conn[i,0],0], x[e_conn[i,3],0], x[e_conn[i,1],0],
                 x[e_conn[i,4],0], x[e_conn[i,2],0], x[e_conn[i,5],0],x[e_conn[i,0],0]],
                    [x[e_conn[i,0],1], x[e_conn[i,3],1], x[e_conn[i,1],1],
                      x[e_conn[i,4],1], x[e_conn[i,2],1], x[e_conn[i,5],1],
                        x[e_conn[i,0],1]],'b')

    elif ndof==3:
        for i in range(n_elems):
            plt.plot([x[e_conn[i,0],0], x[e_conn[i,1],0],
            x[e_conn[i,2],0], x[e_conn[i,0],0]],
            [x[e_conn[i,0],1], x[e_conn[i,1],1], 
            x[e_conn[i,2],1], x[e_conn[i,0],1]],'b')

    n_nodes = max( x.shape );
    for i in range(n_nodes):
        plt.plot(x[i,0],x[i,1],symbol)
        
    plt.show()


# Test twod_plotm2
import numpy as np

x = np.array([[0.  , 0.  ],
               [0.25, 0.  ],
               [0.5 , 0.  ],
               [0.75, 0.  ],
               [1.  , 0.  ],
               [0.  , 0.5 ],
               [0.25, 0.5 ],
               [0.5 , 0.5 ],
               [0.75, 0.5 ],
               [1.  , 0.5 ],
               [0.  , 1.  ],
               [0.25, 1.  ],
               [0.5 , 1.  ],
               [0.75, 1.  ],
               [1.  , 1.  ]])

e_conn = np.array([[ 0, 12, 10,  6, 11,  5],
                   [ 0,  2, 12,  1,  7,  6],
                   [ 2, 14, 12,  8, 13,  7],
                   [ 2,  4, 14,  3,  9,  8]])

    

if __name__ == '__main__':
    twod_plotm2(x, e_conn, 'go')
    
    
