
def twod_f_int(Ff,test,w_g):
    '''
    #----------------------------------------------------------------------#
    #  twod_f_int.py - routine to compute \int{ f*test }                   #
    #                                                                      #
    #  Usage:    F = twod_f_int( Ff, test, w_g )                           #
    #                                                                      #
    #  Variables:     Ff                                                   #
    #                        Function values at the Gauss points           #
    #                                                                      #
    #                 test                                                 #
    #                        matrix of test functions evaluated at the     #
    #                        Gauss points (dim: n_gauss, n_dof)            #
    #                                                                      #
    #                 w_g                                                  #
    #                        Row vector of Gauss weights                   #
    #----------------------------------------------------------------------#
    # n_gauss,n_dof = test.shape                                           #
    #                                                                      #
    # F = np.zeros(n_dof)                                                  #
    # for j in range(n_dof):                                               #
    #       F(j) = test(:,j)' * ( w_g .* Ff )                              #
    #                                                                      #
    #----------------------------------------------------------------------#
    '''
    F = sum(np.dot(test.T, w_g*Ff))
    return F



# Test twod_f_int 
f = lambda x,y:x+y
x_g,w_g,phi,p_x,p_y = twod_shape(x_local, r, s, w)
Ff = f(r,s)

if __name__ == '__main__':
    F = twod_f_int(Ff, phi, w_g)

