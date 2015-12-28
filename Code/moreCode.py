
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import math





        

def dL(C,L,dl,bl):
    global coral_harv
    """
    return the derivative of r evaulated as a
    function of *r* and *f*.  The function should work whether *r* and *f*
    are scalars, 1D arrays or 2D arrays.  The return value should have
    the same dimensionality (shape) as the inputs *r* and *f*.
    """

    return (bl*lion_attack*C*L) - dl*L -  grouper_harv*L  #- hLH(L)  - hLG(L)
        

    
    

def dC(C,L,RC,K):
    """
    return the derivative of C evaulated as a
    function of *r* and *f*.  The function should work whether *r* and *f*
    are scalars, 1D arrays or 2D arrays.  The return value should have
    the same dimensionality (shape) as the inputs *r* and *f*. Calls function to 
    predation rates at each time step. 
    """    
    
    return  C*(RC) - C*RC*C/K - ((C)*L*lion_attack)/(1+C*lion_attack*.4)

def derivs(state,t,dl,RC, K,bl):
    global coral_harv, lion_attack, grouper_harv
    """
    Return the derivatives of R and F, stored in the *state* vector::

       state = [R, F]

    The return data should be [dR, dF] which are the derivatives of R
    and F at position state and time *t*
    """
    
    C, L ,G = state   
    print "C =", C, "L=:", L, "G= ", G, "t= ", t
#########################################################
    """
    return the derivative of r evaulated as a
    function of *r* and *f*.  The function should work whether *r* and *f*
    are scalars, 1D arrays or 2D arrays.  The return value should have
    the same dimensionality (shape) as the inputs *r* and *f*.
    """
    ########################################################### 
    g= C*(RC) - C*RC*C/K - ((C)*L*lion_attack)/(1+C*lion_attack*.4)
    #####################################
    f= (bl*lion_attack*C*L) - dl*L -  (L*(2*G)/(1+.91*G))
    #########################################################################
    h = (G*(.8) - G*.8*G/3000)
    """
    return the derivative of C evaulated as a
    function of *r* and *f*.  The function should work whether *r* and *f*
    are scalars, 1D arrays or 2D arrays.  The return value should have
    the same dimensionality (shape) as the inputs *r* and *f*. Calls function to 
    predation rates at each time step. 
    """    
    
    return g, f,h

def main():
    global grouper_harv
    global  human_attack, human_time
    global prevC
    global lion_attack
    # the initial conditions
    C0 = 212
    L0 = 45
    G0 = 178
    
    K = 800
    RC = .447
    dl = .3 #roughly the lifespan of a fish that live 15-20 years
    bl = .05
    grouper_harv = .349
    lion_attack = 8.6/240

    
    t = np.linspace(0, 5, num=300)
    y0 = [C0, L0, G0]  # the initial 
    
    # integrate your ODE using scipy.integrate.  
    # Lionfish should consume  /year
    y = integrate.odeint(derivs, y0, t, args=(dl,RC,K,bl))
    
      
     
    plt.figure()
    p1, = plt.plot(t, y[:,0], 'k')#coralfish
    p2, = plt.plot(t, y[:,1], 'k' ,lw=2)#Lionfish
    p3, = plt.plot(t, y[:,2], 'k+' ,lw=1)#grouper
    #p3, = plt.plot( grouper_harv, 'k+')
    plt.legend([p1,p2,p3], ["Coral-Fish","Lionfish", "Grouper"], loc="best")
    plt.grid()
#     
#    
#      
#     #the return value from the integration is a Nx2 array.  Extract it
#     # into two 1D arrays caled r and f using numpy slice indexing
#        
#       
#     # time series plot:
#     # funciton of time
#       
#       
#     # phase-plane plot:
#     # make sure you include and xlabel, ylabel and title
#     plt.figure()
#     plt.plot(t, y, 'k--')
#     plt.title('phase plane')
#       
#       
#     #gradient nonsense
#     rmax = 1.01 * t.max()
#     Cmax = 1.01 * t.max()
#     R, F = np.meshgrid(np.arange(-1, rmax), np.arange(-1, Cmax))
#     dR = dL(R, F,dl,bl)
#     dF = dC(R, F, RC,K)
#     plt.quiver(R, F, dR, dF)
#       
#     #Nullcline non-sense
#     R, F = np.meshgrid(np.arange(-1, rmax, 0.1), np.arange(-1, Cmax, 0.1))
#     dR = dL(R, F,dl,bl)
#      
#     dF = dC(R, F,RC,K)
#     plt.contour(R, F, dR, levels=[0], linewidths=.5, colors='k')
#     plt.contour(R, F, dF, levels=[0], linewidths=.5, colors='k')
#     #plt.hlines(0, 0, 30)
#     #plt.vlines(Cbar, -20, 20, linewidths=2, color='k')
#     #plt.vlines(((C0-Cbar)/B), -20, 20, linewidths=2, color='k')
#     plt.legend()
    plt.show()
main()