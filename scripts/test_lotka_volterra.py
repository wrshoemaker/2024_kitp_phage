
import config
import numpy
import matplotlib.pyplot as plt
import matplotlib.colors as colors


from scipy import integrate
# Definition of parameters
a = 2
b = 1e-6
c = 0 #1e-10
d = 50
def dX_dt(X, t=0):
    """ Return the growth rate of fox and rabbit populations. """
    return numpy.array([ a*X[0] -   b*X[0]*X[1] ,
                  -c*X[1] + d*b*X[0]*X[1] ])



#!python
X_f0 = numpy.array([     0. ,  0.])
X_f1 = numpy.array([ c/(d*b), a/b])
all(dX_dt(X_f0) == numpy.zeros(2) ) and all(dX_dt(X_f1) == numpy.zeros(2)) # => True


#!python
def d2X_dt2(X, t=0):
    """ Return the Jacobian matrix evaluated in X. """
    return numpy.array([[a -b*X[1],   -b*X[0]     ],
                  [b*d*X[1] ,   -c +b*d*X[0]] ])


#A_f0 = d2X_dt2(X_f0)
#A_f1 = d2X_dt2(X_f1)                    # >>> array([[ 0.  , -2.  ],
                                        #            [ 0.75,  0.  ]])
# whose eigenvalues are +/- sqrt(c*a).j:
#lambda1, lambda2 = numpy.linalg.eigvals(A_f1) # >>> (1.22474j, -1.22474j)
# They are imaginary numbers. The fox and rabbit populations are periodic as follows from further
# analysis. Their period is given by:
#T_f1 = 2*numpy.pi/abs(lambda1)  


def make_plot():

    t = numpy.linspace(0, 15,  1000)   
    phage_0_array = numpy.logspace(0, 8, base=10, num=20)

    max_cell_t_all = []
    for phage_0 in phage_0_array:

        X0 = numpy.array([609375000*0.01, phage_0])         

        X, infodict = integrate.odeint(dX_dt, X0, t, full_output=True)
        #infodict['message']   



        cell_t, phaeg_t = X.T
        max_cell_t_all.append(max(cell_t[1:]))

        #print(numpy.argmax(cell_t))

        print(min(cell_t[1:]))


    fig = plt.figure(figsize = (4, 4))
    fig.subplots_adjust(bottom= 0.15)

    ax = plt.subplot2grid((1, 1), (0, 0), colspan=1)

    ax.plot(phage_0_array, max_cell_t_all, lw=2)

    ax.set_xscale('log', basex=10)
    #ax.set_yscale('log', basey=10)

    ax.set_xlabel('Initial phage concentration (PFU/mL)', fontsize=10)
    ax.set_ylabel('Biomass at lysis (OD600)', fontsize=12)


    fig.subplots_adjust(hspace=0.35, wspace=0.40)
    fig_name = "%sphage_0_vs_max_od.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



#f1 = p.figure()
#p.plot(t, rabbits, 'r-', label='Rabbits')
#p.plot(t, foxes  , 'b-', label='Foxes')
#p.grid()
#p.legend(loc='best')
##p.xlabel('time')
#p.ylabel('population')
#p.title('Evolution of fox and rabbit populations')
#f1.savefig('%srabbits_and_foxes_1.png' % config.analysis_directory)