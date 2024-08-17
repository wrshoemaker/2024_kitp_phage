import sys
import numpy
from scipy import integrate
from scipy import special
import config
import matplotlib.pyplot as plt




phi_r_max = 0.55


v = 3 # nutritional efficiency
# convert OD600 to cells 
od600_to_cell_density = 609375000
#n_cell_0 = od600_to_cell_density*0.01
n_cell_0 = od600_to_cell_density*0.01

#n_phage_0 = 1e3
I = 1e-6
tau_lysis = 0.88 #(hours)

delta_phage = 1e-8



def lotka_volterra_V(lambda_cell, lambda_phage, n_cell, n_phage):

    V = I*tau_lysis*lambda_phage*n_cell - delta_phage*numpy.log(n_cell) + I*n_phage - lambda_cell*numpy.log(n_phage)

    return V





def solve_n_cell_max(lambda_cell, lambda_phage, n_cell_0, n_phage_0):

    V = lotka_volterra_V(lambda_cell, lambda_phage, n_cell_0, n_phage_0)

    #lambert_w_term = delta_phage*(lambda_phage**-1) * (I**-1) * (tau_lysis**-1) * numpy.exp( (V - lambda_cell*(numpy.log(lambda_cell/I) -1)) / (I*tau_lysis*lambda_phage) )
    lambert_w_term =  (delta_phage/(lambda_phage*I*tau_lysis)) *  numpy.exp( (V - lambda_cell*(numpy.log(lambda_cell/I) -1)) / (I*tau_lysis*lambda_phage) )

    #print((V - lambda_cell*(numpy.log(lambda_cell/I) -1)) / (I*tau_lysis*lambda_phage) )
    #print((delta_phage/(lambda_phage*I*tau_lysis)) )
    #print(lambert_w_term)   
    W = special.lambertw(lambert_w_term, k=0)

    W_real = numpy.real(W)

    n_cell_max = W_real*lambda_phage * I * tau_lysis * (delta_phage**-1)

    return n_cell_max


def lambert_approximation(log_x):
    
    # assume you already have the logarithm
    l_1 = log_x
    l_2 = numpy.log(l_1)

    term_5 = l_2 * (2*(l_2**2) - 9*l_2 + 6)/ (6*(l_1**3))

    return l_1 - l_2 + (l_2/l_1) + (l_2*(l_2 - 2))/(2*(l_1**2)) + term_5



n_cell_0 = 1e7
#n_phage_0 = 1e5

lambda_phage = 10
lambda_cell = 2


def make_figure():

    #n_cell_max = solve_n_cell_max(lambda_cell, lambda_phage, n_cell_0, n_phage_0)

    #x_logrange = numpy.logspace(-3,3, base=10, num=50)

    #for x in x_logrange:

    #    print(special.lambertw(x, k=0))


    n_phage_0_all = numpy.logspace(-1,10, base=10, num=50)
    od_max_all = []

    for n_phage_0_i in n_phage_0_all:


        V = lotka_volterra_V(lambda_cell, lambda_phage, n_cell_0, n_phage_0_i)

        log_term_in_lambert_w = numpy.log(delta_phage) - numpy.log(lambda_phage) - numpy.log(I) - numpy.log(tau_lysis) + ((V - lambda_cell*(numpy.log(lambda_cell/I) -1)) / (I*tau_lysis*lambda_phage))
        lambert_w_approx = lambert_approximation(log_term_in_lambert_w)

        n_cell_max = ((I*tau_lysis*lambda_phage)/delta_phage) * lambert_w_approx

        od_max_all.append(n_cell_max/od600_to_cell_density)
        


    fig = plt.figure(figsize = (4, 4))
    fig.subplots_adjust(bottom= 0.15)

    ax = plt.subplot2grid((1, 1), (0, 0), colspan=1)

    ax.plot(n_phage_0_all, od_max_all)
    ax.set_xscale('log', basex=10)
    ax.set_xlabel('Log of initial phage concentration')
    ax.set_ylabel('Lysis biomass')


    fig.subplots_adjust(hspace=0.3, wspace=0.25)
    fig_name = "%stheory_w_death.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
    plt.close()





