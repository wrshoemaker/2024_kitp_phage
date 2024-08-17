import config
import numpy
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import utils

# density. phage/mL
#n_phage_0 = 1e5
#n_cell_0 = 1e7


lambda_phage = 2
#lambda_phage = 0.01

#lysis_rate_per_phage_per_time = 1e-7

phi_r_max = 0.55



v = 3 # nutritional efficiency
# convert OD600 to cells 
od600_to_cell_density = 609375000
#n_cell_0 = od600_to_cell_density*0.01
n_cell_0 = od600_to_cell_density*0.01

#n_phage_0 = 1e3
I = 1e-8
#I = 0.00000001
tau_lysis = 0.88 #(hours)



def calculate_lambda_cell(phi_r_0):

    lambda_cell = (phi_r_max  - phi_r_0)*v

    return lambda_cell


#def calculate_lambda_phage(phi_r_0):

#    lambda_phage = (phi_r_max  - phi_r_0)*5

#    return lambda_phage



#lambda_cell = 2


def calculate_log_od600_max(lambda_cell, lambda_phage, n_phage_0):


    t_peak = (1/lambda_cell) * numpy.log(1 + (lambda_cell/(I*tau_lysis*lambda_phage*n_cell_0) * numpy.log(lambda_cell/(n_phage_0*I)) ))
    log_n_max = numpy.log(n_cell_0) + lambda_cell * t_peak

    # convert to OD600
    log_od600_max = log_n_max - numpy.log(od600_to_cell_density)

    return log_od600_max



def calculate_log_od600_max_burst_limited(lambda_cell, lambda_phage, n_phage_0):


    t_peak = numpy.log((lambda_cell/I) * ( (lambda_phage**-1) + (n_cell_0/n_phage_0) - (lambda_cell/I) )) / (I*lambda_phage - lambda_cell)
    n_max = n_cell_0*numpy.exp(lambda_cell*t_peak) - (I*n_phage_0*(numpy.exp(I*lambda_phage*t_peak) - numpy.exp(lambda_cell*t_peak) ) / ((I*lambda_phage - lambda_cell))  )

    #print( ( (lambda_phage**-1) + (n_cell_0/n_phage_0) - (lambda_cell/I) ))
    # convert to OD600
    log_od600_max = numpy.log(n_max) - numpy.log(od600_to_cell_density)

    return log_od600_max




def calculate_od600_max_2nd_case(lambda_cell, lambda_phage, n_phage_0, delta=1e-7):

    t_peak = (delta**-1)*numpy.log(I*n_phage_0*(lambda_cell**-1))

    n_max = n_cell_0 * numpy.exp(lambda_cell*t_peak + I*n_phage_0*(delta**-1)*(numpy.exp(-1*delta*t_peak) - 1))

    od600_max = n_max/od600_to_cell_density

    return od600_max


def lotka_volterra_V(lambda_cell, lambda_phage, n_cell, n_phage, delta_phage=1e-7):

    V = I*tau_lysis*lambda_phage*n_cell - delta_phage*numpy.log(n_cell) + I*n_phage - lambda_cell*numpy.log(n_phage)

    return V


def calculate_od600_max_lv_no_death(lambda_cell, lambda_phage, n_phage_0):

    V_0 = lotka_volterra_V(lambda_cell, lambda_phage, n_cell_0, n_phage_0, delta_phage=0)

    n_max = (V_0 - lambda_cell*(1 - numpy.log(lambda_cell/I))) * ((I*tau_lysis*lambda_phage)**-1)

    od600_max = n_max/od600_to_cell_density

    return od600_max









def plot_change_in_cell_growth_only():

    # convert to natural log with numpy.log(10)
    n_phage_0_all =  numpy.logspace(0, 8, base=10, num=100, endpoint=True)
    #n_phage_0_all =  numpy.logspace(0, 6, base=10, num=100, endpoint=True)

    fig = plt.figure(figsize = (4, 4))
    fig.subplots_adjust(bottom= 0.15)

    ax = plt.subplot2grid((1, 1), (0, 0), colspan=1)


    for phi_r_0_i_idx, phi_r_0_i in enumerate([0, 0.1, 0.2]):

        lambda_cell_i = calculate_lambda_cell(phi_r_0_i)
        #lambda_phage_i = calculate_lambda_phage(phi_r_0_i)
        log_od600_max_all = [calculate_log_od600_max(lambda_cell_i, lambda_phage, n) for n in n_phage_0_all]
        #od600_max_2nd_case_all = [calculate_od600_max_2nd_case(lambda_cell_i, lambda_phage, n) for n in n_phage_0_all]
        #log_od600_max_burst_limited_all = [calculate_log_od600_max_burst_limited(lambda_cell_i, lambda_phage, n) for n in n_phage_0_all]
        #od600_max_lv_no_death = [calculate_od600_max_lv_no_death(lambda_cell_i, lambda_phage, n) for n in n_phage_0_all]
        #print(od600_max_lv_no_death)

        print(lambda_cell_i/I)

        fract_nonfunctional = ( 1 - ((phi_r_max - phi_r_0_i)/phi_r_max)) * 100
        #print(fract_nonfunctional)
        #ax.plot(n_phage_0_all, numpy.exp(log_od600_max_all), label=int(fract_nonfunctional), c = utils.rgb_red_antibiotic(phi_r_0_i_idx))
        #ax.plot(n_phage_0_all, od600_max_2nd_case_all, label=int(fract_nonfunctional), c = utils.rgb_red_antibiotic(phi_r_0_i_idx))
        ax.plot(n_phage_0_all, numpy.exp(log_od600_max_all), label=int(fract_nonfunctional), c = utils.rgb_red_antibiotic(phi_r_0_i_idx))
        #ax.plot(n_phage_0_all, numpy.exp(log_od600_max_burst_limited_all), label=int(fract_nonfunctional), c = utils.rgb_red_antibiotic(phi_r_0_i_idx))

        
        #ax.plot(n_phage_0_all, log_od600_max_all, label=int(fract_nonfunctional), c = utils.rgb_red_antibiotic(phi_r_0_i_idx))


    ax.set_xscale('log', basex=10)
    #ax.set_yscale('log', basey=10)
    ax.legend(loc='upper right', title='% non-functional ribosomes', fontsize=8)

    ax.set_xlabel('Initial phage concentration (PFU/mL)', fontsize=10)
    ax.set_ylabel('Biomass at lysis (OD600)', fontsize=12)


    fig.subplots_adjust(hspace=0.35, wspace=0.40)
    fig_name = "%stest_theory.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def calculate_intercept(lambda_cell, lambda_phage, I):

    return n_cell_0 * numpy.log(1 + (lambda_cell/(I*tau_lysis*lambda_phage*n_cell_0)) * numpy.log(lambda_cell/I))
    

@plt.FuncFormatter
def fake_log(x, pos):
    'The two args are the value and tick position'
    return r'$10^{%d}$' % (  x )


@plt.FuncFormatter
def fake_log_float(x, pos):
    'The two args are the value and tick position'

    return r'$10^{%.1f}$' % (  x )

    
def plot_intercept_heatmap():

    lambda_cell_all = numpy.linspace(0.1, 2, num=10, endpoint=True)
    # as burst size
    burst_size_all = numpy.logspace(1, 3, num=10, base=10, endpoint=True)[::-1]

    intecept_matrix = [ [calculate_intercept(lambda_cell_i, burst_size_j/tau_lysis, 1e-11) for lambda_cell_i in lambda_cell_all] for burst_size_j in burst_size_all ]
    intecept_matrix = numpy.asarray(intecept_matrix)
    intecept_matrix_od = intecept_matrix/od600_to_cell_density

    min_, max_ = numpy.min(intecept_matrix_od), numpy.max(intecept_matrix_od)

    fig = plt.figure(figsize = (8, 8))
    fig.subplots_adjust(bottom= 0.15)

    ax = plt.subplot2grid((1, 1), (0, 0))
    #ax_im = ax.imshow(intecept_matrix_od, interpolation="nearest", cmap='Blues',  vmin=min_, vmax=max_)

    burst_size_all_log10 = numpy.log10(burst_size_all)

    pcm = ax.pcolor(lambda_cell_all, burst_size_all_log10, intecept_matrix_od,
                       norm=colors.Normalize(vmin=min_, vmax=max_),
                       cmap='Blues')

    #ax.set_xticks(list(range(len(lambda_cell_all)))[::2])
    #ax.set_xticklabels([round(t, 1) for t in lambda_cell_all[::2]])

    #ax.set_yticks(list(range(len(burst_size_all)))[::2])
    #ax.set_yticklabels([round(t, 1) for t in numpy.log10(burst_size_all)[::2]])
    

    ax.yaxis.set_major_formatter(fake_log_float)
    #ax.yaxis.set_major_formatter(fake_log)

    #cb_ax = fig.add_axes([.82,.16,.03,.70])
    cb_ax = fig.add_axes([.92,.16,.03,.70])
    fig.colorbar(pcm, orientation='vertical',cax=cb_ax)
    cb_ax.tick_params(labelsize=7)
    cb_ax.set_ylabel('Intercept of ' + r'$N_{\mathrm{phage}}(0)$' + ' vs. ' + r'$N_{\mathrm{cell}}^{\mathrm{max}}$', rotation=270, labelpad=25, fontsize=17)

    #print(min_, max_)

    ax.set_xlabel('Cell growth rate, ' + r'$\lambda_{\mathrm{cell}}$', fontsize=20)
    #ax.set_ylabel('Burst size, ' + r'$\lambda_{\mathrm{phage}} \cdot \tau_{\mathrm{lysis}}$', fontsize=20)
    ax.set_ylabel('Burst size, ' + r'$\lambda_{\mathrm{phage}} $', fontsize=20)

    fig.subplots_adjust(hspace=0.3, wspace=0.25)
    fig_name = "%sintercept_heatmap.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
    plt.close()




plot_change_in_cell_growth_only()

#plot_intercept_heatmap()

#plot_change_in_cell_growth_only()

