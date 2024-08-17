import config
import numpy
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as stats

import utils

data_dict = utils.parse_data()




def plot_log_od():

    fig = plt.figure(figsize = (4*len(utils.phage_labels), 4*len(utils.antibiotic_labels)))
    fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

    value_max_all = []
    value_min_all = []
    for key, value in data_dict.items():
        if key != 'hours':
            value_max_all.append(max(value))
            value_min_all.append(min(value))

    for p_idx, p in enumerate(utils.phage_labels):
            
        for a_idx, a in enumerate(utils.antibiotic_labels):

            ax = plt.subplot2grid((len(utils.antibiotic_labels), len(utils.phage_labels)), (a_idx, p_idx))

            #labels_to_plot = [s for s in data_dict.keys() if (p in s) and (a in s)]
            labels_to_plot = [s for s in data_dict.keys() if (p in s) and (a in s) ]

            for l in labels_to_plot:
                ax.plot(data_dict['hours'], data_dict[l], c='k', lw=2, alpha=0.8)

            if a_idx == 0:
                ax.set_title(('Log10 PFU/mL = %d' % (int(p[-1]) -2 )), fontsize=12, weight='bold')

            if a == 'cN':
                label_a = 'Chlor. conc. ' + r'$\mu g/ mL$' + ' = 0'
            else:
                label_a = 'Chlor. conc. ' + r'$\mu g/mL = 25 \cdot 2^{{{}}}  $'.format('-' + a[1:])


            if p_idx == 0:
                ax.text(-0.3, 0.5, label_a, rotation=90, fontsize=12, ha='center', va='center', transform=ax.transAxes)
                # fontweight='bold',

            ax.set_ylim((min(value_min_all), max(value_max_all)))
            ax.set_yscale('log', basey=10)


            if a_idx == len(utils.antibiotic_labels) - 1:
                ax.set_xlabel('Time (h)', fontsize=12)

            if p_idx == 0:
                ax.set_ylabel('Biomass (OD600)', fontsize=12)



    fig.subplots_adjust(hspace=0.25, wspace=0.25)
    fig_name = "%sgrowth_cuves.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
    plt.close()



def plot_deriv_log_od(n=1):

    hours = data_dict['hours']
                
    # get limits
    deriv_all = []
    deriv_dict = {}
    for key, value in data_dict.items():
        if key != 'hours':

            od = data_dict[key]
            log_od = numpy.log(od)
            delta_log_od = utils.fd_derivative(log_od, hours, n=n, m=3)

            clean_idx = ~numpy.isnan(delta_log_od)
            delta_log_od_clean = delta_log_od[clean_idx]
            hours_clean = hours[clean_idx]

            deriv_dict[key] = {}
            deriv_dict[key]['delta_log_od'] = delta_log_od_clean
            deriv_dict[key]['hours'] = hours_clean

            deriv_all.extend(delta_log_od_clean.tolist())

    

    fig = plt.figure(figsize = (4*len(utils.phage_labels), 4*len(utils.antibiotic_labels)))
    fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


    for p_idx, p in enumerate(utils.phage_labels):
            
        for a_idx, a in enumerate(utils.antibiotic_labels):

            ax = plt.subplot2grid((len(utils.antibiotic_labels), len(utils.phage_labels)), (a_idx, p_idx))

            labels_to_plot = [s for s in data_dict.keys() if (p in s) and (a in s)]

            for l in labels_to_plot:
                ax.plot(deriv_dict[l]['hours'], deriv_dict[l]['delta_log_od'], c='dodgerblue', alpha=0.8, lw=2)

            if a_idx == 0:
                ax.set_title(('Log10 PFU/mL = %s' % p[-1]), fontsize=12, weight='bold')

            if a == 'cN':
                label_a = 'Chlor. conc. ' + r'$\mu g/ mL$' + ' = 0'
            else:
                label_a = 'Chlor. conc. ' + r'$\mu g/mL = 25 \cdot 2^{{{}}}  $'.format('-' + a[1:])


            if p_idx == 0:
                ax.text(-0.3, 0.5, label_a, rotation=90, fontsize=12, ha='center', va='center', transform=ax.transAxes)
                # fontweight='bold',

            ax.set_ylim((-1, max(deriv_all)))
            ax.axhline(y=0, ls=':', lw=1, label='No change')

            if a_idx == len(utils.antibiotic_labels) - 1:
                ax.set_xlabel('Time (h)', fontsize=12)

            if n == 1:
                label = '1st order derivative of log OD'

            else:
                label = '2nd order derivative of log OD'

            if p_idx == 0:
                ax.set_ylabel(label, fontsize=12)



    fig.subplots_adjust(hspace=0.25, wspace=0.25)
    fig_name = "%sderiv_log_od_order_%s.png" % (config.analysis_directory, n)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
    plt.close()



def plot_phage_vs_lysis_od():

    fig = plt.figure(figsize = (4.3*len(utils.antibiotic_labels), 4))
    fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

    # load dictionary
    lysis_dict = {}
    for line in open('%ssample_dip_status.csv' % config.data_directory,  encoding='utf-8-sig'):
        line = line.strip().split(',')

        if line[1] == 'TRUE':
            lysis_dict[line[0]] = float(line[-1])

    linregress_dict = {}
    for a_idx, a in enumerate(utils.antibiotic_labels):

        ax = plt.subplot2grid((1, len(utils.antibiotic_labels)), (0, a_idx))

        labels_to_plot = [s for s in lysis_dict.keys() if (a in s)]
        #labels_to_plot = [s for s in lysis_dict.keys() if (a in s) and ('r1' in s)]

        #lysis_od = numpy.asarray([max(data_dict[l][data_dict['hours'] <= lysis_dict[l] ])   for l in labels_to_plot])
        #lysis_od = numpy.asarray([max(data_dict[l][data_dict['hours'] <= lysis_dict[l] ]) - data_dict[l][0]  for l in labels_to_plot])
        lysis_od = numpy.asarray([max(data_dict[l][data_dict['hours'] <= lysis_dict[l] ])  for l in labels_to_plot])
        time_lysis_od = numpy.asarray([data_dict['hours'][(data_dict['hours'] <= lysis_dict[l])][numpy.argmax(data_dict[l][data_dict['hours'] <= lysis_dict[l]])] for l in labels_to_plot])

        initial_od_log10 = numpy.log10([data_dict[l][0] for l in labels_to_plot])
        mean_initial_od_log10 = numpy.mean(initial_od_log10)
        se_initial_od_log10 = numpy.std(initial_od_log10)/numpy.sqrt(len(initial_od_log10))

        phage_conc = [10**(int(s[1])-2) for s in labels_to_plot]

        ax.scatter(phage_conc, lysis_od, s=20, c=utils.rgb_red_antibiotic(a_idx))
        #ax.scatter(phage_conc, time_lysis_od, s=20)
        
        # try regession
        slope, intercept, r_valuer_value, p_value, std_err = stats.linregress(numpy.log10(phage_conc), numpy.log10(lysis_od))

        #std_error_intercept

        #print(len(stats.linregress(numpy.log10(phage_conc), numpy.log10(lysis_od))))

        linregress_dict[a] = {}
        linregress_dict[a]['slope'] = slope
        linregress_dict[a]['intercept'] = intercept

        # https://github.com/scipy/scipy/blob/a3ffdface6d8779ffd91f605e4e102a9fda65a7f/scipy/stats/_stats_py.py
        # Line 10366
        ssxm, ssxym, _, ssym = numpy.cov(numpy.log10(phage_conc), numpy.log10(lysis_od), bias=1).flat
        intercept_stderr = std_err * numpy.sqrt(ssxm + numpy.mean(numpy.log10(phage_conc))**2)

        linregress_dict[a]['slope_se'] = std_err
        linregress_dict[a]['intercept_se'] = intercept_stderr


        x_log10_range =  numpy.linspace(min(numpy.log10(phage_conc)) , max(numpy.log10(phage_conc)) , 10000)
        y_log10_fit_range = 10 ** (slope*x_log10_range + intercept)
        #y_log10_null_range = 10 ** (utils.slope_null*x_log10_range + intercept)

        #print(y_log10_fit_range)

        ax.plot(10**x_log10_range, y_log10_fit_range, c='k', lw=2.5, linestyle='--', zorder=2, label="OLS regression slope")
        
        ax.text(0.7, 0.9, 'Slope = %0.3f' % slope, fontsize=10, ha='center', va='center', transform=ax.transAxes)
        ax.text(0.7, 0.8, 'Intercept = %0.3f' % intercept, fontsize=10, ha='center', va='center', transform=ax.transAxes)


        ax.axhline(y=10**mean_initial_od_log10, lw=2, ls=':', c='k', label='Mean log10 initial OD', zorder=2)
        ax.axhspan(10**(mean_initial_od_log10 - se_initial_od_log10), 10**(mean_initial_od_log10 + se_initial_od_log10), facecolor='grey', alpha=0.8, zorder=1)



        ax.set_xlim((0.7e0, 1.3e5))
        #ax.set_ylim((0, 0.9))
        ax.set_ylim((8e-2, 2e0))

        ax.set_xscale('log', basex=10)
        ax.set_yscale('log', basey=10)

        if a == 'cN':
            label_a = 'Chlor. conc. ' + r'$\mu g/ mL$' + ' = 0'
        else:
            label_a = 'Chlor. conc. ' + r'$\mu g/mL = 25 \cdot 2^{{{}}}  $'.format('-' + a[1:])

        ax.set_title(label_a, fontsize=12)
        ax.set_xlabel('Initial phage concentration (PFU/mL)', fontsize=10)
        ax.set_ylabel('Biomass at lysis (OD600)', fontsize=12)
        #ax.set_ylim((0.1, 1.1))



    fig.subplots_adjust(hspace=0.25, wspace=0.25)
    fig_name = "%sphage_vs_lysis_od.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
    plt.close()


    # plot chlpramphenicol vs. slope
    fig = plt.figure(figsize = (8.5, 4))
    fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

    ax_intercept = plt.subplot2grid((1, 2), (0, 0))
    ax_slope = plt.subplot2grid((1, 2), (0, 1))

    intercept_all = [linregress_dict[a]['intercept'] for a in utils.antibiotic_labels]
    intercept_se_all = [linregress_dict[a]['intercept_se'] for a in utils.antibiotic_labels]

    slope_all = [linregress_dict[a]['slope'] for a in utils.antibiotic_labels]
    slope_se_all = [linregress_dict[a]['slope_se'] for a in utils.antibiotic_labels]

    print(intercept_all)


    colors = [utils.rgb_red_antibiotic(i) for i in range(len(utils.antibiotic_labels))]

    ax_intercept.errorbar(utils.antibiotic_conc, intercept_all, yerr=intercept_se_all, linestyle='-', marker='o', c='k', elinewidth=2, alpha=1, zorder=1)
    #ax_intercept.scatter(utils.antibiotic_conc, intercept_all, color=colors, zorder=2)
    ax_intercept.scatter(utils.antibiotic_conc, intercept_all, s=30, facecolors=colors, edgecolors='k', linewidth=1, alpha=0.9, zorder=3)


    ax_slope.errorbar(utils.antibiotic_conc, slope_all, yerr=slope_se_all, linestyle='-', marker='o', c='k', elinewidth=2, alpha=1, zorder=1)
    #ax_slope.scatter(utils.antibiotic_conc, slope_all, zorder=2, c)
    ax_slope.scatter(utils.antibiotic_conc, slope_all, s=30, facecolors=colors, edgecolors='k', linewidth=1, alpha=0.9, zorder=3)



    ax_intercept.set_xlabel('Chloramphenicol conc., ' + r'$\mu \mathrm{g}/ \mathrm{mL}$',  fontsize=12)
    #ax_intercept.set_ylabel("Intercept of log-log lysis OD vs. PFU/mL regression")
    ax_intercept.set_ylabel("Intercept", fontsize=12)


    ax_slope.set_xlabel('Chloramphenicol conc., ' + r'$\mu \mathrm{g}/ \mathrm{mL}$',  fontsize=12)
    #ax_slope.set_ylabel("Slope of log-log lysis OD vs. PFU/mL regression")
    ax_slope.set_ylabel("Slope",  fontsize=12)


    fig.subplots_adjust(hspace=0.3, wspace=0.35)
    fig_name = "%sslope_intercept.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
    plt.close()




    # plot slope


plot_phage_vs_lysis_od()

#plot_log_od()
#plot_deriv_log_od()
