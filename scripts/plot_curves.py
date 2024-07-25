import config
import matplotlib as mpl
import matplotlib.pyplot as plt

import utils

data_dict = utils.parse_data()

fig = plt.figure(figsize = (4*len(utils.phage_labels), 4*len(utils.antibiotic_labels)))
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)



value_max_all = []
for key, value in data_dict.items():
    if key != 'hours':
        value_max_all.append(max(value))







for p_idx, p in enumerate(utils.phage_labels):
        
    for a_idx, a in enumerate(utils.antibiotic_labels):

        ax = plt.subplot2grid((len(utils.antibiotic_labels), len(utils.phage_labels)), (a_idx, p_idx))

        labels_to_plot = [s for s in data_dict.keys() if (p in s) and (a in s)]

        for l in labels_to_plot:
            ax.plot(data_dict['hours'], data_dict[l], c='k', lw=2)

        if a_idx == 0:
            ax.set_title(('Log10 PFU/mL = %s' % p[-1]), fontsize=12, weight='bold')

        if a == 'cN':
            label_a = 'Chlor. conc. ' + r'$\mu g/ mL$' + ' = 0'
        else:
            label_a = 'Chlor. conc. ' + r'$\mu g/mL = 25 \cdot 2^{{{}}}  $'.format('-' + a[1:])


        if p_idx == 0:
            ax.text(-0.3, 0.5, label_a, rotation=90, fontsize=12, ha='center', va='center', transform=ax.transAxes)
            # fontweight='bold',


        ax.set_ylim((0, max(value_max_all)))


        if a_idx == len(utils.antibiotic_labels) - 1:
            ax.set_xlabel('Time (h)', fontsize=12)

        if p_idx == 0:
            ax.set_ylabel('Biomass (OD600)', fontsize=12)



fig.subplots_adjust(hspace=0.25, wspace=0.25)
fig_name = "%sgrowth_cuves.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
plt.close()