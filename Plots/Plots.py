from __future__ import print_function
import glob
import os
import numpy as np
import importlib_resources
import pathlib
import shutil
import h5py

# Load nuPyProp data & models modules
from nupyprop import models
from nupyprop import data
# Load astropy modules
from astropy.table import Table
from astropy.io import ascii

# Plotting modules
import matplotlib as mpl
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
from itertools import cycle
import seaborn as sns
import itertools



sns.set_theme()
sns.set_context("paper")
sns.set_style("white", {"font.family": "STIXGeneral"})

# Plotting defaults
mpl.rcParams['xtick.labelsize'] = 26
mpl.rcParams['ytick.labelsize'] = 26
mpl.rcParams['axes.labelsize'] = 30
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.minor.width'] = 1.5
mpl.rcParams['xtick.major.pad'] = 8
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['ytick.minor.size'] = 5
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.width'] = 1.5
mpl.rcParams['ytick.major.pad'] = 8
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['legend.fontsize'] = 23
plt.rcParams['legend.title_fontsize'] = 23
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['axes.titlesize'] = 30
mpl.rcParams['figure.titlesize'] = 30
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['axes.linewidth'] = 2


#hdf5_file = 'output_nu_tau_4km_ct18nlo_bdhm_stochastic-FINAL_1e8.h5' # select your output file
hdf5_file = 'output_nu_tau_4km_ct18nlo_allm_stochastic_1e6.h5'
with h5py.File(hdf5_file, 'r') as hf:
    for attr in hf.attrs:
        print(str(attr) + " = " + str(hf.attrs[attr])) # print the file attributes
        
nu_type = 'neutrino' # nu for neutrino & anu for anti-neutrino
ch_lepton = 'tau' # type of charged lepton
cross_section_model = 'ct18nlo' # neutrino cross-section model
pn_model = 'allm' # photonuclear energy loss model
#pn_model = 'bdhm'
idepth = 4 # depth of water layer in km
stats = 1e8 # no. of ingoing neutrinos ie., statistics (also mentioned in the file name)
prop_type = 'stochastic' # type of energy loss; can be stochastic or continuous
prop_types = 'stochastic-old'
#energies = [6.0,6.25,6.5,6.75,7.0,7.25,7.5,7.75,
#            8.0,8.25,8.5,8.75,9.0,9.25,9.5,9.75,10.0,10.25,10.5,10.75,11.0,11.25,
#            11.5,11.75,12.0] 
energies = [7.0,8.0,9.0]
energies = np.array(energies)
energies = 10**energies# default energies in nuPyProp
angles1 = np.arange(0.1,1,0.1).astype('float64') # default angles in nuPyProp
angles2 = np.arange(1,43).astype('float64') 
angles3 = np.concatenate((angles1,angles2), axis=None)
#angles3 = np.asarray([1,2,3,4,5,10])
angles = np.arange(1,43).astype('float64')
fig, axs = plt.subplots(1, figsize=(10/1.1, 8/1.1))

energy_grp = [] # so that 7, 7.25, etc. are grouped into the same legend

lines = ["--","-.",":",(0,(3,1,1,1))]
linecycler = cycle(lines)

# colors = ['r','g','b','k','m']
# color_map = cycle(colors)
color_map = cycle(sns.color_palette("colorblind"))

for i,energy in enumerate(energies):
    print(i,float(np.log10(energy)))
    pexit_no_regen = data.get_pexit(nu_type, ch_lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats)[2]
    #pexit_no_regen_old = data.get_pexit(nu_type, ch_lepton, energy, idepth, cross_section_model, pn_model, prop_types, stats)[2]# without regen
    pexit_regen = data.get_pexit(nu_type, ch_lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats)[1] # with regen
    #pexit_regen_old = data.get_pexit(nu_type, ch_lepton, energy, idepth, cross_section_model, pn_model, prop_types, stats)[1]

    pexit_regen_ratio = pexit_regen[8:len(pexit_regen)]
    ratio_list = []
    #ratio = np.zeros_like(pexit_regen_old)
    '''for i in range(len(pexit_regen_ratio) -1):
        
        if pexit_regen_old[i] == 0.0:
            print('if')
            ratio[i] = 0
            print(ratio)
        else:
            print('else')
            ratio[i]= (pexit_regen_ratio[i]/pexit_regen_old[i])
        '''

    energy_log = float(np.log10(energy))

    if int(np.log10(energy)) not in energy_grp:
        color = next(color_map)
        energy_grp.append(int(np.log10(energy)))

        if energy_log%2 == 0 or (energy_log+1)%2 == 0:
            axs.loglog(angles3, pexit_regen, c = color, label = '%.0d' % energy_log)
            #axs.loglog(angles, pexit_regen_old, c = color,ls = '--', label = '%.0d' % energy_log)
            #axs.loglog(angles3, pexit_no_regen, c = color, label = '%.0d' % energy_log)
            #axs.loglog(angles, pexit_no_regen_old, c = color,ls = '--', label = '%.0d' % energy_log)
        else:
            axs.loglog(angles3, pexit_regen, c = color, label = '%.2f' % energy_log)
            #axs.loglog(angles, pexit_regen_old, c = color,ls = '--', label = '%.0d' % energy_log)
            #axs.loglog(angles3, pexit_no_regen, c = color, label = '%.2f' % energy_log)
            #axs.loglog(angles, pexit_no_regen_old, c = color,ls = '--', label = '%.0d' % energy_log)
            
    else:
        axs.loglog(angles3, pexit_regen, ls = next(linecycler), color = color)
        #axs.loglog(angles, pexit_regen_old, c = color,ls = '--', label = '%.0d' % energy_log)
        #axs.loglog(angles3, pexit_no_regen, c = color, label = '%.2f' % energy_log)
        #axs.loglog(angles, pexit_no_regen_old, c = color,ls = '--', label = '%.0d' % energy_log)


#axs.legend(loc='best', ncol=2, title = r'log$_{10}(E_{\nu}/GeV)$', frameon=False, framealpha=0.5, fontsize=23)

ticks = (1,2,3,4,5,10,15,20,30,40)
axs.set_xticks(ticks)
axs.set_xticklabels([i for i in ticks])
axs.grid(which='major', axis='both', linestyle='--')
axs.grid(which='minor', linestyle=':', linewidth='0.2', color='black')
axs.xaxis.set_ticks_position('both')
axs.yaxis.set_ticks_position('both')
axs.tick_params(axis='x', which='both', labelbottom = True, labeltop = False)
axs.tick_params(axis='y', which='both', left = True, labelleft = True, labelright= True)
for axis in [axs.xaxis]:
    axis.set_major_formatter(ScalarFormatter())
    
axs.set_xlabel(r"$\beta_{tr}$ [degrees]")
if ch_lepton=='tau':axs.set_ylabel(r'$P^{(\tau)}_{exit}$')
if ch_lepton=='muon':axs.set_ylabel(r'$P^{(\mu)}_{exit}$')
#axs.set_title(r"$P_{exit}$ Vs. Earth Emergence Angles for %ss" % str.capitalize(ch_lepton), fontsize=27)
axs.set_ylim(1e-6,1e-1)
axs.set_xlim(1,42)
plt.tight_layout()
plt.savefig('pexit-rainbow-allm_1e8_new_noreg.pdf', bbox_inches='tight', dpi=400)
#plt.show()


print(len(ratio_list))
for i in range(len(ratio_list)):
    for j in range(len(ratio_list[i])):
        
        plt.plot(angles,ratio_list[i][j])
    
plt.show()
