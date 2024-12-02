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

hdf5_file = 'output_nu_tau_4km_ct18nlo_allm_stochastic_1e6.h5'
real_file = 'output_nu_tau_4km_ct18nlo_allm_stochastic_1e8.h5'


nu_type = 'neutrino'
ch_lepton = 'tau'
cross_section_model = 'ct18nlo'
pn_model = 'allm'
prop_type = 'stochastic'
idepth = 4
stats_small = 1e5
stats_big = 1e8

energies = [9.0]
energies = np.array(energies)
energies = 10**energies

angles_lto = np.arange(0.1,1,0.1).astype('float')
angles_gto = np.arange(1,43).astype('float')
angles_full = np.concatenate((angles_lto,angles_gto), axis = None)
angles = np.arange(1,10).astype('float')
angles_short = np.concatenate((angles_lto,angles), axis = None)


fig,axs =plt.subplots(1,figsize=(10/1.1, 8/1.1))

pxt_nr_big = data.get_pexit(nu_type,ch_lepton,energies[0],idepth,cross_section_model, pn_model,prop_type,stats_big)[2]
pxt_wr_big = data.get_pexit(nu_type,ch_lepton,energies[0],idepth,cross_section_model, pn_model,prop_type,stats_big)[1]


pxt_nr_small = data.get_pexit(nu_type,ch_lepton,energies[0],idepth,cross_section_model, pn_model,prop_type,stats_small)[2]
pxt_wr_small = data.get_pexit(nu_type,ch_lepton,energies[0],idepth,cross_section_model, pn_model,prop_type,stats_small)[1]

print(pxt_wr_big,pxt_wr_small)


axs.loglog(angles_full,pxt_nr_big, linestyle = '--', label = 'Fortran')
axs.loglog(angles_short,pxt_nr_small, label = 'Python')
ticks = (1,2,3,4,5,10,15,20,30,40)
axs.set_xticks(ticks)
axs.set_xticklabels([i for i in ticks])
axs.grid(which='major', axis='both', linestyle='--')
axs.grid(which='minor', linestyle=':', linewidth='0.2', color='black')
axs.xaxis.set_ticks_position('both')
for axis in [axs.xaxis]:
    axis.set_major_formatter(ScalarFormatter())
    
plt.title('Pexit no regen compare python (1e5) vs old (1e8)')
axs.set_xlabel(r"$\beta_{tr}$ [degrees]")
axs.set_ylabel(r'$P^{(\tau)}_{exit}$')
plt.legend()
plt.savefig('python_fortran_compare_nr.pdf',bbox_inches = 'tight', dpi = 400)

plt.clf()

fig,axs =plt.subplots(1,figsize=(10/1.1, 8/1.1))

axs.loglog(angles_full,pxt_wr_big, linestyle = '--', label = 'Fortran')
axs.loglog(angles_short,pxt_wr_small, label = 'Python')
ticks = (1,2,3,4,5,10,15,20,30,40)
axs.set_xticks(ticks)
axs.set_xticklabels([i for i in ticks])
axs.grid(which='major', axis='both', linestyle='--')
axs.grid(which='minor', linestyle=':', linewidth='0.2', color='black')
axs.xaxis.set_ticks_position('both')
for axis in [axs.xaxis]:
    axis.set_major_formatter(ScalarFormatter())
plt.title('Pexit regen compare python (1e5) vs old (1e8)')
axs.set_xlabel(r"$\beta_{tr}$ [degrees]")
axs.set_ylabel(r'$P^{(\tau)}_{exit}$')
plt.legend()
plt.savefig('python_fortran_compare_wr.pdf',bbox_inches = 'tight', dpi = 400)