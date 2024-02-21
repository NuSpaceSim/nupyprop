#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 10:27:38 2022

@author: dikgarg
"""
import numpy as np
import pandas as pd

#Load nuPyProp data & models modules
from nupyprop import data
from nupyprop import models
import nupyprop.geometry as Geometry

#Load astropy modules
from astropy.table import Table
from astropy.io import ascii

# Plotting modules
import matplotlib as mpl
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
from matplotlib import pyplot
from itertools import cycle
import seaborn as sns
import itertools
import matplotlib.gridspec as gridspec
import random
from scipy import interpolate

nu_type = 'neutrino' # nu for neutrino & anu for anti-neutrino
ch_lepton = 'tau' # type of charged lepton
ch_leptons = 'taus' # type of charged lepton
cross_section_model = 'ct18nlo' # neutrino cross-section model
pn_model = 'allm' # photonuclear energy loss model
pn_models = 'bdhm'
idepth = 4 # depth of water layer in km
idepths1 = 1
idepths2 = 0.5
idepths3 = 0.1
idepths4 = 0
stats = 1e8 # no. of ingoing neutrinos ie., statistics, Diksha code
stats1 = 1e5
stats2 = 1e6
stats3 = 1e7
stats4 = 1e9
prop_type = 'stochastic' # type of energy loss; can be stochastic or continuous

mpl.style.use('/Users/dikgarg/Desktop/Research/Neutrinos/paper.mplstyle')

files_path = "/Users/dikgarg/Desktop/Research/Neutrinos/Tau_polarization/nupyprop/"
#energy = [8,10]

#data.process_htc_out(files_path, nu_type, ch_lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats, elep_mode=True, arg=None)

def data_extract(str):
    data = pd.read_csv(str,sep='\s+',header=None)
    #data = pd.DataFrame(data)
    return data;

#===============================================#
#========== to make plots from h5 ==============#
#===============================================#
# plt.rcParams['xtick.major.pad']='5'
# plt.rcParams['xtick.labelsize']='24'
# plt.rcParams['ytick.labelsize']='24'

energies = [9] #np.arange(7,12,1)
angles1 = [0.2,0.3,0.6,0.7,0.8,0.9,1,2,3,4,5,10,15,20] #np.arange(1,43).astype('float64') # default angles in nuPyProp
angles = np.arange(0.1,1,0.1)
angles = np.append(angles, np.arange(1,43).astype('float64')) # default angles in nuPyProp
angles2 = [1,2,3,4,5,10,15,20]
# print(angles)

fig= plt.figure(figsize=[8,6], constrained_layout=True)
spec=gridspec.GridSpec(ncols=1,nrows=2,figure=fig)  # set up subplot grid

plt1 = fig.add_subplot(spec[0:2,0])

lss = ['-','--']
colors = sns.color_palette('colorblind')[:len(energies)]

for e, c in zip(energies, colors):
    pexit_regen1 = data.get_pexit(nu_type, ch_lepton, 10**e, idepth, cross_section_model, pn_model, prop_type, stats)[1] # without regen
    pexit_regen2 = data.get_pexit(nu_type, ch_lepton, 10**e, idepths1, cross_section_model, pn_model, prop_type, stats3)[1] # with regen
    pexit_regen3 = data.get_pexit(nu_type, ch_lepton, 10**e, idepths2, cross_section_model, pn_model, prop_type, stats3)[1] # with regen
    pexit_regen4 = data.get_pexit(nu_type, ch_lepton, 10**e, idepths3, cross_section_model, pn_model, prop_type, stats)[1] # with regen
    pexit_regen5 = data.get_pexit(nu_type, ch_lepton, 10**e, idepths4, cross_section_model, pn_model, prop_type, stats3)[1] # with regen

    plt1.loglog(angles, pexit_regen1, linestyle=lss[0], marker='o', color = 'black', label='4 km')
    plt1.loglog(angles1, pexit_regen2, linestyle=lss[0], marker='o', color = 'orange', label='1 km')
    plt1.loglog(angles1, pexit_regen3, linestyle='dashdot', marker='o', color = 'green', label='0.5 km')
    plt1.loglog(angles2, pexit_regen4, linestyle=lss[1], marker='o', color = 'red', label='0.1 km')
    plt1.loglog(angles1, pexit_regen5, linestyle='dotted', marker='o', color = 'blue', label='0 km')

# # #plt.text(8.8, 0.07, r'$\log_{10}(E_{\nu}$/GeV)',fontsize=18)
# # m = plt1.plot( [], [], label= r'$\log_{10}(E_{\nu}$/GeV)', ls=lss[0])
# # j = [plt1.plot( [],[], color=colors[i], label= str(energies[i]), ls=lss[0])[0] for i in range(len(energies))]
# # first_legend = plt1.legend(handles=j, frameon=False, loc='best', ncol=3, fontsize=19) # (0.38,0.65) (25.08, 0.015))
# # ax = plt.gca().add_artist(first_legend)

# # titles = ['Old', 'New'] #['Without regeneration', 'With regeneration']
# # #titles = ['ALLM', 'BDHM']
# # linestyles = ['-', '--']

# # h = [plt1.plot( [],[], color='black', label=titles[i], ls=linestyles[i])[0] for i in range(len(titles)) ]
# # plt1.legend(handles=(h), frameon=False, loc='lower left', fontsize=18)

# plt1.set_xlim(0.2,21)
# plt1.set_ylim(1*1e-6, 2*1e-2)
# plt1.set_ylabel(r'$P^{(\tau)}_{exit}$', fontsize=24)
# plt1.grid(which='major', axis='both', linestyle=':', linewidth=0.5)
# plt1.grid(which='minor', axis='both', linestyle=':', linewidth=0.5)
# plt1.set_xticklabels([])

# plt2 = fig.add_subplot(spec[2:3,0])

# for e, c in zip(energies, colors):
#     pexit_no_regen = data.get_pexit(nu_type, ch_lepton, 10**e, idepth, cross_section_model, pn_model, prop_type, stats)[1] # without regen
#     pexit_regen = data.get_pexit(nu_type, ch_lepton, 10**e, idepth, cross_section_model, pn_model, prop_type, stats)[2] # with regen

#     plt2.semilogx(angles, pexit_no_regen[9:]/pexit_regen[9:], linestyle=lss[1], color = c)

# plt2.set_xlim(1,42)
# ticks = (1,2,3,4,5,10,15,20,30,40)
# plt2.set_xticks(ticks)
# plt2.set_xticklabels(ticks)
# plt2.set_ylabel('Ratio', fontsize=24)
# #plt2.set_ylabel('Ratio (ALLM/BDHM)', fontsize=19)
#plt1.set_xlabel(r"$\beta_{tr}$ [degrees]", fontsize=24)
# plt2.grid(which='major', axis='both', linestyle=':', linewidth=0.5)
# plt2.grid(which='minor', axis='both', linestyle=':', linewidth=0.5)

plt.legend()
#plt.savefig('pexit_waterdepth8.png', bbox_inches='tight', dpi=400)
plt.show()

#===============================================#
#============ to modify h5files ================#
#===============================================#
#data.add_cdf(nu_type, ch_lepton, idepth, cross_section_model, pn_model, prop_type, stats, bins=np.logspace(-7,0,71))
# energies = [6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0]
# for energy in energies:
#     data.get_pexit(nu_type, ch_lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats, out=False, out_h5=True, arg=None)
#     data.get_cdf(nu_type, ch_lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats, out=False, out_h5=True)
    # data.get_pexit(nu_type, ch_lepton, ch_leptons, energy, idepth, cross_section_model, pn_model, prop_type, stats, stats1, stats2, stats3, stats4, out=False, out_h5=True, arg=None)

    # data.get_clep_out(nu_type, ch_lepton, ch_leptons, energy, idepth, cross_section_model, pn_model, prop_type, stats, stats1, stats2, stats3, stats4, out=False, out_h5=False, arg=None)

#===============================================#
#=========== to check polarization =============#
#===============================================#
# energies=[1e10]
# angle = data.get_polarization(nu_type, ch_lepton, energies[0], idepth, cross_section_model, pn_models, prop_type, stats)[0]
# #pol = data.get_polarization(nu_type, ch_lepton, energies[0], idepth, cross_section_model, pn_model, prop_type, stats)[1]
# pol1 = data.get_polarization(nu_type, ch_lepton, energies[0], idepth, cross_section_model, pn_models, prop_type, stats)[1]

# #plt.plot(angle, -1*pol, label='old pol')
# plt.plot(angle, -1*pol1, label = 'new pol')
# plt.legend()
# plt.show()

#===============================================#
#=========== to check cdf ======================#
#===============================================#
#energies = [1e8] #, 11.0, 12.0]

# z = data.get_cdf(nu_type, ch_lepton, energies[0], idepth, cross_section_model, pn_model, prop_type, stats)[0] # get z bins
# angles_in_data = data.get_cdf(nu_type, ch_lepton, energies[0], idepth, cross_section_model, pn_model, prop_type, stats)[1] # get the list of earth emergence angles in the output file
# angle_indices = [int(np.where(angles_in_data == i)[0]) for i in angles_in_data] # gets the indices of the angles we need based on that in the output file
# e_out_cdf = data.get_cdf(nu_type, ch_lepton, energies[0], idepth, cross_section_model, pn_model, prop_type, stats)[2] # binned CDF values for all angles

# z1 = data.get_cdf(nu_type, ch_lepton, energies[0], idepth, cross_section_model, pn_models, prop_type, stats)[0] # get z bins
# angles_in_data1 = data.get_cdf(nu_type, ch_lepton, energies[0], idepth, cross_section_model, pn_models, prop_type, stats)[1]
# angle_indices1 = [int(np.where(angles_in_data1 == i)[0]) for i in angles_in_data1]
# e_out_cdf1 = data.get_cdf(nu_type, ch_lepton, energies[0], idepth, cross_section_model, pn_models, prop_type, stats)[2] # binned CDF values for all angles

# an = [0, 9, 19, 34] #1, 10, 20, 30, 40 degrees
# c = sns.color_palette('colorblind')[:len(an)] #['blue', 'black', 'green', 'orange', 'red']
# #ls = ['-', '--']
# ls = ['-', '--', ':', '-.']

# for i in range(len(an)):
#     print("eout = ", 1e8*e_out_cdf[an[i]+9])
#     plt.semilogx(z, e_out_cdf[an[i]+9], color=c[i], linestyle=ls[i]) #, label=str(an[i]+1))
#     #plt.semilogx(z1, e_out_cdf1[an[i]], color=c[i], linestyle=ls[1])

# j = [plt.plot( [],[], color=c[i], label=str(an[i]+1)+r"$^\circ$", ls=ls[0])[0] for i in range(len(an))]
# first_legend = plt.legend(handles=j, title=r"$\beta=$", frameon=False, loc='best')
# ax = plt.gca().add_artist(first_legend)

# labels = ['New CDF', 'Old CDF']
# linestyles = ['-', '--']

# h = [plt.plot( [],[], color='black', label=labels[i], ls=linestyles[i])[0] for i in range(len(labels)) ]
# first_legend1 = plt.legend(handles=(h), frameon=False, loc='lower right') #,fontsize=20)
# ax = plt.gca().add_artist(first_legend1)

# plt.title(r"$log(E_\nu/GeV)=$"+str(energies[0]))
# plt.ylabel("CDF", fontsize=18)
# plt.xlabel(r"$z=E_{\tau}/E_{\nu}$", fontsize=18)
#plt.xlim(1e-7,1e0)

# # plt.savefig('CDF_compare_Enu'+str(energies[0])+'.jpeg', bbox_inches='tight', dpi=400)
#plt.show()


# plt.legend()
# plt.show()

#===============================================#
#=========== to plot energy ====================#
#===============================================#
# energies = [1e10] # [1e6,1e7,1e8,1e9,1e10,1e11,1e12]
# angles = np.array([1,10,20]).astype('float64')
# # e_out_type = 'cdf' # for plot type CDF
# e_out_type = 'no_cdf' # for plot type histogram

# E_nu = data.E_nu # neutrino energy bins
# color_map = cycle(sns.color_palette("colorblind")[0:len(angles)])
# lss = [':', '-.', '--', '-']

# for energy in energies:
#     fig, axs = plt.subplots(1, figsize=(10/1.1,8/1.1))
#     energy_log = float(np.log10(energy))

#     for angle in range(len(angles)):
#         e_out = data.get_clep_out(nu_type, ch_lepton, energy, angles[angle], idepth, cross_section_model, pn_model, prop_type, stats3)

#         print(len(e_out))
#         bins = E_nu # binned according to E_nu
#         count, bins_count = np.histogram(e_out, bins)
#         pdf = count / sum(count)

#         dist = count / max(count) # normalized
#         cdf = np.cumsum(dist)/np.cumsum(dist)[-1] # normalized
#         dist = np.insert(dist,0,0)
#         cdf = np.insert(cdf,0,0)

#         if e_out_type == 'no_cdf':
#             axs.step(E_nu, dist, where='pre', ls = lss[angle], color = next(color_map), label = r'$%.f^\circ$' % angles[angle], linewidth=2)
#             axs.set_xscale('log')
#             if ch_lepton=='muon':axs.set_ylabel(r"Scaled $N_{\mu}$", fontsize = 30)
#             if ch_lepton=='tau':axs.set_ylabel(r"Scaled $N_{\tau}$", fontsize = 30)

#         else:
#             axs.semilogx(E_nu, cdf, ls = '-', color = next(color_map), label = r'$%.f^\circ$' % angles[angle]) # normalized CDF
#             axs.set_ylabel("CDF", fontsize = 30)

#         if ch_lepton=='tau':
#             axs.set_xlabel(r"$E_{\tau}$ [GeV]", fontsize = 30)
#             axs.set_xlim(1e5,energy*10)
#         if ch_lepton=='muon':
#             axs.set_xlabel(r"$E_{\mu}$ [GeV]", fontsize = 30)
#             axs.set_xlim(100,energy*10)

#         axs.set_title(r"log$_{10}(E_{\nu}$/GeV) = %.2f" % energy_log +" for %.0f km depth" % idepth)

#         axs.legend(loc='upper left', ncol=int(np.ceil(len(angles)/8)), title = r'$\beta_{tr}$', frameon=False, framealpha=0.5, fontsize=27)
#         #axs.grid(which='major', axis='both', linestyle='--')
#         #axs.grid(which='minor', linestyle=':', linewidth='0.2', color='black')
#         axs.xaxis.set_ticks_position('both')
#         axs.yaxis.set_ticks_position('both')
#         axs.tick_params(axis='x', which='both', labelbottom = True, labeltop = False)
#         axs.tick_params(axis='y', which='both', left = True, labelleft = True, labelright= False)
#         axs.set_ylim(0,1.2)
#     plt.savefig('energy-dist-8-0km.pdf', bbox_inches='tight', dpi=400)
#     plt.tight_layout()
#     plt.show()

#===============================================#
#=========== to plot energy in 2x2 =============#
#===============================================#
# energies = [1e10] # [1e6,1e7,1e8,1e9,1e10,1e11,1e12]
# angles = np.array([1,10,20,35]).astype('float64')
# # e_out_type = 'cdf' # for plot type CDF
# e_out_type = 'no_cdf' # for plot type histogram

# E_nu = data.E_nu # neutrino energy bins
# color_map = sns.color_palette("colorblind")[:2]
# lss = ['-', '--']

# fig, axs = plt.subplots(2, 2, sharex=True, sharey=True)
# plt.tight_layout(pad=0.6, w_pad=0.01, h_pad=0.01)
# #plt.subplots_adjust(wspace=0.1, hspace=0.1)
# #fig.suptitle(r"log$_{10}(E_{\nu}/GeV) = %.2f$" % np.log10(energies[0]), fontsize=15)

# #=====
# eout_allm = data.get_clep_out(nu_type, ch_lepton, energies[0], angles[0], idepth, cross_section_model, pn_model, prop_type, stats)
# eout_bdhm = data.get_clep_out(nu_type, ch_lepton, energies[0], angles[0], idepth, cross_section_model, pn_models, prop_type, stats)

# bins = E_nu # binned according to E_nu
# count, bins_count = np.histogram(eout_allm, bins)
# pdf = count / sum(count)
# dist = count / max(count) # normalized
# dist = np.insert(dist,0,0)
# axs[0,0].step(E_nu, dist, where='pre', ls = lss[0], color = color_map[0], label = r'$\texttt{allm}$', linewidth=2)

# count, bins_count = np.histogram(eout_bdhm, bins)
# pdf = count / sum(count)
# dist = count / max(count) # normalized
# dist = np.insert(dist,0,0)
# axs[0,0].step(E_nu, dist, where='pre', ls = lss[1], color = color_map[1], label = r'$\texttt{bdhm}$', linewidth=2)

# axs[0,0].set_ylabel(r"Scaled $N_{\tau}$", fontsize = 22)
# axs[0,0].set_xscale('log')
# axs[0,0].legend(loc='best', title = r'$\beta_{tr} = %.f$' % (angles[0])+r'$^\circ$', frameon=False, framealpha=0.5, fontsize=18)
# axs[0,0].text(137879, 0.099, r"log$_{10}(E_{\nu}$/GeV)$=%.2f$" % np.log10(energies[0]), fontsize=20)

# #=====
# eout_allm = data.get_clep_out(nu_type, ch_lepton, energies[0], angles[1], idepth, cross_section_model, pn_model, prop_type, stats)
# eout_bdhm = data.get_clep_out(nu_type, ch_lepton, energies[0], angles[1], idepth, cross_section_model, pn_models, prop_type, stats)

# bins = E_nu # binned according to E_nu
# count, bins_count = np.histogram(eout_allm, bins)
# pdf = count / sum(count)
# dist = count / max(count) # normalized
# dist = np.insert(dist,0,0)
# axs[0,1].step(E_nu, dist, where='pre', ls = lss[0], color = color_map[0], label = r'$\texttt{allm}$', linewidth=2)

# count, bins_count = np.histogram(eout_bdhm, bins)
# pdf = count / sum(count)
# dist = count / max(count) # normalized
# dist = np.insert(dist,0,0)
# axs[0,1].step(E_nu, dist, where='pre', ls = lss[1], color = color_map[1], label = r'$\texttt{bdhm}$', linewidth=2)

# axs[0,1].set_xscale('log')
# axs[0,1].legend(loc='best', title = r'$\beta_{tr} = %.f$' % (angles[1])+r'$^\circ$', frameon=False, framealpha=0.5, fontsize=18)

# #=====
# eout_allm = data.get_clep_out(nu_type, ch_lepton, energies[0], angles[2], idepth, cross_section_model, pn_model, prop_type, stats)
# eout_bdhm = data.get_clep_out(nu_type, ch_lepton, energies[0], angles[2], idepth, cross_section_model, pn_models, prop_type, stats)

# bins = E_nu # binned according to E_nu
# count, bins_count = np.histogram(eout_allm, bins)
# pdf = count / sum(count)
# dist = count / max(count) # normalized
# dist = np.insert(dist,0,0)
# axs[1,0].step(E_nu, dist, where='pre', ls = lss[0], color = color_map[0], label = r'$\texttt{allm}$', linewidth=2)

# count, bins_count = np.histogram(eout_bdhm, bins)
# pdf = count / sum(count)
# dist = count / max(count) # normalized
# dist = np.insert(dist,0,0)
# axs[1,0].step(E_nu, dist, where='pre', ls = lss[1], color = color_map[1], label = r'$\texttt{bdhm}$', linewidth=2)

# axs[1,0].set_xlabel(r"$E_{\tau}$ [GeV]", fontsize = 22)
# axs[1,0].set_ylabel(r"Scaled $N_{\tau}$", fontsize = 22)
# axs[1,0].set_xscale('log')
# axs[1,0].legend(loc='best', title = r'$\beta_{tr} = %.f$' % (angles[2])+r'$^\circ$', frameon=False, framealpha=0.5, fontsize=18)

# #=====
# eout_allm = data.get_clep_out(nu_type, ch_lepton, energies[0], angles[3], idepth, cross_section_model, pn_model, prop_type, stats)
# eout_bdhm = data.get_clep_out(nu_type, ch_lepton, energies[0], angles[3], idepth, cross_section_model, pn_models, prop_type, stats)

# bins = E_nu # binned according to E_nu
# count, bins_count = np.histogram(eout_allm, bins)
# pdf = count / sum(count)
# dist = count / max(count) # normalized
# dist = np.insert(dist,0,0)
# axs[1,1].step(E_nu, dist, where='pre', ls = lss[0], color = color_map[0], label = r'$\texttt{allm}$', linewidth=2)

# count, bins_count = np.histogram(eout_bdhm, bins)
# pdf = count / sum(count)
# dist = count / max(count) # normalized
# dist = np.insert(dist,0,0)
# axs[1,1].step(E_nu, dist, where='pre', ls = lss[1], color = color_map[1], label = r'$\texttt{bdhm}$', linewidth=2)

# axs[1,1].set_xscale('log')
# axs[1,1].legend(loc='best', title = r'$\beta_{tr} = %.f$' % (angles[3])+r'$^\circ$', frameon=False, framealpha=0.5, fontsize=18)

# axs[1,1].set_xlabel(r"$E_{\tau}$ [GeV]", fontsize = 22)
# axs[1,1].set_ylim(0,1.05)
# axs[1,1].set_xlim(1e5,energies[0]*10)
# plt.show()

#===============================================#
#=========== to check pexit ====================#
#===============================================#
# energies=[1e7, 1e8, 1e9, 1e10, 1e11]
# #angles = np.arange(1.5,3,0.1)
# colors = sns.color_palette("colorblind")[:len(energies)]
# lss = ['-', '--']

# fig = plt.figure(figsize=[8,6], constrained_layout=True)
# spec=gridspec.GridSpec(ncols=1,nrows=3,figure=fig)  # set up subplot grid
# plt.tight_layout()

# plt.rcParams['legend.title_fontsize'] = 'small'
# plt.rcParams['xtick.labelsize']='22'
# plt.rcParams['ytick.labelsize']='22'

# plt1 = fig.add_subplot(spec[0:2,0])

# for e,c in zip(energies, colors):
#     p_no_regen = data.get_pexit(nu_type, ch_lepton, e, idepth, cross_section_model, pn_model, prop_type, stats)[1]
#     p_regen = data.get_pexit(nu_type, ch_lepton, e, idepth, cross_section_model, pn_model, prop_type, stats)[2]

#     angles = data.get_pexit(nu_type, ch_lepton, e, idepth, cross_section_model, pn_models, prop_type, stats)[0]
#     p_no_regen1 = data.get_pexit(nu_type, ch_lepton, e, idepth, cross_section_model, pn_models, prop_type, stats)[1]
#     p_regen1 = data.get_pexit(nu_type, ch_lepton, e, idepth, cross_section_model, pn_models, prop_type, stats)[2]

#     #plt.semilogy(angles,p_no_regen,color=c,label='Old no_regen')
#     plt1.loglog(angles, p_regen[9:], color=c, ls='-', marker='.', label=np.log10(e)) #'Water depth = 4 km')

#     #plt.semilogy(angles,p_no_regen1,color='green',ls='--',label='New No_regen')
#     plt1.loglog(angles, p_regen1, color=c, ls='--', marker='.', label=np.log10(e)) #'Water depth = 3 km')

# m = plt1.plot( [], [], label= r'$\log_{10}(E_{\nu}$/GeV)', ls=lss[0])
# j = [plt1.plot( [],[], color=colors[i], label= str(int(np.log10(energies[i]))), ls=lss[0])[0] for i in range(len(energies))]
# first_legend = plt1.legend(handles=j, frameon=False, loc='upper right', ncol=3, fontsize=19) #(25.08, 0.015))
# ax = plt.gca().add_artist(first_legend)

# titles = [r'$\texttt{allm}$', r'$\texttt{bdhm}$'] #['PREM-4', 'PREM-3']
# linestyles = ['-', '--']

# h = [plt1.plot( [],[], color='black', label=titles[i], ls=linestyles[i])[0] for i in range(len(titles)) ]
# plt1.legend(handles=(h), frameon=False, ncol=2, loc='lower left', fontsize=22)

# plt1.set_xlim(angles[0], angles[-1])
# plt1.set_ylabel(r'$P^{(\tau)}_{exit}$', fontsize=24)
# plt1.set_xticklabels([])
# plt1.grid(which='major', axis='both', linestyle=':', linewidth=0.5)
# plt1.grid(which='minor', axis='both', linestyle=':', linewidth=0.5)

# plt2 = fig.add_subplot(spec[2:3,0])

# for e,c in zip(energies, colors):
#     #angles1 = data.get_pexit(nu_type, ch_lepton, e, idepth, cross_section_model, pn_model, prop_type, stats3)[0]
#     p_regen = data.get_pexit(nu_type, ch_lepton, e, idepth, cross_section_model, pn_model, prop_type, stats)[2]
#     p_regen1 = data.get_pexit(nu_type, ch_lepton, e, idepth, cross_section_model, pn_models, prop_type, stats)[2]

#     plt2.semilogx(angles, p_regen[9:]/p_regen1, color=c, ls=lss[1])

# # j = [plt2.plot( [],[], color=colors[i], label= str(np.log10(energies[i])), ls=lss[1])[0] for i in range(len(energies))]
# # first_legend = plt2.legend(handles=j, title= r'$\log_{10}(E_\nu/{\rm GeV})$', ncol=2, frameon=False, loc='upper right', fontsize=18) #(25.08, 0.015))
# # ax = plt.gca().add_artist(first_legend)

# xticks = [1, 2, 3, 4, 5, 10, 15, 20, 30, 40] #5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0]
# plt2.set_ylabel(r'Ratio ($\texttt{allm}/\texttt{bdhm}$)', fontsize=24)
# plt2.set_xlabel(r'$\beta_{tr}$ [degrees]', fontsize=24)
# plt2.set_xlim(angles[0], angles[-1])
# plt2.set_xticks(xticks)
# plt2.set_xticklabels(xticks)
# #plt2.set_yticklabels(fontsize=18)
# plt2.grid(which='major', axis='both', linestyle=':', linewidth=0.5)
# plt2.grid(which='minor', axis='both', linestyle=':', linewidth=0.5)
# plt.savefig('pexit_allm_bdhm.pdf', bbox_inches='tight', dpi=400)
# plt.show()

#===============================================#
# To look at the interactions of tau, y*Ei>1PeV #
#===============================================#
# data = np.genfromtxt("lepint_9.00_ 1.0.dat")
# data = np.asarray(sorted(data, key = lambda x: x[0]))
# print(data)

# stats = [1,3,4] #np.arange(1,11,1)
# print(stats)

# colors = ['red', 'green', 'blue', 'orange', 'black', 'grey', 'cyan', 'pink', 'brown', 'purple']
# ls = ['dashed', '-', '-.']
# mk = ['o', 'x', '.']
# for j in range(len(stats)):
#     dist,energy = [], []
#     for i in range(len(data)):
#         if data[i][0] == stats[j]:
#             dist.append(data[i][1])
#             energy.append(data[i][2])

#     plt.plot(dist, energy, colors[j], marker=mk[j], linestyle=ls[j], linewidth=0.5)
#     dist = np.insert(dist, 0, 0)
#     delx=[abs(dist[k]-dist[k+1]) for k in range(len(dist)-1)]

#     print(delx)

#     #plt.plot(delx, energy, marker='o', linestyle='None') #, colors[j], marker=mk[j], linestyle=ls[j], linewidth=0.5)


# plt.xlabel("Dist travelled before interacting (km)")
# plt.ylabel("Shower energy, y*Ei (GeV)")
# plt.title(r"$E_{i,\tau}=10^9$ GeV")
# plt.show()

#===================================================#
# To plot col depth for different water layer depth #
#===================================================#
# cols = data.get_trajs('col', 1, 1, out=False)
# cols1 = data.get_trajs('col', 1, 0, out=False)
# data_1m = data_extract("1m_data.dat")
# #data_0km = data_extract("0km_data.dat")

# print((cols[1][0]))

# xalong_1m, cdalong_1m = [], []
# for i in range(len(data_1m)):
#     if data_1m[0][i]==1:
#         xalong_1m.append(data_1m[1][i])
#         cdalong_1m.append(data_1m[2][i])

# print(cdalong_1m)

# for j in range(len(cols)):
#     if abs(cols[1][j] - data_1m[2][j]) < 0.01:
#         print("true")

# print(len(cdalong_1m))

# xalong_0km, cdalong_0km = [], []
# for i in range(len(data_0km)):
#     if data_0km[0][i]==10:
#         xalong_0km.append(data_0km[1][i])
#         cdalong_0km.append(data_0km[2][i])

# angles = np.arange(0.1,90.1,0.1)
# cols=[]
# for a in angles:
#     print("%.1f" % a)
#     chord, water = (data.get_trajs('water', "%.1f"%a, 1, out=False))
    #cols.append((chord, water))
#plt.plot(cols[0], cols[1], label='1km nupy')
# plt.plot(cols1[0], cols1[1], label='0km nupy')
# plt.plot(cols1[0], cols1[1], ls='dotted', label='1km')
# plt.plot(xalong_1m, cdalong_1m, ls='--', label='1m')
# plt.plot(xalong_0km, cdalong_0km, ls='-', label='0km')

# plt.legend()
# plt.show()

#print(Geometry.sagitta_deg(1))

# beta, x, cd = Geometry.gen_col_trajs(4)
# # beta, x, cd = Geometry.gen_water_trajs(0.8)

# datas = np.stack([beta, x, cd], axis=1)
# np.savetxt("col_4.dat", datas, delimiter="\t")

# #plt.plot(x, cd, marker='.', label='homemade')
# plt.legend()
# plt.show()

# angle = 1;
# hallsie = data_extract("traj-columndepth0p1.dat")
# diksha = data_extract("col_0p5.dat")
# diksha0 = data_extract("col_0.dat")
# diksha0p5 = data_extract("col_4.dat")

# x_dik, col_dik = [], []
# x_dik0, col_dik0 = [], []
# x_dik0p5, col_dik0p5 = [], []
# for a in range(len(diksha)):
#     if diksha[0][a]==angle:
#         #print(diksha[2][a])
#         x_dik.append(diksha[1][a])
#         col_dik.append(diksha[2][a])
#         x_dik0.append(diksha0[1][a])
#         col_dik0.append(diksha0[2][a])
#         x_dik0p5.append(diksha0p5[1][a])
#         col_dik0p5.append(diksha0p5[2][a])
# plt.semilogy(x_dik, col_dik, label="Diksha 0.1 km")
# plt.semilogy(x_dik0, col_dik0, ls='--', label="Diksha 0 km")
# plt.semilogy(x_dik0p5, col_dik0p5, marker='.', ls='none', label="Diksha 4 km")

# x_hal, col_hal = [], []
# for a1 in range(len(hallsie)):
#     if hallsie[0][a1]==angle:
#         x_hal.append(hallsie[1][a1])
#         col_hal.append(hallsie[2][a1])
# plt.plot(x_hal, col_hal, marker='.', ls='none', label="Hallsie 0.5")

# plt.legend()
# plt.show()

#a = Geometry.create_traj_table(0)
#data_1m = np.stack([a[0], a[1], a[2]], axis=1)
#np.savetxt("0km_data.dat", data_1m, delimiter="\t")
#print(a[0])

#print( Geometry.find_interface(0.8) )