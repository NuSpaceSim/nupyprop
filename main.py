#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 9 23:04:35 2020

@author: sam
"""
import os

os.environ['KMP_DUPLICATE_LIB_OK']='True'

import data as Data
import geometry as Geometry
import energy_loss as Energy_loss
import my_interpolation as Interpolation
import transport as Transport
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

import deepdiff
import copy
import numpy as np
import time
import numba as nb
# from numba.typed import Dict
# from numba import njit, prange
import functools
print = functools.partial(print, flush=True)

E_nu = Data.E_nu
E_lep = Data.E_lep
rho_water = 1.02 # g/cm^2
rho_rock = 2.65 # g/cm^2

def init_xc(lepton, nu_model, pn_model, prop_type='stochastic'):
    nu_xc = Data.get_xc('nu', nu_model, particle='neutrino')
    xc_water = Data.get_xc(lepton, pn_model, material='water')
    xc_rock = Data.get_xc(lepton, pn_model, material='rock')
    alpha_water = Data.get_alpha(lepton, 'water')
    alpha_rock = Data.get_alpha(lepton, 'rock')
    if prop_type == 'stochastic':
        beta_water = Data.get_beta(lepton, 'water', 'continuous', pn_model) # change continuous to cut
        beta_rock = Data.get_beta(lepton, 'rock', 'continuous', pn_model)
    elif prop_type == 'continuous':
        beta_water = Data.get_beta(lepton, 'water', 'total', pn_model)
        beta_rock = Data.get_beta(lepton, 'rock', 'total', pn_model)
    return nu_xc,xc_water,xc_rock,alpha_water,alpha_rock,beta_water,beta_rock

def ixc_nb(ixc_dict):
    if 'cc' in ixc_dict.keys():
        models = ['cc', 'nc']
        energies = E_nu
    else:
        models = ['brem', 'pair', 'pn']
        energies = E_lep

    ind_dict = {}
    en_dict = {}
    ixc = {}
    for model in models:
        en_dict = {}
        for j in energies:
            ind_dict = {}
            for i in range(1,31):
                ind_dict[i] = ixc_dict[model][j][i]
            en_dict[j] = ind_dict
        ixc[model] = en_dict
    return ixc

def init_ixc(lepton, nu_model, pn_model):
    ixc_nu = Data.get_nu_iixc(nu_model, particle='neutrino')
    ixc_water = Data.get_lep_iixc(pn_model, 'water')
    ixc_rock = Data.get_lep_iixc(pn_model, 'rock')
    return ixc_nu, ixc_water, ixc_rock

# @njit(nogil=True)
def bin_data(angle, energy, eb_no_regen, eb_regen):

    emid = 10**np.asarray([(float(i+30)-0.5)/10.0 for i in range(1,92)])
    prob_no_regen = np.cumsum(eb_no_regen)[-1]/float(stat)
    prob_regen = np.cumsum(eb_regen)[-1]/float(stat)

    return energy, angle, prob_no_regen, prob_regen, emid

# @njit(nogil=True, parallel=True)
def run_stat(energy, angle, nu_xc, nu_ixc, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird): # depthE is the total column depth from PREM at an angle
    regen_cnt = 0
    no_regen_tot = 0
    regen_tot = 0

    cd2distd = Interpolation.f_cd2distd(xalong, cdalong)
    densityatx = Geometry.f_densityatx(angle)

    # for i in range(1,stat+1):

    depth = np.full(stat, depthE)
    # depth0 = np.zeros(stat) # start with this each time

    # 80 continue

    # tnu goes until neutrino either goes to dtot, or converts to a tau
    ip, dtr, ef = Transport.propagate_nu(energy, nu_xc, nu_ixc, depth)

    # if ip == 'nu': # still a neutrino at the end of the road
    #     # go to 10
    #     continue # break outside stat; continue is correct here
    # continue here: we have a tau
    nu_mask = ip == 'nu'
    tau_mask = ~nu_mask
    taus_count = np.count_nonzero(tau_mask)
    e_out = np.zeros(taus_count)

    depth = depth[tau_mask]
    dtr = dtr[tau_mask]
    etauin = ef[tau_mask]

    # how far did the neutrino go? dtr is how far traveled
    # depth0 += dtr # how far is the neutrino on trajectory?
    dleft = depthE - dtr # how far is left for the neutrino to travel?


    # still need to propagate the tau, column depth to go

    ipp, dfinal, etauf = Transport.tau_thru_layers(angle, depth, dwater, dtr, etauin, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, cd2distd, densityatx) # note: angle is now in betad

    dleft = depth-dfinal

    # if ipp == 'not_decayed' and dleft <= 0: # a tau has emerged through column depth
    #     no_regen_tot += 1
    #     regen_tot += 1 # update the regen tau array once
    #     e_out.append(etauf)
    #     # go to 10; we are done with the loop
    #     continue # break outside stat; continue is correct here
    full_emerged_mask = (ipp == 'not_decayed') & (dleft <= 0)
    no_regen_tot = np.count_nonzero(full_emerged_mask)
    regen_tot = no_regen_tot

    e_out[full_emerged_mask] = etauf[full_emerged_mask]

    unemerged_mask = ~full_emerged_mask
    unemerged_count = np.count_nonzero(unemerged_mask)

    # 11 continue; beginning of regeneration loop
    # must be a neutrino. Is there still column depth to propagate?
    etauf = etauf[unemerged_mask]
    depth = depth[unemerged_mask]
    dfinal = dfinal[unemerged_mask]

    regen_cnt = 1 # tau out after first interaction
    ipp3 = np.full('dummy_value', unemerged_count)
    # regen_emerged_mask = unemerged_mask
    while np.any(dfinal < depthE) and np.any(ipp3!='not_decayed') and regen_cnt<=10: # tau has decayed before the end

        etauin = etauf # regen finds neutrino energy

        ipp3, dtau2, ef2 = Transport.regen(angle, etauin, depth, dwater, dfinal, nu_xc, nu_ixc, ithird, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, cd2distd, densityatx) # note: angle is now in betad

        regen_cnt +=1

        regen_emerged_mask = (ipp3 == 'not_decayed') # size of subset, true where not decayed
        regen_tot += np.count_nonzero(regen_emerged_mask)

        e_out[unemerged_mask][regen_emerged_mask] = etauf[regen_emerged_mask]

        regen_unemerged_mask = ~regen_emerged_mask
        unemerged_mask[unemerged_mask] = regen_unemerged_mask

        unemerged_count = np.count_nonzero(regen_unemerged_mask)
        ipp3 = np.full('dummy_value', unemerged_count)

        etauf = ef2[unemerged_mask]
        depth = depth[unemerged_mask]
        dfinal = dtau2[unemerged_mask]

    e_out = e_out[~unemerged_mask]

    return no_regen_tot, regen_tot, e_out

def main():

    nu_xc, xc_water, xc_rock, alpha_water, alpha_rock, beta_water, beta_rock = init_xc(lepton, cross_section_model, pn_model, prop_type='stochastic')

    nu_ixc, lep_ixc_water, lep_ixc_rock = init_ixc(lepton, cross_section_model, pn_model)

    ithird = 0 # use dn/dy in tau to neutrino

    prob_dict = {}
    angle_arr = []
    no_regen_arr = []
    regen_arr = []

    print("The water -> rock transition occurs at %.2f degrees" % Geometry.find_interface()[0])

    for energy in E_prop:
        for angle in angles:
            xalong, cdalong = Data.get_trajs('col', angle, idepth) # initialize arrays here for each angle, to reduce a ton of overhead when tauthrulayers & regen are called

            print("Energy = %.0e, Angle = %.d" % (energy, angle))

            chord, water = Data.get_trajs('water', angle, idepth)
            dwater = water*rho_water # depth in water [kmwe] in last or only section
            depthE = Geometry.columndepth(angle)*1e-5 # column depth in g/cm^2


            no_regen, regen, e_out = run_stat(energy, angle, nu_xc, nu_ixc, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird)

            prob_no_regen = no_regen/float(stat)
            prob_regen = regen/float(stat)

            angle_arr.append(angle)
            no_regen_arr.append(prob_no_regen)
            regen_arr.append(prob_regen)

            lep_dict = {'lep_energy':np.asarray(e_out)}
            # Data.add_lep_out(energy, angle, lep_dict)

        # # end of for loop for angles

        prob_dict_single = {'angle':np.asarray(angle_arr),'no_regen':np.asarray(no_regen_arr),'regen':np.asarray(regen_arr)}
        # Data.add_pexit(energy, prob_dict_single)

        angle_arr = []
        no_regen_arr = []
        regen_arr = []

        prob_dict[energy]={'no_regen':prob_dict_single['no_regen'],'regen':prob_dict_single['regen']}

        print("Exit Probability and lepton energy out lookup tables successfully created")


    # close of for loop for energy

    return prob_dict, e_out


# =============================================================================
#
# =============================================================================
if __name__ == "__main__":
    # ray.init()
    # random.seed(30)
    start_time = time.time()
    angles = np.array([5])
    # angles = np.array([1,3,5,7,10,12,15,17,20,25,30,35])
    E_prop = np.array([1e7])

    idepth = 4
    Geometry.idepth = idepth
    fac_nu = 1
    m_le = Transport.m_le = Energy_loss.m_tau
    lepton = 'tau'
    # material = 'rock'
    cross_section_model = 'ncteq15'
    pn_model = 'allm'
    stat = int(1e6)
    Transport.fac_nu = fac_nu
    c_tau = Transport.c_tau = 8.703e-3

    nu_xc, xc_water, xc_rock, alpha_water, alpha_rock, beta_water, beta_rock = init_xc(lepton, cross_section_model, pn_model, prop_type='stochastic')

    # nu_ixc, lep_ixc_water, lep_ixc_rock = init_ixc(lepton, cross_section_model, pn_model)
    prob_dict, e_out = main()
    
    end_time = time.time()
    print(f"It took {end_time-start_time:.2f} seconds to compute")
