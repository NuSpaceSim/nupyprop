#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 14:04:46 2021

@author: sam
"""

import data as Data
import geometry as geom_py

from propagate import interpolation as Interpolation
from propagate import geometry as Geometry
from propagate import transport as Transport
from propagate import run as Run

import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import time
import os
import glob
import functools
print = functools.partial(print, flush=True)

rho_water = 1.02 # g/cm^3
rho_rock = 2.65 # g/cm^3

def init_xc(nu_type, lepton, nu_model, pn_model, prop_type):
    pn_model = 'pn_' + pn_model

    nu_xc = Data.get_xc('nu', nu_model, particle=nu_type)

    xc_water = Data.combine_lep('xc', lepton, 'water', pn_model)

    xc_rock = Data.combine_lep('xc', lepton, 'rock', pn_model)

    alpha_water = Data.get_alpha(lepton, 'water')

    alpha_rock = Data.get_alpha(lepton, 'rock')

    if prop_type == 'stochastic':
        beta_water = Data.combine_lep('beta', lepton, 'water', pn_model, beta_type='cut')
        beta_rock = Data.combine_lep('beta', lepton, 'rock', pn_model, beta_type='cut')
    elif prop_type == 'continuous':
        beta_water = Data.combine_lep('beta', lepton, 'water', pn_model, beta_type='total')
        beta_rock = Data.combine_lep('beta', lepton, 'rock', pn_model, beta_type='total')

    return nu_xc,xc_water,xc_rock,alpha_water,alpha_rock,beta_water,beta_rock

def init_ixc(nu_type, lepton, nu_model, pn_model):
    pn_model = 'pn_' + pn_model

    ixc_nu = Data.get_ixc('nu', nu_model, particle=nu_type)

    ixc_water = Data.combine_lep('ixc', lepton, 'water', pn_model)

    ixc_rock = Data.combine_lep('ixc', lepton, 'rock', pn_model)

    return ixc_nu, ixc_water, ixc_rock

def create_lep_out_dict(energy, angle):

    fnm = str("e_out_%.2E_%.2f" % (energy, angle))

    e_out = np.genfromtxt(fnm)

    lep_dict = {'lep_energy': e_out}

    return lep_dict

def main(E_prop, angles, nu_type, cross_section_model, pn_model, idepth, lepton, fac_nu, stat, prop_type):



    nu_xc, xc_water, xc_rock, alpha_water, alpha_rock, beta_water, beta_rock = init_xc(nu_type, lepton, cross_section_model, pn_model, prop_type)


    nu_ixc, lep_ixc_water, lep_ixc_rock = init_ixc(nu_type, lepton, cross_section_model, pn_model)


    ithird = 0 # use dn/dy in tau to neutrino

    prob_dict = {}
    angle_arr = []
    no_regen_arr = []
    regen_arr = []

    if prop_type == 'stochastic':
        prop_type_int = 1
    else:
        prop_type_int = 2

    if lepton == 'tau':
        lepton_int = 1
    else:
        lepton_int = 2

    if nu_type == 'neutrino':
        nu_type = 'nu'
    elif nu_type == 'antii-neutrino':
        nu_type = 'anu'

    chk_flag = Data.chk_file(nu_type, lepton, idepth, cross_section_model, pn_model, prop_type, stat)

    if chk_flag == 0:
        return print("No changes made")

    print("The water -> rock transition occurs at %.2f degrees" % geom_py.find_interface(idepth)[0])

    for energy in E_prop:
        start_time = time.time()
        for angle in angles:
            xalong, cdalong = Data.get_trajs('col', angle, idepth) # initialize arrays here for each angle, to reduce a ton of overhead when tauthrulayers & regen are called

            print("Energy = %.2e, Angle = %.1f" % (energy, angle))

            chord, water = Data.get_trajs('water', angle, idepth)
            dwater = water*rho_water # depth in water [kmwe] in last or only section
            depthE = geom_py.columndepth(angle, idepth)*1e-5 # column depth in kmwe?


            no_regen, regen = Run.run_stat_single(energy, angle, nu_xc, nu_ixc, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton_int, fac_nu, stat, prop_type_int)


            prob_no_regen = no_regen/float(stat)
            prob_regen = regen/float(stat)

            angle_arr.append(angle)
            no_regen_arr.append(prob_no_regen)
            regen_arr.append(prob_regen)

            lep_dict = create_lep_out_dict(energy, angle)
            Data.add_lep_out(nu_type, lepton, energy, angle, lep_dict, idepth, cross_section_model, pn_model, prop_type, stat)

        # # end of for loop for angles

        prob_dict_single = {'angle':np.asarray(angle_arr),'no_regen':np.asarray(no_regen_arr),'regen':np.asarray(regen_arr)}
        Data.add_pexit(nu_type, lepton, energy, prob_dict_single, idepth, cross_section_model, pn_model, prop_type, stat)

        angle_arr = []
        no_regen_arr = []
        regen_arr = []

        prob_dict[energy]={'no_regen':prob_dict_single['no_regen'],'regen':prob_dict_single['regen']}


        print("Exit Probability and lepton energy out lookup tables successfully created")
        end_time = time.time()
        print(f"It took {end_time-start_time:.2f} seconds to compute")


    # close of for loop for energy

    Data.add_cdf(E_prop, angles, nu_type, lepton, idepth, cross_section_model, pn_model, prop_type, stat)

    return prob_dict, lep_dict

# =============================================================================
#
# =============================================================================
if __name__ == "__main__":
    # for i in range(10):
    # random.seed(30)
    start_time = time.time()
    angles = np.array([1])
    # angles = np.arange(1,36)
    # angles = np.array([1,2,3,4,5])
    # angles = np.array([1,2,3,5,7,10,12,15,17,20,25,30,35])
    # angles = np.array([17,20,25,30,35])
    E_prop = np.array([10**7])

    idepth = 4
    fac_nu = 1
    lepton = 'tau'
    cross_section_model = 'ct18nlo'
    pn_model = 'bb'
    prop_type = 'stochastic'
    stat = int(1e8)
    nu_type = 'neutrino'

    # nu_xc, xc_water, xc_rock, alpha_water, alpha_rock, beta_water, beta_rock = init_xc(lepton, cross_section_model, pn_model, prop_type)

    # nu_ixc, lep_ixc_water, lep_ixc_rock = init_ixc(lepton, cross_section_model, pn_model)

    # prob_dict, lep_dict = main(E_prop, angles, nu_type, cross_section_model, pn_model, idepth, lepton, fac_nu, stat, prop_type)
    main(E_prop, angles, nu_type, cross_section_model, pn_model, idepth, lepton, fac_nu, stat, prop_type)

    # print(prob_dict)

    e_out_files = glob.glob("e_out_*")

    for file in e_out_files:
        os.remove(file)

    # with ProcessPoolExecutor(max_workers = 4) as executor:
    #     prob_dict, e_out = executor.map(main, angles, cross_section_model, pn_model, [idepth], lepton, [fac_nu], [stat], prop_type)

    end_time = time.time()
    print(f"It took a total of {end_time-start_time:.2f} seconds to compute")
    # print("str(%.2E)" % 1e7.5)