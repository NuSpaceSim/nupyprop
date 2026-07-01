#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 14:04:46 2021

@author: Diksha Garg and Sameer Patel
Comments: This is the main code where calls happen for propagating 
neutrinos, charged leptons for each incoming neutrino energy and angle. 
It also creates the final output file. 
"""

import nupyprop.data as Data
import nupyprop.geometry as Geometry
import nupyprop.run as Run
import nupyprop.constants as const

import numpy as np
from collections.abc import Iterable
import time
import os
import glob

rho_water = const.rho_water # g/cm^3

def init_xc(nu_type, ch_lepton, nu_model, pn_model, prop_type):
    '''

    Parameters
    ----------
    nu_type : str
        Type of neutrino particle. Can be neutrino or anti_neutrino.
    ch_lepton : str
        Type of charged lepton. Can be tau or muon.
    nu_model : str
        Neutrino cross-section model.
    pn_model : str
        Photonuclear energy loss model.
    prop_type : str
        Energy loss propagation type. Can be stochastic or continuous.

    Returns
    -------
    nu_xc : ndarray
        2D array containing CC & NC neutrino cross-section values, in cm^2.
    xc_water : ndarray
        2D array containing bremmstrahlung, pair production & photonuclear energy loss cross-section values for water, in cm^2.
    xc_rock : ndarray
        2D array containing bremmstrahlung, pair production & photonuclear energy loss cross-section values for rock, in cm^2.
    alpha_water : ndarray
        1D array containing ionization energy loss values for water, in (GeV*cm^2)/g.
    alpha_rock : ndarray
        1D array containing ionization energy loss values for rock, in (GeV*cm^2)/g.
    beta_water : ndarray
        2D array containing bremmstrahlung, pair production & photonuclear energy loss values for water, in cm^2/g.
    beta_rock : ndarray
        2D array containing bremmstrahlung, pair production & photonuclear energy loss values for rock, in cm^2/g.

    '''
    if prop_type == 'stochastic':beta_type = 'cut'
    elif prop_type == 'continuous':beta_type = 'total'

    nu_xc = Data.get_xc('nu', nu_model, nu_type)

    xc_water = Data.get_xc(ch_lepton, pn_model, 'water')
    alpha_water = Data.get_alpha(ch_lepton, 'water')
    beta_water = Data.get_beta(ch_lepton, pn_model, 'water', beta_type)

    xc_rock = Data.get_xc(ch_lepton, pn_model, 'rock')
    alpha_rock = Data.get_alpha(ch_lepton, 'rock')
    beta_rock = Data.get_beta(ch_lepton, pn_model, 'rock', beta_type)

    return nu_xc,xc_water,xc_rock,alpha_water,alpha_rock,beta_water,beta_rock

def init_ixc(nu_type, ch_lepton, nu_model, pn_model):
    '''

    Parameters
    ----------
    nu_type : str
        Type of neutrino particle. Can be neutrino or anti-neutrino.
    ch_lepton : str
        Type of charged lepton. Can be tau or muon.
    nu_model : str
        Neutrino cross-section model.
    pn_model : str
        Photonuclear energy loss model.

    Returns
    -------
    nu_ixc : ndarray
        3D array containing neutrino integrated cross-section CDF values.
    lep_ixc_water : ndarray
        3D array containing lepton integrated cross-section CDF values for water.
    lep_ixc_rock : ndarray
        3D array containing lepton integrated cross-section CDF values for rock.

    '''
    nu_ixc = Data.get_ixc('nu', nu_model, nu_type)

    lep_ixc_water = Data.get_ixc(ch_lepton, pn_model, 'water')

    lep_ixc_rock = Data.get_ixc(ch_lepton, pn_model, 'rock')

    return nu_ixc, lep_ixc_water, lep_ixc_rock


def main(E_prop, angles, nu_type, cross_section_model, pn_model, earth_model, idepth, ch_lepton, fac_nu, stats, prop_type, elep_mode, htc_mode, job_num = 0):
    '''
    Parameters
    ----------
    E_prop : ndarray
        Ingoing neutrino energies, in GeV.
    angles : ndarray
        Earth emergence angles (beta), in degrees.
    nu_type : str
        Type of neutrino particle. Can be neutrino or anti_neutrino.
    cross_section_model : str
        Neutrino cross-section model.
    pn_model : str
        Photonuclear energy loss model.
    earth_model : str
        Earth density model.
    idepth : int
        Depth of water layer in km.
    ch_lepton : str
        Type of charged lepton. Can be tau or muon.
    fac_nu : float
        Rescaling factor for BSM cross-sections.
    stats : float
        Statistics; no. of ingoing neutrinos.
    prop_type : str
        Energy loss propagation type. Can be stochastic or continuous.
    elep_mode : str
        Option to print exiting charged lepton's final energy. Can be yes or no.
    htc_mode : str
        Opion to create files such that they can be used for 
        high throughput computing (htc) mode, Can be yes or no.
    job_num : int
        Gives the job number running on cluster (only when htc mode is on. Default is 0.

    Returns
    -------
    NONE
        Adds exit probability, tau polarization & outgoing lepton energy results to output_x.h5

    '''
    make_array = lambda x : x if isinstance(x, Iterable) else np.array([x]) # to avoid errors with single or no columns

    nu_xc, xc_water, xc_rock, alpha_water, alpha_rock, beta_water, beta_rock = init_xc(nu_type, ch_lepton, cross_section_model, pn_model, prop_type)

    nu_ixc, lep_ixc_water, lep_ixc_rock = init_ixc(nu_type, ch_lepton, cross_section_model, pn_model)

    ithird = const.ithird # use dn/dy in tau to neutrino

    if prop_type == 'stochastic':
        prop_type_int = 1
    else:
        prop_type_int = 2

    if ch_lepton == 'tau':
        lepton_int = 1
    else: # muon
        lepton_int = 2

    out_file = Data.output_file(nu_type,ch_lepton,idepth,cross_section_model,pn_model,earth_model,prop_type,stats,arg=None)

    print("The water -> rock transition occurs at %.2f degrees" % Geometry.find_interface(idepth)[0])
    
    for energy in sorted(E_prop):
        eout_list = []
        for angle in sorted(angles):
            start_time = time.time() # start time of simulations
            print("Neutrino Energy = 10^(%.2f) GeV, Earth Emergence Angle = %.2f degrees" % (energy, angle), flush=True)

            xalong, cdalong = Data.get_trajs('col', angle, idepth, earth_model)
            
            _, water = Data.get_trajs('water', angle, idepth, earth_model)
            dwater = water*rho_water # depth in water [kmwe] in last or only section
            depthE = Geometry.columndepth(angle, idepth, model_name=earth_model)*1e-5 # column depth in kmwe
            
            results_nr, results_r, e_out, p_out = Run.run_stat(
                10**energy, angle, nu_xc, nu_ixc, depthE, dwater, xc_water, xc_rock,
                lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock,
                xalong, cdalong, ithird, idepth, lepton_int, fac_nu, prop_type_int, stats, earth_model
            )
            
            # polarization
            P_out_exit = make_array(p_out)
            P_avg = -1 if P_out_exit.size == 0 else Data.sign(np.mean(P_out_exit))

            if htc_mode == 'no': #htc mode off
                prob_no_regen = float(np.sum(results_nr)/stats) # exit prob. with no regeneration
                prob_regen = float(np.sum(results_r)/stats) # exit prob. with regeneration

                Data.add_pexit(ch_lepton, energy,
                               np.array([angle]),
                               np.array([prob_no_regen]),
                               np.array([prob_regen]),
                               out_file)

                Data.add_polarization(ch_lepton, energy,
                                      np.array([angle]),
                                      np.array([P_avg]),
                                      out_file)
                    
                #exiting lepton's final energy
                eout_vals = np.insert(e_out[e_out > 0.0], 0, angle)
                eout_list.append(eout_vals)
                
                if elep_mode == 'yes':
                    print("I am going to print final energies in output file.")

                    e_out_exit = Data.patch_for_astropy(make_array(e_out[e_out > 0.0]))
                    Data.add_clep_out(ch_lepton, energy, angle, e_out_exit, out_file)

                Data.add_cdf(ch_lepton, energy, eout_list, out_file, htc_mode=False, arg=None)
            
            else: # save p_exit files with single energy and angle; 
                  # no post-processing of files here in the main execution (for that, use data.process_htc_out)
                no_regen = np.sum(results_nr)
                regen = np.sum(results_r)
            
                with open("pexit_{:.2f}_{:.1f}_{:d}.dat".format(energy, angle, job_num), "w") as pexit_file:
                    pexit_file.write("%.5e\t%.5e\t%.5e\t%.5e\n" % (10**energy, angle, no_regen, regen))

                eout_vals = e_out[e_out > 0.0]
                with open("eout_{:.2f}_{:.1f}_{:d}.dat".format(energy, angle, job_num), "w") as eout_file:
                    for e in eout_vals:
                        eout_file.write("%.5e\n" % e)

                with open("Pout_{:.2f}_{:.1f}_{:d}.dat".format(energy, angle, job_num), "w") as pout_file:
                    pout_file.write("%.5e\n" % P_avg)

            end_time = time.time() # end time for simulations
            print(f"{end_time-start_time:.2f} seconds to compute \n"
                  f" for Neutrino Energy = 10^({energy:.2f}) GeV,"
                  f" Earth Emergence Angle = {angle:.2f} degrees")
# =============================================================================
#
# =============================================================================
if __name__ == "__main__":
    pass

    # for i in range(10):
    # random.seed(30)
    start_time = time.time()
    e_nu = np.array([7])
    angles = np.array([1])

    idepth = 4
    fac_nu = 1
    ch_lepton = 'tau'
    cross_section_model = 'ctw'
    pn_model = 'allm'
    prop_type = 'stochastic'
    stats = int(1e6)
    nu_type = 'neutrino'
    htc_mode = 'no'
    earth_model = 'prem'
    elep_mode = 'no'

    main(e_nu, angles, nu_type, cross_section_model, pn_model,earth_model, idepth, ch_lepton, fac_nu, stats, 
         prop_type, elep_mode, htc_mode)

    end_time = time.time()
    print(f"It took a total of {end_time-start_time:.2f} seconds to compute")
