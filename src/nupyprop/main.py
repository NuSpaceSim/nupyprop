#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 14:04:46 2021

@author: sam
"""

import nupyprop.data as Data
import nupyprop.geometry as Geometry
from nupyprop.propagate import run as Run

import numpy as np
from astropy.table import Table
from collections import OrderedDict
from collections.abc import Iterable
import time
import os
import glob

rho_water = 1.02 # g/cm^3

def file_cleaner(output_type):
    '''

    Parameters
    ----------
    output_type : str
        Type of output file to clean. Can be e_out, P_out, polarization_* or p_exit.

    Returns
    -------
    NONE
        Cleans temporary generated output files.

    '''
    if output_type == 'e_out':
        files = glob.glob("eout_*") # cleanup of Fortran e_out files
    elif output_type == 'P_out':
        files = glob.glob("Pout_*") # cleanup of Fortran polarization files
    elif output_type == 'polarization':
        files = glob.glob("polarization_*")
    elif output_type == 'p_exit':
        files = glob.glob("pexit_*") # cleanup of p_exit files

    for file in files:
            os.remove(file)

    return None

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


def main(E_prop, angles, nu_type, cross_section_model, pn_model, idepth, ch_lepton, fac_nu, stats, prop_type, elep_mode, htc_mode):
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

    Returns
    -------
    NONE
        Adds exit probability & outgoing lepton energy results to output_x.h5

    '''
    make_array = lambda x : x if isinstance(x, Iterable) else np.array([x]) # to avoid errors with single or no columns

    nu_xc, xc_water, xc_rock, alpha_water, alpha_rock, beta_water, beta_rock = init_xc(nu_type, ch_lepton, cross_section_model, pn_model, prop_type)


    nu_ixc, lep_ixc_water, lep_ixc_rock = init_ixc(nu_type, ch_lepton, cross_section_model, pn_model)

    ithird = 0 # use dn/dy in tau to neutrino

    if prop_type == 'stochastic':
        prop_type_int = 1
    else:
        prop_type_int = 2

    if ch_lepton == 'tau':
        lepton_int = 1
    else: # muon
        lepton_int = 2

    out_file = Data.output_file(nu_type,ch_lepton,idepth,cross_section_model,pn_model,prop_type,stats,arg=None)

    start_time = time.time()

    print("The water -> rock transition occurs at %.2f degrees" % Geometry.find_interface(idepth)[0])

    for energy in sorted(E_prop):
        # log_energy = np.log10(energy)
        eout_list = [] #to store final energies of exiting charged lepton
        for angle in sorted(angles):

            xalong, cdalong = Data.get_trajs('col', angle, idepth) # initialize arrays here for each angle, to reduce a ton of overhead when tauthrulayers & regen are called

            print("Neutrino Energy = 10^(%.2f) GeV, Earth Emergence Angle = %.2f degrees" % (energy, angle), flush=True)

            chord, water = Data.get_trajs('water', angle, idepth)
            dwater = water*rho_water # depth in water [kmwe] in last or only section
            depthE = Geometry.columndepth(angle, idepth)*1e-5 # column depth in kmwe

            no_regen, regen = Run.run_stat_single(10**energy, angle, nu_xc, nu_ixc, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton_int, fac_nu, stats, prop_type_int)

            prob_no_regen = no_regen/float(stats)
            prob_regen = regen/float(stats)

            if htc_mode == 'no': # HTC mode off
                with open("pexit_%.2f.dat" % energy, "a") as pexit_file:
                    pexit_file.write("%.5e\t%.5e\t%.5e\t%.5e\n" % (10**energy,angle,prob_no_regen,prob_regen))

                P_out = make_array(np.genfromtxt(str("Pout_{:.2f}_{:4.1f}.dat".format(energy, angle))))  #polarization of the exiting charged leptons
                if P_out.size==0:
                    P_avg = -1
                else:
                    P_avg = Data.sign(np.mean(P_out))  # avg of polarization of the exiting charged leptons

                with open("polarization_%.2f.dat" % energy, "a") as pola_file:
                    pola_file.write("%.5e\t%.5e\t%.5e\n" % (10**energy,angle,P_avg))

                if elep_mode == 'yes':
                    print("Printing final lepton energies in the output file.")
                    e_out = make_array(np.genfromtxt(str("eout_{:.2f}_{:4.1f}.dat".format(energy, angle))))
                    e_out = Data.patch_for_astropy(e_out)
                    Data.add_clep_out(ch_lepton, energy, angle, e_out, out_file)

                # exiting lepton's final energy for CDFs
                eout_vals = np.genfromtxt("eout_{:.2f}_{:4.1f}.dat".format(energy, angle))
                eout_vals = np.insert(eout_vals, 0, angle)
                eout_list.append(eout_vals)

            else: # HTC mode on. Post-processing can be done using data.process_htc_out
                with open("pexit_%.2f_%.1f.dat" % (energy,angle), "w") as pexit_file:
                    pexit_file.write("%.5e\t%.5e\t%.5e\t%.5e\n" % (10**energy,angle,prob_no_regen,prob_regen))

                    eout_name = "eout_{:.2f}_{:4.1f}.dat".format(energy, angle)
                    Pout_name = "Pout_{:.2f}_{:4.1f}.dat".format(energy, angle)
                    os.rename(eout_name, "eout_{:.2f}_{:.1f}.dat".format(energy, angle))
                    os.rename(Pout_name, "Pout_{:.2f}_{:.1f}.dat".format(energy, angle))

        # end of for loop for angles
        if htc_mode == 'no': # HTC mode off; do some post-processing
            # Pexit of charged leptons
            p_angle, p_noregen, p_regen = np.genfromtxt("pexit_%.2f.dat" % energy, usecols=(1,2,3), unpack=True)

            p_angle = Data.patch_for_astropy(p_angle)
            p_noregen = Data.patch_for_astropy(p_noregen)
            p_regen = Data.patch_for_astropy(p_regen)

            Data.add_pexit(ch_lepton, energy, p_angle, p_noregen, p_regen, out_file)
            file_cleaner('p_exit')  # remove p_exit files

            # Avg polarization of the exiting charged leptons
            pola_angle, avg_pola = np.genfromtxt("polarization_%.2f.dat" % energy, usecols=(1,2), unpack=True)

            pola_angle = Data.patch_for_astropy(pola_angle)
            avg_pola = Data.patch_for_astropy(avg_pola)

            Data.add_polarization(ch_lepton, energy, pola_angle, avg_pola, out_file)
            file_cleaner('P_out') # remove P_out files
            file_cleaner('polarization') #remove polarization_* files

            # energy CDFs for the exiting charged leptons
            Data.add_cdf(ch_lepton, energy, eout_list, out_file, htc_mode=False, arg=None) # adds the binned cdf values for all neutrino energies and angles to the output file.
            file_cleaner('e_out') # remove e_out files

        end_time = time.time()
        print(f"It took {end_time-start_time:.2f} seconds to compute")

        if htc_mode== 'no': print("Done!")
        else: print("Done!") # for HTC mode on

    return None

# =============================================================================
#
# =============================================================================
if __name__ == "__main__":
    # for i in range(10):
    # random.seed(30)
    start_time = time.time()
    e_nu = np.array([7])
    angles = np.array([1,2,10])

    idepth = 4
    fac_nu = 1
    ch_lepton = 'tau'
    cross_section_model = 'ctw'
    pn_model = 'allm'
    prop_type = 'stochastic'
    stats = int(1e6)
    nu_type = 'neutrino'
    htc_mode = 'no'

    main(e_nu, angles, nu_type, cross_section_model, pn_model, idepth, ch_lepton, fac_nu, stats, prop_type, htc_mode)

    end_time = time.time()
    print(f"It took a total of {end_time-start_time:.2f} seconds to compute")
