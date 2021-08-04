#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 14:04:46 2021

@author: sam
"""

import nupyprop.data as Data
import nupyprop.geometry as Geometry
from nupyprop.propagate import run as Run

# import data as Data
# import geometry as Geometry
# from propagate import run as Run

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
        Type of output file to clean. Can be e_out or p_exit.

    Returns
    -------
    NONE
        Cleans temporary generated output files.

    '''
    if output_type == 'e_out':
        files = glob.glob("eout_*") # cleanup of Fortran e_out files
    elif output_type == 'p_exit':
        files = glob.glob("pexit_*") # cleanup of p_exit files

    for file in files:
            os.remove(file)
    return None

def init_xc(nu_type, lepton, nu_model, pn_model, prop_type):
    '''

    Parameters
    ----------
    nu_type : str
        Type of neutrino particle. Can be neutrino or anti-neutrino.
    lepton : str
        Type of lepton. Can be tau or muon.
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
    # pn_model = 'pn_' + pn_model

    if prop_type == 'stochastic':beta_type = 'cut'
    elif prop_type == 'continuous':beta_type = 'total'

    nu_xc = Data.get_xc('nu', nu_model, nu_type)

    xc_water = Data.get_xc(lepton, pn_model, 'water')
    alpha_water = Data.get_alpha(lepton, 'water')
    beta_water = Data.get_beta(lepton, pn_model, 'water', beta_type)

    xc_rock = Data.get_xc(lepton, pn_model, 'rock')
    alpha_rock = Data.get_alpha(lepton, 'rock')
    beta_rock = Data.get_beta(lepton, pn_model, 'rock', beta_type)

    return nu_xc,xc_water,xc_rock,alpha_water,alpha_rock,beta_water,beta_rock

def init_ixc(nu_type, lepton, nu_model, pn_model):
    '''

    Parameters
    ----------
    nu_type : str
        Type of neutrino particle. Can be neutrino or anti-neutrino.
    lepton : str
        Type of lepton. Can be tau or muon.
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
    # pn_model = 'pn_' + pn_model

    nu_ixc = Data.get_ixc('nu', nu_model, nu_type)

    lep_ixc_water = Data.get_ixc(lepton, pn_model, 'water')

    lep_ixc_rock = Data.get_ixc(lepton, pn_model, 'rock')

    return nu_ixc, lep_ixc_water, lep_ixc_rock


def main(E_prop, angles, nu_type, cross_section_model, pn_model, idepth, lepton, fac_nu, stats, prop_type, cdf_only):
    '''

    Parameters
    ----------
    E_prop : ndarray
        Ingoing neutrino energies, in GeV.
    angles : ndarray
        Earth emergence angles (beta), in degrees.
    nu_type : str
        Type of neutrino particle. Can be neutrino or anti-neutrino.
    cross_section_model : str
        Neutrino cross-section model.
    pn_model : str
        Photonuclear energy loss model.
    idepth : int
        Depth of water layer in km.
    lepton : str
        Type of lepton. Can be tau or muon.
    fac_nu : float
        Rescaling factor for BSM cross-sections.
    stats : float
        Statistics; no. of ingoing neutrinos.
    prop_type : str
        Energy loss propagation type. Can be stochastic or continuous.
    cdf_only : str
        If set to yes, the output file will NOT contain outgoing lepton energies.

    Returns
    -------
    NONE
        Adds exit probability & outgoing lepton energy results to output_x.h5

    '''

    make_array = lambda x : x if isinstance(x, Iterable) else np.array([x]) # to avoid errors with single or no columns

    nu_xc, xc_water, xc_rock, alpha_water, alpha_rock, beta_water, beta_rock = init_xc(nu_type, lepton, cross_section_model, pn_model, prop_type)


    nu_ixc, lep_ixc_water, lep_ixc_rock = init_ixc(nu_type, lepton, cross_section_model, pn_model)


    ithird = 0 # use dn/dy in tau to neutrino

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
    elif nu_type == 'anti-neutrino':
        nu_type = 'anu'

    start_time = time.time()
    print("The water -> rock transition occurs at %.2f degrees" % Geometry.find_interface(idepth)[0])

    for energy in E_prop:

        with open("pexit_%.2f.dat" % np.log10(energy), "a") as pexit_file:
            pexit_file.write("#\tEnergy\tAngle\tWithout Regeneration\tWith Regeneration\n")

        for angle in angles:

            xalong, cdalong = Data.get_trajs('col', angle, idepth) # initialize arrays here for each angle, to reduce a ton of overhead when tauthrulayers & regen are called

            print("Energy = %.2e, Angle = %.1f" % (energy, angle), flush=True)

            chord, water = Data.get_trajs('water', angle, idepth)
            dwater = water*rho_water # depth in water [kmwe] in last or only section
            depthE = Geometry.columndepth(angle, idepth)*1e-5 # column depth in kmwe

            no_regen, regen = Run.run_stat_single(energy, angle, nu_xc, nu_ixc, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton_int, fac_nu, stats, prop_type_int)

            prob_no_regen = no_regen/float(stats)
            prob_regen = regen/float(stats)

            with open("pexit_%.2f.dat" % np.log10(energy), "a") as pexit_file:
                pexit_file.write("%.5e\t%.5e\t%.5e\t%.5e\n" % (energy,angle,prob_no_regen,prob_regen))

            e_out = make_array(np.genfromtxt(str("eout_%.2E_%.2f.dat" % (energy, angle))))


            lep_meta = OrderedDict({'Description':'Outgoing %s energies' % lepton,
                                    'lep_energy':'Outgoing %s energy, in log_10(E) GeV'})

            lep_table = Table([e_out], names=('lep_energy',), meta=lep_meta)


            Data.add_cdf(nu_type, lepton, energy, angle, idepth, cross_section_model, pn_model, prop_type, stats, lep_table) # adds the binned cdf values for each energy and angle to output file


            file_cleaner('e_out') # remove e_out files

            if cdf_only == 'no': # adds lep_out energies to output file
                Data.add_lep_out(nu_type, lepton, energy, angle, idepth, cross_section_model, pn_model, prop_type, stats, lep_table)


        # end of for loop for angles
        p_angle, p_noregen, p_regen = np.genfromtxt("pexit_%.2f.dat" % np.log10(energy), usecols=(1,2,3), skip_header=1, unpack=True)

        pexit_meta = OrderedDict({'Description':'Exit probability for %s' % lepton,
                                  'angle':'Earth emergence angle, in degrees',
                                  'no_regen':'Exit probability without including any %s regeneration' % lepton,
                                  'regen':'Exit probability including %s regeneration' % lepton})

        pexit_table = Table([make_array(p_angle), make_array(p_noregen), make_array(p_regen)], names=('angle','no_regen','regen'), meta=pexit_meta)

        Data.add_pexit(nu_type, lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats, pexit_table) # adds p_exit results to output file

        file_cleaner('p_exit')  # remove p_exit files

    # close of for loop for energy
    end_time = time.time()
    print(f"It took {end_time-start_time:.2f} seconds to compute")
    print("Done!")

    return None

def main_htc(E_prop, angles, nu_type, cross_section_model, pn_model, idepth, lepton, fac_nu, stats, prop_type):
    '''

    Parameters
    ----------
    E_prop : ndarray
        Ingoing neutrino energies, in GeV.
    angles : ndarray
        Earth emergence angles (beta), in degrees.
    nu_type : str
        Type of neutrino particle. Can be neutrino or anti-neutrino.
    cross_section_model : str
        Neutrino cross-section model.
    pn_model : str
        Photonuclear energy loss model.
    idepth : int
        Depth of water layer in km.
    lepton : str
        Type of lepton. Can be tau or muon.
    fac_nu : float
        Rescaling factor for BSM cross-sections.
    stats : float
        Statistics; no. of ingoing neutrinos.
    prop_type : str
        Energy loss propagation type. Can be stochastic or continuous.
    cdf_only : str
        If set to yes, the output file will NOT contain outgoing lepton energies.

    Returns
    -------
    NONE
        Adds exit probability & outgoing lepton energy results to output_x.h5

    '''

    nu_xc, xc_water, xc_rock, alpha_water, alpha_rock, beta_water, beta_rock = init_xc(nu_type, lepton, cross_section_model, pn_model, prop_type)

    nu_ixc, lep_ixc_water, lep_ixc_rock = init_ixc(nu_type, lepton, cross_section_model, pn_model)

    ithird = 0 # use dn/dy in tau to neutrino

    if prop_type == 'stochastic':
        prop_type_int = 1
    else:
        prop_type_int = 2

    if lepton == 'tau':
        lepton_int = 1
    else:
        lepton_int = 2

    print("The water -> rock transition occurs at %.2f degrees" % Geometry.find_interface(idepth)[0])

    for energy in E_prop: # needs to be a 1-element array
        for angle in angles: # needs to be a 1-element array

            xalong, cdalong = Data.get_trajs('col', angle, idepth) # initialize arrays here for each angle, to reduce a ton of overhead when tauthrulayers & regen are called

            print("Energy = %.2e, Angle = %.1f" % (energy, angle), flush=True)

            chord, water = Data.get_trajs('water', angle, idepth)
            dwater = water*rho_water # depth in water [kmwe] in last or only section
            depthE = Geometry.columndepth(angle, idepth)*1e-5 # column depth in kmwe

            no_regen, regen = Run.run_stat_single(energy, angle, nu_xc, nu_ixc, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton_int, fac_nu, stats, prop_type_int)

            prob_no_regen = no_regen/float(stats)
            prob_regen = regen/float(stats)

            with open("pexit_%.2E_%.2f.dat" % (energy,angle), "w") as pexit_file:
                pexit_file.write("%.5e\t%.5e\t%.5e\t%.5e\n" % (energy,angle,prob_no_regen,prob_regen))

    print("Done!")

    return None

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
    pn_model = 'allm' # do not need to give the 'pn_' prefix if using run.py
    prop_type = 'stochastic'
    stat = int(1e7)
    nu_type = 'neutrino'
    cdf_only = 'yes'

    # nu_xc, xc_water, xc_rock, alpha_water, alpha_rock, beta_water, beta_rock = init_xc(lepton, cross_section_model, pn_model, prop_type)

    # nu_ixc, lep_ixc_water, lep_ixc_rock = init_ixc(lepton, cross_section_model, pn_model)

    # prob_dict, lep_dict = main(E_prop, angles, nu_type, cross_section_model, pn_model, idepth, lepton, fac_nu, stat, prop_type)
    main(E_prop, angles, nu_type, cross_section_model, pn_model, idepth, lepton, fac_nu, stat, prop_type, cdf_only)

    # print(prob_dict)

    # e_out_files = glob.glob("e_out_*")

    # for file in e_out_files:
    #     os.remove(file)

    # with ProcessPoolExecutor(max_workers = 4) as executor:
    #     prob_dict, e_out = executor.map(main, angles, cross_section_model, pn_model, [idepth], lepton, [fac_nu], [stat], prop_type)

    end_time = time.time()
    print(f"It took a total of {end_time-start_time:.2f} seconds to compute")
    # print("str(%.2E)" % 1e7.5)
