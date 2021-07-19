#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 16:44:13 2020

@author: sam
"""

import numpy as np
import pandas as pd
from pandas import HDFStore
from decimal import Decimal
from collections.abc import Iterable
import time
import os
# import h5py
import importlib_resources

E_nu = np.logspace(3,12,91,base=10).astype(np.float64)
E_lep = np.logspace(0,12,121,base=10).astype(np.float64)

ref = importlib_resources.files('nupyprop.datafiles') / 'lookup_tables.h5' # path for lookup_tables

def output_file(nu_type, lepton, idepth, cross_section_model, pn_model, prop_type, stats):
    '''

    Parameters
    ----------
    nu_type : str
        Type of neutrino particle. Can be nu (neutrino) or anu (anti-neutrino).
    lepton : str
        Type of lepton. Can be tau or muon.
    idepth : int
        Depth of water layer in km.
    cross_section_model : str
        Neutrino cross-section model.
    pn_model : str
        Photonuclear energy loss model.
    prop_type : str
        Type of energy loss mechanism. Can be stochastic or continuous.
    stats : float
        Statistics or number of neutrinos injected.

    Returns
    -------
    fnm : str
        Gets the name of the output file based on input parameters.

    '''
    idepth_str = str(idepth) + 'km'
    stats_str = sci_str(stats)
    fnm = "output_%s_%s_%s_%s_%s_%s_%s.h5" % (nu_type,lepton,idepth_str,cross_section_model,pn_model,prop_type,stats_str)
    return fnm

def sci_str(exp_value):
    dec = Decimal(exp_value)
    str_val = ('{:.' + str(len(dec.normalize().as_tuple().digits) - 1) + 'e}').format(dec).replace('+', '')
    return str_val

def chk_file(nu_type, lepton, idepth, cross_section_model, pn_model, prop_type, stats):
    '''

    Parameters
    ----------
    nu_type : str
        Type of neutrino particle. Can be nu (neutrino) or anu (anti-neutrino).
    lepton : str
        Type of lepton. Can be tau or muon.
    idepth : int
        Depth of water layer in km.
    cross_section_model : str
        Neutrino cross-section model.
    pn_model : str
        Photonuclear energy loss model.
    prop_type : str
        Type of energy loss mechanism. Can be stochastic or continuous.
    stats : float
        Statistics or number of neutrinos injected.

    Returns
    -------
    int
        0 is for stop execution (save previous file); 1 is for delete the old one and create a new output file; 2 is to overwrite old file (only do this if you know what you're doing or else you'll end up with mixed results!).
        Option no. 2 can be used if you need to 'add' more results for the same set of parameters or in case of abrupt code termination.

    '''
    fnm = output_file(nu_type,lepton,idepth,cross_section_model,pn_model,prop_type,stats)
    if os.path.exists(fnm):
        choice = input('There already exists a output file with these set of parameters (%s). Press \'d\' for deleting the output old file and creating a new output file, \'s\' for keeping the old output file or \'o\' for overwriting the old output file: ' % fnm)

        if choice not in {"d", "s", "o"}:
            print("Invalid option. Please enter \'d\', \'s\' or \'o\'")
            return chk_file(nu_type, lepton, idepth, cross_section_model, pn_model, prop_type, stats)
        elif choice == 's':
            return 0
        elif choice == 'd':
            os.remove(fnm)
            return 1

    else: # so basically choice = 'o'
        return 2 # output file non existant or overwrite enabled
    # return out_val


def add_trajs(type_traj, idepth, traj_array):
    '''

    Parameters
    ----------
    type_traj : str
        Type of trajectory. Can be col (for column depth) or water (for water depth).
    idepth : int
        Depth of water layer in km.
    traj_array : arr
        Trajectory numpy array.

    Returns
    -------
    None
        Adds trajectory lookup tables to lookup_tables.h5.

    '''
    with importlib_resources.as_file(ref) as lookup_tables:
        hdf = HDFStore(lookup_tables,'a')
        if type_traj == 'col':branch = 'Column_Trajectories' # sub-sub branch inside the Earth/traj_idepth branch
        elif type_traj == 'water':branch = 'Water_Trajectories'
        hdf.put('Earth/traj_%s/%s' % (str(idepth),branch), traj_array, format='t', data_columns=True)
        hdf.close()
    return print("%s lookup table successfully created for idepth = %s" % (branch,str(idepth)))

def get_trajs(type_traj, beta, idepth, out='n'):# returns {xalong:cdalong} for beta if type=col or returns chord, water for beta if type=water
    '''

    Parameters
    ----------
    type_traj : str
        Type of trajectory. Can be col (for column depth) or water (for water depth).
    beta : int
        Earth emergence angle in degrees.
    idepth : int
        Depth of water layer in km.
    out : str, optional
        'y' to write output to file. The default is 'n'.

    Returns
    -------
    ndarray
        2D array - xalong & cdalong (for col) and chord & water (for water).
        xalong - distance in water, in km.
        cdalong - column depth at xalong, in g/cm^2.
        chord - chord length, in km.
        water - final water layer distance, in km.

    '''
    if type_traj == 'col':
        with importlib_resources.as_file(ref) as lookup_tables:
            dataset = pd.read_hdf(lookup_tables,'Earth/traj_%s/Column_Trajectories' % str(idepth))
        dataset_sliced = dataset[dataset['beta']==beta]
        xalong = np.asfortranarray(dataset_sliced.xalong.T)
        cdalong = np.asfortranarray(dataset_sliced.cdalong.T)
        if out=='y':
            fnm = "%s_%s_%skm.dat" % (type_traj,beta,idepth)
            np.savetxt(fnm, np.transpose([xalong,cdalong]), header="water_distance[km]" + "\t" + "cd_at_x[g/cm^2]")
            print('Column trajectory data saved to file %s' % fnm)
        return xalong, cdalong

    elif type_traj == 'water':
        with importlib_resources.as_file(ref) as lookup_tables:
            dataset = pd.read_hdf(lookup_tables,'Earth/traj_%s/Water_Trajectories' % str(idepth))
        chord = float(dataset.chord[dataset['beta']==beta])
        water = float(dataset.water[dataset['beta']==beta])
        if out=='y':
            fnm = "%s_%s_%skm.dat" % (type_traj,beta,idepth)
            np.savetxt(fnm, np.transpose([chord,water]), header="chord_length[km]" + "\t" + "water_layer_dist[km]")
            print('Water trajectory data saved to file %s' % fnm)
        return chord, water
    return "Error in get_trajs in Data"

def add_xc(part_type, xc_obj, model, **kwargs):
    '''

    Parameters
    ----------
    part_type : str
        Neutrino or lepton? Can be nu or tau or muon.
    xc_obj : dict or ndarray
        Dictionary containing neutrino cross-section values or 1D array containing lepton cross-section values, in cm^2.
    model : str
        Neutrino cross section model.
    **kwargs
        Material: material of propagation, for leptons.

    Returns
    -------
    None
        Creates neutrino/lepton cross-section lookup entries in lookup_tables.h5.

    '''

    with importlib_resources.as_file(ref) as lookup_tables:
        hdf = HDFStore(lookup_tables,'a')
        if part_type=='nu': # here, xc_obj is a dict
            particle_type = ['nu','anu']
            for particle in particle_type:
                cc = xc_obj[particle]['cc']
                nc = xc_obj[particle]['nc']
                dframe = pd.DataFrame({'energy':E_nu, 'sigma_cc':cc, 'sigma_nc':nc})
                if 'a' in particle: # anti-neutrino group
                    hdf.put('Neutrino_Cross_Sections/anti_neutrino/xc/%s' % model, dframe, format='t', data_columns=True)
                else: # neutrino group
                    hdf.put('Neutrino_Cross_Sections/neutrino/xc/%s' % model, dframe, format='t', data_columns=True)
            hdf.close()
            return print("%s_sigma CC & NC lookup tables successfully created for %s model" % (part_type, model))
        else: # energy loss XC; here, xc_obj is an array
            material = kwargs['material']
            dframe = pd.DataFrame({'energy':E_lep, 'sigma_%s' % model:xc_obj})
            hdf.put('Energy_Loss/%s/%s/xc/%s' % (part_type,material,model), dframe, format='t', data_columns=True)
            hdf.close()
            return print("%s_sigma lookup table successfully created for %s in %s" % (part_type, model, material))
    return None

def get_xc(part_type, model, out='n', **kwargs):
    '''

    Parameters
    ----------
    part_type : str
        Neutrino or lepton? Can be nu or tau or muon.
    model : str
        Neutrino cross-section/lepton photonuclear energy loss model.
    out : str, optional
        'y' to write output to file. The default is 'n'.
    **kwargs
        particle: Type of neutrino particle. Can be neutrino or anti-neutrino.
        material: material of propagation, for leptons.

    Returns
    -------
    ndarray - 2D (neutrino) or 1D (lepton) cross-section array.
    if part_type = nu; cscc/csnc = neutrino-nucleon charged/neutral current cross sections, in cm^2.
    if part_type = tau/muon; out_arr = lepton-nucleon cross section, in cm^2.

    '''
    if part_type=='nu':
        particle = kwargs['particle']
        if particle=='anti-neutrino':particle='anti_neutrino'
        try:
            with importlib_resources.as_file(ref) as lookup_tables:
                dataset_xc = pd.read_hdf(lookup_tables,'Neutrino_Cross_Sections/%s/xc/%s' % (particle,model))
            cscc = dataset_xc.sigma_cc
            csnc = dataset_xc.sigma_nc
            out_arr = np.asarray([cscc,csnc])

            if out=='y':
                fnm = "xc_%s_%s.dat" % (particle,model)
                np.savetxt(fnm, np.transpose([E_nu,cscc,csnc]), header="Energy[GeV]" + "\t" + "sigma_cc[cm^2]" + "\t" + "sigma_nc[cm^2]")
                print('Cross-section data saved to file %s' % fnm)

            return np.asfortranarray(out_arr.T)
        except KeyError:
            model = str(input(("Error finding cross-section values for %s model, please enter a valid model name: " % model)))
            return None
    else: # energy loss; part_type == 'tau' or 'muon'
        try:
            material = kwargs['material']

            with importlib_resources.as_file(ref) as lookup_tables:
                dataset_xc = pd.read_hdf(lookup_tables,'Energy_Loss/%s/%s/xc/%s' % (part_type,material,model))
            out_arr = np.asarray(dataset_xc['sigma_%s' % model])

            if out=='y':
                fnm = "xc_%s_%s_%s.dat" % (part_type,material,model)
                np.savetxt(fnm, np.transpose([E_lep,dataset_xc['sigma_%s' % model]]), header="Energy[GeV]" + "\t" + str("sigma_%s[cm^2]") % model)
                print('Cross-section data saved to file %s' % fnm)

            return np.asfortranarray(out_arr.T)

        except KeyError:
            model = str(input(("Error finding cross-section values for %s model, please enter a valid model name: " % model)))
            return None
    return None

def add_ixc(part_type, ixc_dict, model, **kwargs):
    '''

    Parameters
    ----------
    part_type : str
        Neutrino or lepton? Can be nu or tau or muon.
    ixc_dict : dict
        Integrated cross-section CDF value dictionary for neutrinos or leptons.
    model : str
        Neutrino cross-section/lepton photonuclear energy loss model.
    **kwargs
        Material: material of propagation, for leptons.

    Returns
    -------
    None
        Creates neutrino/lepton integrated cross-section lookup entries in lookup_tables.h5.

    '''
    with importlib_resources.as_file(ref) as lookup_tables:
        hdf = HDFStore(lookup_tables,'a')
        if part_type == 'nu':
            particle_current = ['anucc','anunc','nucc','nunc']
            for particle in particle_current:
                if 'a' in particle: # anti-neutrino group
                    hdf.put('Neutrino_Cross_Sections/anti_neutrino/ixc/%s_%s' % (model,particle[3:]),ixc_dict[particle], format='t')
                else: # neutrino group
                    hdf.put('Neutrino_Cross_Sections/neutrino/ixc/%s_%s' % (model,particle[2:]),ixc_dict[particle], format='t')
            hdf.close()
            return print("%s_sigma CDF CC & NC lookup tables successfully created for %s model" % (part_type, model))

        else: # energy_loss; ixc_type == 'muon' or 'tau'
            material = kwargs['material']
            hdf.put('Energy_Loss/%s/%s/ixc/%s' % (part_type,material,model),ixc_dict, format='t')
            hdf.close()
            return print("%s_sigma CDF lookup table successfully created for %s model in %s" % (part_type, model, material))
    return None

def get_ixc(part_type, model, out='n', **kwargs):
    '''

    Parameters
    ----------
    part_type : str
        Neutrino or lepton? Can be nu or tau or muon.
    model : str
        Neutrino cross-section/lepton photonuclear energy loss model.
    out : str, optional
        'y' to write output to file. The default is 'n'.
    **kwargs
        material: material of propagation, for leptons.

    Returns
    -------
    ndarray - 2D (neutrino) or 1D (lepton) integrated cross-section CDF array.
    if part_type = nu; cscc/csnc = neutrino-nucleon charged/neutral current integrated cross section CDF values.
    if part_type = tau/muon; out_arr = lepton-nucleon integrated cross section CDF values.

    '''
    if part_type == 'nu':
        particle = kwargs['particle']
        if particle=='anti-neutrino':particle='anti_neutrino'
        try:
            with importlib_resources.as_file(ref) as lookup_tables:
                dataset_ixc_cc = pd.read_hdf(lookup_tables,'Neutrino_Cross_Sections/%s/ixc/%s_cc' % (particle,model))
                dataset_ixc_nc = pd.read_hdf(lookup_tables,'Neutrino_Cross_Sections/%s/ixc/%s_nc' % (particle,model))

            ixc_cc, ixc_nc = [], []

            for energy in range(len(E_nu)):
                ixc_cc.append(np.asarray(dataset_ixc_cc[E_nu[energy]]))
                ixc_nc.append(np.asarray(dataset_ixc_nc[E_nu[energy]]))

            ixc_cc = np.asarray(ixc_cc)
            ixc_nc = np.asarray(ixc_cc)
            # print('ixc_cc shape = ', ixc_cc.shape)

            out_arr = np.asarray([ixc_cc, ixc_nc])

            if out=='y':
                v2 = -np.linspace(0.1,3,num=30)
                yvals = 10**v2 # The integrated cross-section values should go from y = 10^(-0.1) to y = 10^(-3). This is a convention we chose to adopt. (This is also why we don't use the CTW ixc_low)
                yvals = np.insert(yvals,0,0)
                fnm = "ixc_%s_%s.dat" % (particle,model)
                np.savetxt(fnm, np.transpose([np.repeat(E_nu,31),np.tile(yvals,91),ixc_cc.flatten(),ixc_nc.flatten()]), header="Energy[GeV]" + "\t" + "y" + "\t" + "cc_cdf" + "\t" + "nc_cdf")
                print('Cross-section CDF data saved to file %s' % fnm)

            return np.asfortranarray(out_arr.T)

        except KeyError or TypeError:
            model = str(input(("Error finding integrated cross-section values for %s model, please enter a valid model name." % str(model))))
            return None
    else: # energy loss; ixc_type == 'tau' or 'muon'
        try:
            material = kwargs['material']
            with importlib_resources.as_file(ref) as lookup_tables:
                dataset_ixc = pd.read_hdf(lookup_tables,'Energy_Loss/%s/%s/ixc/%s' % (part_type,material,model))

            ixc = []
            for energy in range(len(E_lep)):
                ixc.append(np.asarray(dataset_ixc[E_lep[energy]]))

            ixc = np.asarray(ixc)

            out_arr = ixc

            if out=='y':
                v2 = -np.linspace(0.1,3,num=30)
                yvals = 10**v2 # The integrated cross-section values should go from y = 0, 10^(-0.1), 10^(-0.2).., 10^(-3). This is a convention we chose to adopt.
                yvals = np.insert(yvals,0,0) # pad with 0 at the beginning
                fnm = "ixc_%s_%s_%s.dat" % (part_type,material,model)
                np.savetxt(fnm, np.transpose([np.repeat(E_lep,31),np.tile(yvals,121),ixc.flatten()]), header="Energy[GeV]" + "\t" + "y" + "\t" + "%s_cdf" % model)
                print('Cross-section CDF data saved to file %s' % fnm)

            return np.asfortranarray(out_arr)

        except KeyError or TypeError:
            model = str(input("Error finding energy loss cross-section values for %s model, please enter a valid model name: " % str(model)))
            return None
    return None

def add_alpha(alpha, particle, material):
    '''

    Parameters
    ----------
    alpha : dict
        Ionization energy loss dictionary.
    particle : str
        Lepton. Can be tau or muon.
    material : str
        Material of propagation.

    Returns
    -------
    None
        Creates lepton ionization energy loss lookup entries in lookup_tables.h5.

    '''
    with importlib_resources.as_file(ref) as lookup_tables:
        hdf = HDFStore(lookup_tables,'a')
        alpha_df = pd.DataFrame({'energy':E_lep,'alpha':alpha})
        hdf.put('Energy_Loss/%s/%s/alpha' % (particle,material),alpha_df, format='t', data_columns=True)
        hdf.close()
        return print("%s_alpha lookup table successfully created in %s" % (particle,material))

def get_alpha(particle, material, out='n'):
    '''

    Parameters
    ----------
    particle : str
        Lepton. Can be tau or muon.
    material : str
        Material of propagation.
    out : str, optional
        'y' to write output to file. The default is 'n'.

    Returns
    -------
    ndarray
        1D array of lepton ionization energy loss, in (GeV*cm^2)/g.

    '''
    with importlib_resources.as_file(ref) as lookup_tables:
        alpha_df = pd.read_hdf(lookup_tables,'Energy_Loss/%s/%s/alpha' % (particle,material))
    alpha_arr = alpha_df.alpha
    out_arr = np.asarray(alpha_arr)

    if out=='y':
        fnm = "alpha_%s_%s.dat" % (particle,material)
        np.savetxt(fnm, np.transpose([E_lep,out_arr]), header="Energy[GeV]" + "\t" + str("alpha[(GeV*cm^2)/g]"))
        print('Alpha data saved to file %s' % fnm)

    return np.asfortranarray(out_arr.T)

def add_beta(beta_arr, particle, material, model, beta_type):
    '''

    Parameters
    ----------
    beta_arr : ndarray
        1D array containing energy loss parameter (beta), in cm^2/g.
    particle : str
        Lepton. Can be tau or muon.
    material : str
        Material of propagation.
    model : str
        Lepton energy loss model/process.
    beta_type : str
        Can be cut (for stochastic energy loss) or full (for continuous energy loss).

    Returns
    -------
    None
        Creates lepton (non-ionization) energy loss lookup entries in lookup_tables.h5.

    '''
    with importlib_resources.as_file(ref) as lookup_tables:
        hdf = HDFStore(lookup_tables,'a')
        beta_df = pd.DataFrame({'energy':E_lep, 'beta_%s' % model:beta_arr})
        hdf.put('Energy_Loss/%s/%s/beta_%s/%s' % (particle,material,beta_type,model), beta_df, format='t', data_columns=True)
        hdf.close()
        return print("%s_beta_%s lookup table successfully created in %s" % (particle,beta_type,material))

def get_beta(particle, material, model, beta_type, out='n'):
    '''

    Parameters
    ----------
    particle : str
        Lepton. Can be tau or muon.
    material : str
        Material of propagation.
    model : str
        Lepton energy loss model/process.
    beta_type : str
        Can be cut (for stochastic energy loss) or full (for continuous energy loss).
    out : str, optional
        'y' to write output to file. The default is 'n'.

    Returns
    -------
    ndarray
        1D array containing lepton energy loss model/process beta value, in cm^2/g.

    '''
    with importlib_resources.as_file(ref) as lookup_tables:
        beta_df = pd.read_hdf(lookup_tables,'Energy_Loss/%s/%s/beta_%s/%s' % (particle,material,beta_type,model))
    beta_arr = beta_df['beta_%s' % model]
    out_arr = np.asarray(beta_arr)

    if out=='y':
        fnm = "beta_%s_%s_%s_%s.dat" % (particle,material,model,beta_type)
        np.savetxt(fnm, np.transpose([E_lep,out_arr]), header="Energy[GeV]" + "\t" + str("beta[cm^2/g]"))
        print('Beta data saved to file %s' % fnm)

    return np.asfortranarray(out_arr.T)

def combine_lep(data_type, particle, material, pn_model, **kwargs):
    '''

    Parameters
    ----------
    data_type : str
        Can be xc (lepton cross-section type), beta (lepton energy loss) or ixc (lepton integrated cross-section CDFs).
    particle : str
        Lepton. Can be tau or muon.
    material : str
        Material of propagation.
    pn_model : str
        Lepton photonuclear energy loss model.
    **kwargs
        beta_type : Can be cut or full, for data_type = beta only.

    Returns
    -------
    ndarray
        Returns a combined 3D array depending on data_type.

    '''

    if data_type == 'xc':
        xc_brem = get_xc(particle, 'brem', material=material)
        xc_pair = get_xc(particle, 'pair', material=material)
        xc_pn = get_xc(particle, pn_model, material=material)
        xc_arr = np.asarray([xc_brem, xc_pair, xc_pn])
        xc = np.asfortranarray(xc_arr.T)
        return xc

    elif data_type == 'beta':
        beta_type = kwargs['beta_type']
        beta_brem = get_beta(particle, material, 'brem', beta_type)
        beta_pair = get_beta(particle, material, 'pair', beta_type)
        beta_pn = get_beta(particle, material, pn_model, beta_type)
        beta_arr = np.asarray([beta_brem, beta_pair, beta_pn])
        beta = np.asfortranarray(beta_arr.T)
        return beta

    elif data_type == 'ixc':
        ixc_brem = get_ixc(particle, 'brem', material=material)
        ixc_pair = get_ixc(particle, 'pair', material=material)
        ixc_pn = get_ixc(particle, pn_model, material=material)
        ixc_arr = np.asarray([ixc_brem, ixc_pair, ixc_pn])
        ixc = np.asfortranarray(ixc_arr.T)
        return ixc

    return None

def add_pexit(nu_type, lepton, energy, prob_dict, idepth, cross_section_model, pn_model, prop_type, stats):
    '''

    Parameters
    ----------
    nu_type : str
        Type of neutrino particle. Can be nu (neutrino) or anu (anti-neutrino).
    lepton : str
        Type of lepton. Can be tau or muon.
    energy : float
        Neutrino energy, in GeV.
    prob_dict : dict
        Exit probability dictionary.
    idepth : int
        Depth of water layer in km.
    cross_section_model : str
        Neutrino cross-section model.
    pn_model : str
        Photonuclear energy loss model.
    prop_type : str
        Type of energy loss mechanism. Can be stochastic or continuous.
    stats : float
        Statistics or number of neutrinos injected.

    Returns
    -------
    Adds exit probability table to output_x.h5.

    '''
    log_energy = np.log10(energy)
    energy_str = str(log_energy)

    make_array = lambda x : x if isinstance(x, Iterable) else np.array([x])

    angle = make_array(prob_dict['angle'])
    no_regen = make_array(prob_dict['no_regen'])
    regen = make_array(prob_dict['regen'])

    hdf = HDFStore(output_file(nu_type,lepton,idepth,cross_section_model,pn_model,prop_type,stats),'a')
    prob_df = pd.DataFrame({'angle':angle, 'no_regen':no_regen,'regen':regen})
    prob_df.set_index("angle", inplace = True)

    hdf.put('Exit_Probability/%s' % energy_str,prob_df, format='t', data_columns=True)
    hdf.close()
    return None

def get_pexit(nu_type, lepton, energy, p_type, idepth, cross_section_model, pn_model, prop_type, stats, out='n'):
    '''

    Parameters
    ----------
    nu_type : str
        Type of neutrino particle. Can be neutrino or anti-neutrino.
    lepton : str
        Type of lepton. Can be tau or muon.
    energy : float
        Neutrino energy, in GeV.
    p_type : str
        Exit probability type. Can be regen (with regeneration) or no_regen (no regeneration).
    idepth : int
        Depth of water layer in km.
    cross_section_model : str
        Neutrino cross-section model.
    pn_model : str
        Photonuclear energy loss model.
    prop_type : str
        Type of energy loss mechanism. Can be stochastic or continuous.
    stats : float
        Statistics or number of neutrinos injected.
    out : str, optional
        'y' to write output to file. The default is 'n'.

    Returns
    -------
    dict
        Returns a exit probability results dictionary with {angle:p_type}.

    '''
    log_energy = np.log10(energy)
    energy_str = str(log_energy)

    p_exit = pd.read_hdf(output_file(nu_type,lepton,idepth,cross_section_model,pn_model,prop_type,stats),'Exit_Probability/%s' % energy_str)
    angle = p_exit.index
    no_regen = p_exit.no_regen
    regen = p_exit.regen
    out_arr = np.asarray([no_regen, regen])

    if out=='y':
        fnm = "pexit_%s_%s_%s_%skm_%s_%s_%s_%s.dat" % (nu_type, lepton, energy_str, idepth, cross_section_model, pn_model, prop_type, sci_str(stats))
        np.savetxt(fnm, np.transpose([angle,no_regen,regen]), header="Earth_emergence_angle[deg]" + "\t" + "no_regen" + "\t" + "regen")
        print('Exit probability data saved to file %s' % fnm)

    return out_arr

def add_lep_out(nu_type, lepton, energy, angle, lep_dict, idepth, cross_section_model, pn_model, prop_type, stats):
    '''

    Parameters
    ----------
    nu_type : str
        Type of neutrino particle. Can be neutrino or anti-neutrino.
    lepton : str
        Type of lepton. Can be tau or muon.
    energy : float
        Neutrino energy, in GeV.
    angle : float
        Earth emergence angle (beta), in degrees.
    lep_dict : dict
        Lepton out energy dictionary, with energy in log10(GeV).
    idepth : int
        Depth of water layer in km.
    cross_section_model : str
        Neutrino cross-section model.
    pn_model : str
        Photonuclear energy loss model.
    prop_type : str
        Type of energy loss mechanism. Can be stochastic or continuous.
    stats : float
        Statistics or number of neutrinos injected.

    Returns
    -------
    None
        Adds lepton out energy tables to output_x.h5

    '''
    log_energy = np.log10(energy)
    energy_str = str(log_energy)

    lep_energies = lep_dict["lep_energy"]

    hdf = HDFStore(output_file(nu_type,lepton,idepth,cross_section_model,pn_model,prop_type,stats),'a')
    try:
        if len(lep_energies) > 1:
            lep_df = pd.DataFrame({'lep_energy':np.around(lep_energies,decimals=5)}) # so output isn't a massive file!
        elif len(lep_energies) == 1:
            lep_df = pd.DataFrame({'lep_energy':np.around(lep_energies,decimals=5)}, index=[0])
        else:
            lep_df = pd.DataFrame({'lep_energy':np.array([0])}, index=[0])
    except TypeError or ValueError:
        lep_df = pd.DataFrame({'lep_energy':np.array([0])}, index=[0])

    hdf.put('Lep_out_energies/%s/%d' % (energy_str,angle),lep_df, format='t', data_columns=True)
    hdf.close()
    return None

def get_lep_out(nu_type, lepton, energy, angle, idepth, cross_section_model, pn_model, prop_type, stats, out='n'):
    '''

    Parameters
    ----------
    nu_type : str
        Type of neutrino particle. Can be neutrino or anti-neutrino.
    lepton : str
        Type of lepton. Can be tau or muon.
    energy : float
        Neutrino energy, in GeV.
    angle : float
        Earth emergence angle (beta), in degrees.
    idepth : int
        Depth of water layer in km.
    cross_section_model : str
        Neutrino cross-section model.
    pn_model : str
        Photonuclear energy loss model.
    prop_type : str
        Type of energy loss mechanism. Can be stochastic or continuous.
    stats : float
        Statistics or number of neutrinos injected.
    out : str, optional
        'y' to write output to file. The default is 'n'.

    Returns
    -------
    out_lep : ndarray
        1D array containing lepton out energies, in GeV.

    '''
    log_energy = np.log10(energy)
    energy_str = str(log_energy)

    e_out = pd.read_hdf(output_file(nu_type,lepton,idepth,cross_section_model,pn_model,prop_type,stats),'Lep_out_energies/%s/%s' % (energy_str,angle))
    out_lep = 10**(np.asarray(e_out.lep_energy)) # changed 13/7/21

    if out=='y':
        fnm = "lep_out_%s_%s_%s_%sdeg_%skm_%s_%s_%s_%s.dat" % (nu_type, lepton, energy_str, angle, idepth, cross_section_model, pn_model, prop_type, sci_str(stats))
        np.savetxt(fnm, np.transpose([out_lep]), header="lep_out_energy[GeV]")
        print('Lepton outgoing energy data saved to file %s' % fnm)

    return out_lep

# def add_pexit_manual(nu_type, energy, angles, idepth, cross_section_model, pn_model, prop_type, stats): # manual will only work for regen (so basically, for muons)

#     log_energy = np.log10(energy)
#     lepton = 'muon'
#     energy_str = str(log_energy)

#     angle, no_regen, regen = [],[],[]
#     for angle in angles:
#         e_out = pd.read_hdf(output_file(nu_type,lepton,idepth,cross_section_model,pn_model,prop_type,stats),'Lep_out_energies/%s/%s' % (energy_str,angle))
#         e_out_len =  len(e_out.lep_energy)
#         if e_out_len == 1:
#             if np.asarray(e_out.lep_energy)[0] == 0:
#                 p_exit = 0.0
#             else:
#                 p_exit = e_out_len/stats
#         else:
#             p_exit = e_out_len/stats
#         angle.append(angle)
#         no_regen.append(p_exit)
#         regen.append(p_exit)

#     hdf = HDFStore(output_file(nu_type,lepton,idepth,cross_section_model,pn_model,prop_type,stats),'a')
#     prob_df = pd.DataFrame({'angle':angle, 'no_regen':no_regen,'regen':regen})
#     prob_df.set_index("angle", inplace = True)

#     hdf.put('Exit_Probability/%s' % energy_str,prob_df, format='t', data_columns=True)
#     hdf.close()
#     return None

def add_cdf(nu_type, lepton, energy, angle, lep_dict, idepth, cross_section_model, pn_model, prop_type, stats):
    '''

    Parameters
    ----------
    nu_type : str
        Type of neutrino particle. Can be neutrino or anti-neutrino.
    lepton : str
        Type of lepton. Can be tau or muon.
    energy : float
        Neutrino energy, in GeV.
    angle : float
        Earth emergence angle (beta), in degrees.
    lep_dict : dict
        Outgoing lepton energy dictionary with {"lep_energy":energy_arr}.
    idepth : int
        Depth of water layer in km.
    cross_section_model : str
        Neutrino cross-section model.
    pn_model : str
        Photonuclear energy loss model.
    prop_type : str
        Type of energy loss mechanism. Can be stochastic or continuous.
    stats : float
        Statistics or number of neutrinos injected.

    Returns
    -------
    NONE
        Adds outgoing lepton energy CDF values to output_x.h5.

    '''
    log_energy = np.log10(energy)
    energy_str = str(log_energy)

    hdf = HDFStore(output_file(nu_type,lepton,idepth,cross_section_model,pn_model,prop_type,stats),'a')

    bins = np.logspace(-5,0,51) # Default binning for use with nuSpaceSim. Change if different binning required.
    lep_out = 10**lep_dict['lep_energy'] # because lep_out energies are in log10(GeV)
    z = lep_out/energy # z = E_lep/E_nu
    binned_z = np.digitize(z, bins)
    bin_counts = np.bincount(binned_z)
    if len(bin_counts) < len(bins):
        zeros_z = np.zeros(len(bins) - len(bin_counts))
        binned_z_fixed = np.concatenate((bin_counts, zeros_z))
    else:
        binned_z_fixed = bin_counts
    z_cumsum = np.cumsum(binned_z_fixed)
    z_cdf = z_cumsum/z_cumsum[-1]
    z_cdf[0] = 0
    z_cdf_df = pd.DataFrame({'z':bins, 'CDF':np.around(z_cdf,decimals=8)})
    hdf.put('Lep_out_cdf/%s/%d' % (energy_str,angle),z_cdf_df, format='t', data_columns=True)
    hdf.close()
    return None

def get_cdf(nu_type, lepton, energy, angle, idepth, cross_section_model, pn_model, prop_type, stats, out='n'):
    '''

    Parameters
    ----------
    nu_type : str
        Type of neutrino particle. Can be neutrino or anti-neutrino.
    lepton : str
        Type of lepton. Can be tau or muon.
    energy : float
        Neutrino energy, in GeV.
    angle : float
        Earth emergence angle (beta), in degrees.
    idepth : int
        Depth of water layer in km.
    cross_section_model : str
        Neutrino cross-section model.
    pn_model : str
        Photonuclear energy loss model.
    prop_type : str
        Type of energy loss mechanism. Can be stochastic or continuous.
    stats : float
        Statistics or number of neutrinos injected.
    out : str, optional
        'y' to write output to file. The default is 'n'.

    Returns
    -------
    cdf : dict
        Dictionary of lepton outgoing value CDFs binned in z=E_nu/E_lep.

    '''

    log_energy = np.log10(energy)
    energy_str = str(log_energy)

    df = pd.read_hdf(output_file(nu_type,lepton,idepth,cross_section_model,pn_model,prop_type,stats),'Lep_out_cdf/%s/%s' % (energy_str,angle))

    cdf = dict(zip(df.z, df.cdf))
    # cdf = np.asarray(df.cdf) # uncomment this line if you want a cdf array instead of a dict

    if out=='y':
        fnm = "cdf_%s_%s_%s_%sdeg_%skm_%s_%s_%s_%s.dat" % (nu_type, lepton, energy_str, angle, idepth, cross_section_model, pn_model, prop_type, sci_str(stats))
        np.savetxt(fnm, np.transpose([df.z, df.cdf]), header="z=E_lep/E_nu" + "\t" + "cdf")
        print('Lepton outgoing energy CDF data saved to file %s' % fnm)
    return cdf

# def add_header(nu_type, lepton, idepth, nu_cs, lep_pn, loss_type, stats):
#     f = h5py.File('output_nu_tau_4km_ct18nlo_allm_stochastic_1e8.h5','a')
#     if nu_type=='nu':f.attrs['Matter of Neutrino'] = "Neutrino"
#     elif nu_type=='anu':f.attrs['Matter of Neutrino'] = "Anti-Neutrino"
#     f.attrs['lepton'] = lepton
#     f.attrs['water_depth'] = str(idepth) + "km"
#     f.attrs['neutrino_model'] = nu_cs
#     f.attrs['photonuclear_model'] = lep_pn
#     f.attrs['loss_type'] = loss_type
#     f.attrs['neutrinos_in'] = stats
#     cdf_bin_size = 51 # this is fixed from the variable 'bins' in function add_cdf(...)
#     f.attrs['cdf_bin_size'] = cdf_bin_size
#     f.close()
#     return None

# =============================================================================
# Test
# =============================================================================
if __name__ == "__main__":
    # arr = get_beta('tau', 'rock', 'total', 'allm', True)
    # add_pexit_manual(1e9, np.arange(1,36), 1e8)
    pass
    # pexit = get_pexit('tau', 1e7, p_type='regen', loss_type='stochastic')
    # add_cdf('tau', 'full')
    # write_xc('nu', 'allm', particle='neutrino')
    # write_xc('nu', 'allm', particle='anti-neutrino')
    # write_xc('tau', 'pn_bb', material='rock')
    # get_ixc('tau', 'pn_allm', out='y', material='rock')