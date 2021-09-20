#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 16:44:13 2020

@author: sam
"""

import numpy as np
from decimal import Decimal
import importlib_resources
from astropy.table import Table
from astropy.io import ascii
from collections import OrderedDict
import glob
import h5py
import os
from collections.abc import Iterable
from scipy.interpolate import interpn

E_nu = np.logspace(3,12,91,base=10).astype(np.float64)
E_lep = np.logspace(0,12,121,base=10).astype(np.float64)

nu_models = ['allm', 'bdhm', 'ct18nlo', 'nct15']
pn_models = ['allm', 'bb']

ref = importlib_resources.files('nupyprop.datafiles') / 'lookup_tables.h5' # path for lookup_tables

class ModelError(Exception):
    """Exception raised for errors in the input custom model file.

    Attributes:
        fnm -- input filename which caused the error
        message -- explanation of the error
    """

    def __init__(self, fnm, message="This is either an incorrectly formatted model file or file not found"):
        self.fnm = fnm
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.fnm} -> {self.message}'

def patch_for_astropy(arr):
    """makes a patch for astropy to avoid errors with '0 len' arrays

    Args:
        arr (ndarray): 1D array containing data

    Returns:
        ndarray: astropy patched array containing the same data
    """
    try:
        len(arr)
    except TypeError:
        arr = arr.reshape(1)
    if arr.size == 0:arr = np.zeros(1) # set empty array to 0
    return arr

def sci_str(exp_value):
    """converts exponential value to scientific string format

    Args:
        exp_value (float): exponential value to be converted

    Returns:
        str: string equivalent of input with some custom formatting
    """
    dec = Decimal(exp_value)
    str_val = ('{:.' + str(len(dec.normalize().as_tuple().digits) - 1) + 'e}').format(dec).replace('+', '')
    return str_val

def get_custom_path(data_type, part_type, model, arg):
    """gets the path of the model file to be used for loading data from custom model files

    Args:
        data_type (str): type of data; can be xc for cross-section values, ixc for cross-section CDF
            values and beta for energy loss parameter values
        part_type (str): type of particle; can be nu for neutrinos or tau for tau leptons or muon for muons
        model (str): name of the model
        arg (str): type of neutrino for part_type=nu or the propagation material for part_type=tau and part_type=muon

    Raises:
        ModelError: either an incorrectly formatted model file or file not found
        ModelError: either an incorrectly formatted model file or file not found

    Returns:
        PosixPath: path of the csutom mode file
    """
    if part_type == 'nu':
        nu_type = arg
        fnm = data_type + '_%s_%s.ecsv' % (nu_type,model)
        file = importlib_resources.files('nupyprop.models') / fnm
        if not os.path.exists(file):
            raise ModelError(fnm)
        return file

    elif part_type == 'tau' or part_type == 'muon':
        material = arg
        fnm = data_type + '_%s_pn_%s_%s.ecsv' % (part_type,model,material)
        file = importlib_resources.files('nupyprop.models') / fnm
        if not os.path.exists(file):
            raise ModelError(fnm)
        return file

def output_file(nu_type, ch_lepton, idepth, cross_section_model, pn_model, prop_type, stats, arg=None):
    """gets the name of the output file based on input args

    Args:
        nu_type (str): type of neutrino particle; can be neutrino or anti_neutrino
        ch_lepton (str): type of charged lepton; can be tau or muon
        idepth (int): depth of water layer, in km
        cross_section_model (str): neutrino cross-section model
        pn_model (str): photonuclear energy loss model
        prop_type (str): type of energy loss mechanism; can be stochastic or continuous
        stats (int): statistics or number of neutrinos injected
        arg (str, optional): additional arguments at the end of the file name. Defaults to None.

    Returns:
        str: name of the output file
    """
    idepth_str = str(idepth) + 'km'
    stats_str = sci_str(stats)
    pn_model = pn_model.replace("pn_","")
    if nu_type=='neutrino':nu_type='nu' # shorten output filename
    elif nu_type=='anti_neutrino':nu_type='anu' # shorten output filename
    if arg is None:
        fnm = "output_%s_%s_%s_%s_%s_%s_%s.h5" % (nu_type,ch_lepton,idepth_str,cross_section_model,pn_model,prop_type,stats_str)
    else: # if a custom argument is present, it should be in the end of the filename
        fnm = "output_%s_%s_%s_%s_%s_%s_%s_%s.h5" % (nu_type,ch_lepton,idepth_str,cross_section_model,pn_model,prop_type,stats_str,arg)
    return fnm

def add_trajs(type_traj, idepth, traj_table):
    """adds trajectory values to lookup_tables.h5

    Args:
        type_traj (str): type of trajectory; can be col (for column depth) or water (for water depth)
        idepth (int): depth of water layer, in km
        traj_table (`~astropy.table.Table`): astropy table containing the trajectories
            Each table will need the following columns:
            - ``beta``: Earth emergence angle, in degrees [``beta``]
            - ``xalong``: Distance in water, in km [`xalong``]
            - ``cdalong``: Column depth at xalong, in g/cm^2 [``cdalong``]

    Returns:
        None
    """
    if type_traj == 'col':branch = 'Column_Trajectories' # sub-sub branch inside the Earth/traj_idepth branch
    elif type_traj == 'water':branch = 'Water_Trajectories'

    with importlib_resources.as_file(ref) as lookup_tables:
        traj_table.write(lookup_tables, path='Earth/%s/%skm' % (branch,str(idepth)), append=True, overwrite=True)
    return print("%s lookup table successfully created for idepth = %s" % (branch,str(idepth)))

def get_trajs(type_traj, angle, idepth, out=False):
    """get trajectory values

    Args:
        type_traj (str): type of trajectory; can be col (for column depth) or water (for water depth)
        angle (float): earth emergence angle, in degrees
        idepth (int): depth of water layer, in km
        out (bool, optional): saves the data as an ASCII ecsv file if set to True; returns the array
            value if set to False. Defaults to False.

    Returns:
        None/tuple: None if out=True otherwise tuple containing:
            xalong (ndarray): 1D Fortran array of shape (100,) containing distance in water, in km
            cdalong (ndarray): 1D Fortran array of shape (100,) containing column depth at xalong, in g/cm^2
            or
            chord (float): Chord length, in km
            water (float): Final water layer distance, in km
    """
    if type_traj == 'col':
        with importlib_resources.as_file(ref) as lookup_tables:
            traj_table = Table.read(lookup_tables,path='Earth/Column_Trajectories/%skm' % str(idepth))

        sliced_table = traj_table[traj_table['beta']==angle]
        xalong = np.asfortranarray(sliced_table['xalong'].T)
        cdalong = np.asfortranarray(sliced_table['cdalong'].T)

        if out:
            fnm = "%s_%.2fdeg_%skm.ecsv" % (type_traj,angle,idepth)
            ascii.write(traj_table, fnm, format='ecsv', fast_writer=True, overwrite=True)
            return print('Column trajectory data saved to file %s' % fnm)
        return xalong, cdalong

    elif type_traj == 'water':
        with importlib_resources.as_file(ref) as lookup_tables:
            traj_table = Table.read(lookup_tables,path='Earth/Water_Trajectories/%skm' % str(idepth))

        chord = float(traj_table['chord'][traj_table['beta']==angle])
        water = float(traj_table['water'][traj_table['beta']==angle])

        if out:
            fnm = "%s_%.2fdeg_%skm.ecsv" % (type_traj,angle,idepth)
            ascii.write(traj_table, fnm, format='ecsv', fast_writer=True, overwrite=True)
            return print('Water trajectory data saved to file %s' % fnm)
        return chord, water
    return "Error in get_trajs in Data"

def add_xc(part_type, xc_table, arg):
    """adds cross-section values to lookup_tables.h5

    Args:
        part_type (str): type of particle; can be nu for neutrinos or tau for tau leptons or muon for muons
        xc_table (`~astropy.table.Table`): astropy table containing cross-section values
            Each table will need the following columns:
            - if part_type=nu:
                - ``energy``: Neutrino (or anti-neutrino) energy, in GeV [``energy``]
                - ``sigma_cc_x``: Charged current cross-section for x model, in cm^2; x is the name
                    of the model [`sigma_cc_x``]
                - ``sigma_nc_x``: Neutral current cross-section for x model, in cm^2; x is the name
                    of the model [``sigma_nc_x``]
            - if part_type=tau or part_type=muon:
                - ``energy``: Charged lepton energy, in GeV [``energy``]
                - ``sigma_brem``: N_A/A * cross-section for bremmstrahlung in material, in cm^2/g [`sigma_brem``]
                - ``sigma_pair``: N_A/A * cross-section for pair production in material, in cm^2/g [``sigma_pair``]
                - ``sigma_pn_x``: N_A/A * cross-section for x model in material, in cm^2/g; x is the name of the
                    photonuclear energy loss model [``sigma_pn_x``]
        arg (str): type of neutrino for part_type=nu or the propagation material for part_type=tau and part_type=muon

    Returns:
        None
    """
    if part_type=='nu':
        nu_type = arg
        with importlib_resources.as_file(ref) as lookup_tables:
            xc_table.write(lookup_tables, path='Neutrinos/%s/xc' % nu_type, append=True, overwrite=True)
        return print("%s_sigma CC & NC lookup tables successfully created" % part_type)
    else: # charged lepton energy loss XC
        material = arg
        with importlib_resources.as_file(ref) as lookup_tables:
            xc_table.write(lookup_tables, path='Charged_Leptons/%s/%s/xc' % (part_type,material), append=True, overwrite=True)
        return print("%s_sigma lookup table successfully created in %s" % (part_type, material))
    return None

def get_xc(part_type, model, arg, out=False):
    """get cross-section values; works with custom PN models

    Args:
        part_type (str): type of particle; can be nu for neutrinos or tau for tau leptons or muon for muons
        model (str): name of the model
        arg (str): type of neutrino for part_type=nu or the propagation material for part_type=tau and part_type=muon
        out (bool, optional): saves the data as an ASCII ecsv file if set to True; returns the array
            value if set to False. Defaults to False.

    Returns:
        None/ndarray: None if out=True otherwise otherwise:
            2D Fortran array of shape (91,2) containing charged current and neutral
            current neutrino (or anti-neutrino) cross-section values, in cm^2, or
            3D Fortran array of shape (121,3) containing Bremmstrahlung, pair-production
            and photonuclear cross-section values for charged lepton, all multiplied by N_A/A,
            in cm^2/g
    """
    if part_type=='nu':
        nu_type = arg
        if model in nu_models: # default nu_xc model selection
            with importlib_resources.as_file(ref) as lookup_tables:
                xc_table = Table.read(lookup_tables,path='Neutrinos/%s/xc' % nu_type)
        else: # custom nu_xc model selection
            file = get_custom_path('xc', part_type, model, nu_type)
            with importlib_resources.as_file(file) as lookup_table:
                xc_table = Table.read(lookup_table, format='ascii.ecsv')

        cscc = xc_table['sigma_cc_%s' % model]
        csnc = xc_table['sigma_nc_%s' % model]
        xc_arr = np.asarray([cscc,csnc])

        if out:
            fnm = "xc_%s_%s.ecsv" % (nu_type,model)
            xc_meta = OrderedDict({'Description':'%s-nucleon cross-section values for %s' % (nu_type.capitalize(),str.upper(model)),
                                   'energy':'%s energy, in GeV' % nu_type.capitalize(),
                                   'sigma_cc_%s' % model:'Charged current cross-section for %s, in cm^2' % str.upper(model),
                                   'sigma_nc_%s' % model:'Neutral current cross-section for %s, in cm^2' % str.upper(model)})
            out_table = Table([xc_table['energy'], cscc, csnc], names=('energy','sigma_cc_%s' % model,'sigma_nc_%s' % model), meta=xc_meta)
            ascii.write(out_table, fnm, format='ecsv', fast_writer=False, overwrite=True)
            return print('%s cross-section data saved to file %s' % (nu_type,fnm))

        return np.asfortranarray(xc_arr.T)

    else: # energy loss; part_type == 'tau' or 'muon'
        material = arg
        with importlib_resources.as_file(ref) as lookup_tables:
            xc_table = Table.read(lookup_tables,path='Charged_Leptons/%s/%s/xc' % (part_type,material))

        xc_brem = np.asarray(xc_table['sigma_brem']) # brem is not a custom model
        xc_pair = np.asarray(xc_table['sigma_pair']) # pair is not a custom model

        if model in pn_models: # default PN model selection

            xc_pn = np.asarray(xc_table['sigma_pn_%s' % model])

        else: # custom PN model
            with importlib_resources.as_file(get_custom_path('xc',part_type,model,material)) as lookup_table:
                pn_table = Table.read(lookup_table,format='ascii.ecsv')
            xc_pn = np.asarray(pn_table['sigma_pn_%s' % model])

        xc_arr = np.asarray([xc_brem, xc_pair, xc_pn])

        if out:
            fnm = "xc_%s_pn_%s_%s.ecsv" % (part_type,model,material)
            xc_meta = OrderedDict({'Description':'%s-nucleon cross-section values for %s in %s' % (part_type.capitalize(),str.upper(model),material),
                                   'energy':'%s energy, in GeV' % part_type.capitalize(),
                                   'sigma_brem':'N_A/A * cross-section for bremmstrahlung in %s, in cm^2/g' % material,
                                   'sigma_pair':'N_A/A * cross-section for pair production in %s, in cm^2/g' % material,
                                   'sigma_pn_%s' % model:'N_A/A * cross-section for PN_%s in %s, in cm^2/g' % (str.upper(model),material)})
            out_table = Table([xc_table['energy'], xc_brem, xc_pair, xc_pn], names=('energy','sigma_brem','sigma_pair','sigma_pn_%s' % model), meta=xc_meta)
            ascii.write(out_table, fnm, format='ecsv', fast_writer=False, overwrite=True)
            return print('%s cross-section data saved to file %s' % (part_type,fnm))

        return np.asfortranarray(xc_arr.T)

def add_ixc(part_type, ixc_table, arg):
    """adds cross-section CDF values to lookup_tables.h5

    Args:
        part_type (str): type of particle; can be nu for neutrinos or tau for tau leptons or muon for muons
        xc_table (`~astropy.table.Table`): astropy table containing cross-section CDF values
            Each table will need the following columns:
            - if part_type=nu:
                - ``energy``: Neutrino (or anti-neutrino) energy, in GeV [``energy``]
                - ``y``: Inelasticity, y = (E_init-E_final)/E_initial [``y``]
                - ``cc_cdf_x``: Charged current cross-section CDF values for x model; x is the name
                    of the model [`cc_cdf_x``]
                - ``nc_cdf_x``: Neutral current cross-section CDF values for x model; x is the name
                    of the model [``nc_cdf_x``]
            - if part_type=tau or part_type=muon:
                - ``energy``: Charged lepton energy, in GeV [``energy``]
                - ``y``: Inelasticity, y = (E_init-E_final)/E_initial [``y``]
                - ``cdf_brem``: Cross-section CDF values for Bremmstrahlung in material, in cm^2 [`cdf_brem``]
                - ``cdf_pair``: Cross-section CDF values for pair-production in material, in cm^2 [`cdf_pair``]
                - ``cdf_pn_x``: Cross-section CDF values for x model in material, in cm^2/g; x is the name
                    of the photonuclear energy loss model [`cdf_pn_x``]
        arg (str): type of neutrino for part_type=nu or the propagation material for part_type=tau and part_type=muon

    Returns:
        None
    """
    if part_type == 'nu':
        nu_type = arg
        with importlib_resources.as_file(ref) as lookup_tables:
            ixc_table.write(lookup_tables, path='Neutrinos/%s/ixc' % nu_type, append=True, overwrite=True)

        return print("%s_sigma CDF CC & NC lookup tables successfully created for" % nu_type)

    else: # energy_loss; part_type == 'muon' or 'tau'
        material = arg

        with importlib_resources.as_file(ref) as lookup_tables:
            ixc_table.write(lookup_tables, path='Charged_Leptons/%s/%s/ixc' % (part_type,material), append=True, overwrite=True)

        return print("%s_sigma CDF lookup table successfully created in %s" % (part_type, material))
    return None

def get_ixc(part_type, model, arg, out=False):
    """get integrated cross-section CDF values; works with custom PN models

    Args:
        part_type (str): type of particle; can be nu for neutrinos or tau for tau leptons or muon for muons
        model (str): name of the model
        arg (str): type of neutrino for part_type=nu or the propagation material for part_type=tau and part_type=muon
        out (bool, optional): saves the data as an ASCII ecsv file if set to True; returns the array
            value if set to False. Defaults to False.

    Returns:
        None/ndarray: None if out=True otherwise otherwise:
            3D Fortran array of shape (31,91,2) containing integrated charged current and neutral
            current neutrino (or anti-neutrino) cross-section CDF values, or
            3D Fortran array of shape (31,121,3) containing integrated Bremmstrahlung, pair-production
            and photonuclear cross-section CDF values for charged lepton
    """
    if part_type == 'nu':
        nu_type = arg
        if model in nu_models: # default nu_ixc model selection
            with importlib_resources.as_file(ref) as lookup_tables:
                ixc_table = Table.read(lookup_tables,path='Neutrinos/%s/ixc' % nu_type)
        else: # custom nu_ixc model selection
            file = get_custom_path('ixc', part_type, model, nu_type)
            with importlib_resources.as_file(file) as lookup_table:
               ixc_table = Table.read(lookup_table, format='ascii.ecsv')

        ixc_cc = np.asarray([ixc_table['cc_cdf_%s' % model][ixc_table['energy']==i] for i in E_nu])
        ixc_nc = np.asarray([ixc_table['nc_cdf_%s' % model][ixc_table['energy']==i] for i in E_nu])

        ixc_arr = np.asarray([ixc_cc, ixc_nc])

        if out:
            fnm = "ixc_%s_%s.ecsv" % (nu_type,model)
            ixc_meta = OrderedDict({'Description':'%s-nucleon cross-section CDF values for %s' % (nu_type.capitalize(),str.upper(model)),
                                    'energy':'%s energy, in GeV' % nu_type.capitalize(),
                                    'y':'Inelasticity; y = (E_initial-E_final)/E_initial',
                                    'cc_cdf_%s' % model:'Charged current cross-section CDF values for %s' % str.upper(model),
                                    'nc_cdf_%s' % model:'Neutral current cross-section CDF values for %s' % str.upper(model)})
            out_table = Table([ixc_table['energy'], ixc_table['y'], ixc_table['cc_cdf_%s' % model], ixc_table['nc_cdf_%s' % model]], meta=ixc_meta)
            ascii.write(out_table, fnm, format='ecsv', fast_writer=False, overwrite=True)
            return print('%s cross-section CDF data saved to file %s' % (nu_type,fnm))

        return np.asfortranarray(ixc_arr.T)


    else: # energy loss; ixc_type == 'tau' or 'muon'
        material = arg

        with importlib_resources.as_file(ref) as lookup_tables:
            ixc_table = Table.read(lookup_tables,path='Charged_Leptons/%s/%s/ixc' % (part_type,material))

        ixc_brem = np.asarray([ixc_table['cdf_brem'][ixc_table['energy']==i] for i in E_lep]) # brem is not a custom model
        ixc_pair = np.asarray([ixc_table['cdf_pair'][ixc_table['energy']==i] for i in E_lep]) # pair is not a custom model

        if model in pn_models: # default PN model selection

            ixc_pn = np.asarray([ixc_table['cdf_pn_%s' % model][ixc_table['energy']==i] for i in E_lep])

        else: # custom PN model
            with importlib_resources.as_file(get_custom_path('ixc',part_type,model,material)) as lookup_table:
                pn_table = Table.read(lookup_table,format='ascii.ecsv')
            ixc_pn = np.asarray([pn_table['cdf_pn_%s' % model][ixc_table['energy']==i] for i in E_lep])

        ixc_arr = np.asarray([ixc_brem, ixc_pair, ixc_pn])

        if out:
            fnm = "ixc_%s_pn_%s_%s.ecsv" % (part_type,model,material)
            ixc_meta = OrderedDict({'Description':'%s-nucleon cross-section CDF values for PN_%s in %s' % (part_type.capitalize(),str.upper(model),material),
                                    'energy':'%s energy, in GeV' % part_type.capitalize(),
                                    'y':'Inelasticity; y = (E_initial-E_final)/E_initial',
                                    'cdf_brem':'Cross-section CDF values for bremmstrahlung in %s' % material,
                                    'cdf_pair':'Cross-section CDF values for pair production in %s' % material,
                                    'cdf_pn_%s' % model:'Cross-section CDF values for PN_%s in %s' % (str.upper(model),material)})
            out_table = Table([ixc_table['energy'], ixc_table['y'], ixc_brem.flatten(), ixc_pair.flatten(), ixc_pn.flatten()], names=('energy','y','cdf_brem','cdf_pair','cdf_pn_%s' % model), meta=ixc_meta)
            ascii.write(out_table, fnm, format='ecsv', fast_writer=True, overwrite=True)
            return print('%s cross-section CDF data saved to file %s' % (part_type,fnm))

        return np.asfortranarray(ixc_arr.T)


def add_alpha(ch_lepton, material, alpha_table):
    """adds charged lepton ionization energy loss values to lookup_tables.h5

    Args:
        ch_lepton (str): type of charged lepton; can be tau or muon
        material (str): material of propagation for charged lepton
        alpha_table (`~astropy.table.Table`): ionization energy loss values
            Each table will need the following columns:
            - ``energy``: Charged lepton energy, in GeV [``energy``]
            - ``alpha``: Ionization energy loss value in material, in (GeV*cm^2)/g [``alpha``]

    Returns:
        None
    """
    with importlib_resources.as_file(ref) as lookup_tables:
        alpha_table.write(lookup_tables, path='Charged_Leptons/%s/%s/alpha' % (ch_lepton,material), append=True, overwrite=True)
    return print("%s_alpha lookup table successfully created for %s" % (ch_lepton,material))

def get_alpha(ch_lepton, material, out=False):
    """get charged lepton ionization energy loss values

    Args:
        ch_lepton (str): type of charged lepton; can be tau or muon
        material (str): material of propagation for charged lepton
        out (bool, optional): saves the data as an ASCII ecsv file if set to True; returns the array
            value if set to False. Defaults to False.

    Returns:
        None/ndarray: None if out=True otherwise otherwise:
            1D Fortran array of shape (121,) containing charged lepton ionization energy loss values in material,
            in (GeV*cm^2)/g
    """
    with importlib_resources.as_file(ref) as lookup_tables:
        alpha_table = Table.read(lookup_tables,path='Charged_Leptons/%s/%s/alpha' % (ch_lepton,material))
    alpha_arr = alpha_table['alpha']

    if out:
        fnm = "alpha_%s_%s.ecsv" % (ch_lepton,material)
        ascii.write(alpha_table, fnm, format='ecsv', fast_writer=True, overwrite=True)
        return print('Alpha data saved to file %s' % fnm)

    return np.asfortranarray(alpha_arr.T)

def add_beta(ch_lepton, material, beta_table):
    """adds charged lepton energy loss parameter values to lookup_tables.h5

    Args:
        ch_lepton (str): type of charged lepton; can be tau or muon
        material (str): material of propagation for charged lepton
        beta_table (`~astropy.table.Table`): energy loss parameter values
            Each table will need the following columns:
            - ``energy``: Charged lepton energy, in GeV [``energy``]
            - ``beta_brem_cut``: Energy loss parameter value for Bremmstrahlung in material (cut value),
            in cm^2/g [``beta_brem_cut``]
            - ``beta_brem_total``: Energy loss parameter value for Bremmstrahlung in material (total value),
            in cm^2/g [``beta_brem_total``]
            - ``beta_pair_cut``: Energy loss parameter value for pair-production in material (cut value),
            in cm^2/g [``beta_pair_cut``]
            - ``beta_pair_total``: Energy loss parameter value for pair-production in material (total value),
            in cm^2/g [``beta_pair_total``]
            - ``beta_pn_x_cut``: Energy loss parameter value in material (cut value), in cm^2/g;
            x is the name of the photonuclear energy loss model  [``beta_pn_x_cut``]
            - ``beta_pn_x_total``: Energy loss parameter value in material (total value), in cm^2/g;
            x is the name of the photonuclear energy loss model  [``beta_pn_x_total``]

    Returns:
        None
    """
    with importlib_resources.as_file(ref) as lookup_tables:
        beta_table.write(lookup_tables, path='Charged_Leptons/%s/%s/beta' % (ch_lepton,material), append=True, overwrite=True)
    return print("%s_beta in %s lookup table successfully created" % (ch_lepton,material))

def get_beta(ch_lepton, model, material, arg, out=False):
    """get charged lepton ionization energy loss values; works with custom PN models

    Args:
        ch_lepton (str): type of charged lepton; can be tau or muon
        material (str): material of propagation for charged lepton
        arg (str): type of beta parameter; can be cut or total
        out (bool, optional): saves the data as an ASCII ecsv file if set to True; returns the array
            value if set to False. Defaults to False.

    Returns:
        None/ndarray: None if out=True otherwise otherwise:
            2D Fortran array of shape (121,3) containing charged lepton energy loss parameter values in material,
            for Bremmstrahlung, pair-production and photonuclear energy loss model, in cm^2/g
    """
    beta_type = arg

    if beta_type=='cut':beta_type_str = 'y_max=1e-3'
    elif beta_type=='total':beta_type_str = 'y_max'

    with importlib_resources.as_file(ref) as lookup_tables:
        beta_table = Table.read(lookup_tables,path='Charged_Leptons/%s/%s/beta' % (ch_lepton,material))

    beta_brem = beta_table['beta_brem_%s' % beta_type] # brem is not a custom model
    beta_pair = beta_table['beta_pair_%s' % beta_type] # pair is not a custom model

    if model in pn_models: # default PN model selection
        beta_pn = beta_table['beta_pn_%s_%s' % (model,beta_type)]

    else: # custom PN model
        with importlib_resources.as_file(get_custom_path('beta',ch_lepton,model,material)) as lookup_table:
            pn_table = Table.read(lookup_table,format='ascii.ecsv')
        beta_pn = pn_table['beta_pn_%s_%s' % (model,beta_type)]

    beta_arr = np.asarray([beta_brem, beta_pair, beta_pn])

    if out:
        fnm = "beta_%s_pn_%s_%s_%s.ecsv" % (ch_lepton,model,material,beta_type)
        beta_meta = OrderedDict({'Description':'Model/parameterization dependent photonuclear energy loss lookup table for %s in %s' % (ch_lepton.capitalize(),material),
                                  'energy':'%s energy, in GeV' % ch_lepton.capitalize(),
                                  'beta_pn_%s_%s' % (model,beta_type):'Photonuclear %s energy loss model beta values integrated from y_min to %s, in cm^2/g' % (str.upper(model),beta_type_str)})
        out_table = Table([beta_table['energy'], beta_pn], names=('energy','beta_pn_%s_%s' % (model,beta_type)), meta=beta_meta)
        ascii.write(out_table, fnm, format='ecsv', fast_writer=False, overwrite=True)
        return print('Beta data saved to file %s' % fnm)

    return np.asfortranarray(beta_arr.T)


def add_pexit(nu_type, ch_lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats, pexit_table, arg=None):
    """adds exit probability values to output file

    Args:
        nu_type (str): type of neutrino particle; can be neutrino or anti_neutrino
        ch_lepton (str): type of charged lepton; can be tau or muon
        energy (float): ingoing neutrino (or anti-neutrino) energy, in GeV
        idepth (int): depth of water layer, in km
        cross_section_model (str): neutrino cross-section model
        pn_model (str): photonuclear energy loss model
        prop_type (str): type of energy loss mechanism; can be stochastic or continuous
        stats (int): statistics or number of neutrinos injected
        pexit_table (`~astropy.table.Table`): charged lepton exit probability values
            Each table will need the following columns:
            - ``angle``: Earth emergence angle, in degrees [``angle``]
            - ``no_regen``: Exit probability value without any regeneration [``no_regen``]
            - ``regen``: Exit probability value including regeneration (max. of 6 rounds of regen) [``regen``]
        arg (str, optional): additional arguments at the end of the file name. Defaults to None.

    Returns:
        None
    """
    log_energy = np.log10(energy)
    energy_str = str(log_energy)

    out_file = output_file(nu_type,ch_lepton,idepth,cross_section_model,pn_model,prop_type,stats,arg)

    pexit_table.write(out_file, path='Exit_Probability/%s' % energy_str, append=True, overwrite=True)

    return None

def get_pexit(nu_type, ch_lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats, out=False, arg=None):
    """get charged lepton exit probability values

    Args:
        nu_type (str): type of neutrino particle; can be neutrino or anti_neutrino
        ch_lepton (str): type of charged lepton; can be tau or muon
        energy (float): ingoing neutrino (or anti-neutrino) energy, in GeV
        idepth (int): depth of water layer, in km
        cross_section_model (str): neutrino cross-section model
        pn_model (str): photonuclear energy loss model
        prop_type (str): type of energy loss mechanism; can be stochastic or continuous
        stats (int): statistics or number of neutrinos injected
        out (bool, optional): saves the data as an ascii ecsv file if set to True; returns the array
        value if set to False. Defaults to False.
        arg (str, optional): additional arguments at the end of the file name. Defaults to None.

    Returns:
        None/ndarray: None if out=True otherwise otherwise:
            2D numpy array of shape (2,x) containing charged lepton exit probability values without
            regeneration and including regeneration; x is the number of Earth emergence angles.
    """
    log_energy = np.log10(energy)
    energy_str = str(log_energy)

    in_file = output_file(nu_type,ch_lepton,idepth,cross_section_model,pn_model,prop_type,stats,arg)

    pexit_table = Table.read(in_file, path='Exit_Probability/%s' % energy_str)

    no_regen = pexit_table['no_regen']
    regen = pexit_table['regen']
    out_arr = np.asarray([no_regen, regen])

    if out:
        if nu_type=='neutrino':nu_type='nu'
        else:nu_type='anu'
        fnm = "pexit_%s_%s_%sGeV_%skm_%s_%s_%s_%s.ecsv" % (nu_type, ch_lepton, energy_str, idepth, cross_section_model, pn_model, prop_type, sci_str(stats))
        ascii.write(pexit_table, fnm, format='ecsv', fast_writer=True, overwrite=True)
        return print('Exit probability data saved to file %s' % fnm)

    return out_arr

def add_clep_out(nu_type, ch_lepton, energy, angle, idepth, cross_section_model, pn_model, prop_type, stats, clep_table, arg=None):
    """adds outgoing charged lepton energy values to output file

    Args:
        nu_type (str): type of neutrino particle; can be neutrino or anti_neutrino
        ch_lepton (str): type of charged lepton; can be tau or muon
        energy (float): ingoing neutrino (or anti-neutrino) energy, in GeV
        angle (float): earth emergence angle, in degrees
        idepth (int): depth of water layer, in km
        cross_section_model (str): neutrino cross-section model
        pn_model (str): photonuclear energy loss model
        prop_type (str): type of energy loss mechanism; can be stochastic or continuous
        stats (int): statistics or number of neutrinos injected
        clep_table (`~astropy.table.Table`): charged lepton outgoing energy values
            Each table will need the following column:
            - ``lep_energy``: Outgoing charged lepton energy, in log_10(E) GeV [``lep_energy``]
        arg (str, optional): additional arguments at the end of the file name. Defaults to None.

    Returns:
        None
    """
    log_energy = np.log10(energy)
    energy_str = str(log_energy)

    out_file = output_file(nu_type,ch_lepton,idepth,cross_section_model,pn_model,prop_type,stats,arg)

    clep_table.write(out_file, path='CLep_out_energies/%s/%s' % (energy_str,angle), append=True, overwrite=True)
    return None

def get_clep_out(nu_type, ch_lepton, energy, angle, idepth, cross_section_model, pn_model, prop_type, stats, out=False, arg=None):
    """get charged lepton outgoing energy values

    Args:
        nu_type (str): type of neutrino particle; can be neutrino or anti_neutrino
        ch_lepton (str): type of charged lepton; can be tau or muon
        energy (float): ingoing neutrino (or anti-neutrino) energy, in GeV
        angle (float): earth emergence angle, in degrees
        idepth (int): depth of water layer, in km
        cross_section_model (str): neutrino cross-section model
        pn_model (str): photonuclear energy loss model
        prop_type (str): type of energy loss mechanism; can be stochastic or continuous
        stats (int): statistics or number of neutrinos injected
        out (bool, optional): saves the data as an ascii ecsv file if set to True; returns the array
        value if set to False. Defaults to False.
        arg (str, optional): additional arguments at the end of the file name. Defaults to None.

    Returns:
        None/ndarray: None if out=True otherwise otherwise:
            1D numpy array of shape (x,) containing outgoing charged lepton energy values, in GeV;
            x is the number of outgoing charged leptons
    """
    log_energy = np.log10(energy)
    energy_str = str(log_energy)

    in_file = output_file(nu_type,ch_lepton,idepth,cross_section_model,pn_model,prop_type,stats,arg)

    e_out = Table.read(in_file, 'CLep_out_energies/%s/%s' % (energy_str,angle))
    out_lep = 10**(np.asarray(e_out['lep_energy'])) # changed 13/7/21

    if out:
        if nu_type=='neutrino':nu_type='nu'
        else:nu_type='anu'
        clep_meta = OrderedDict({'Description':'Outgoing %s energies' % ch_lepton,
                                'lep_energy':'Outgoing %s energy, in GeV'})
        clep_table = Table([out_lep], names=('lep_energy',), meta=clep_meta)
        fnm = "CLep_out_%s_%s_%sGeV_%.2fdeg_%skm_%s_%s_%s_%s.ecsv" % (nu_type, ch_lepton, energy_str, angle, idepth, cross_section_model, pn_model, prop_type, sci_str(stats))
        ascii.write(clep_table, fnm, format='ecsv', fast_writer=True, overwrite=True)
        return print('Outgoing %s energy data saved to file %s' % (ch_lepton,fnm))

    return out_lep

def add_cdf(nu_type, ch_lepton, idepth, cross_section_model, pn_model, prop_type, stats, bins=np.logspace(-5,0,51), arg=None):
    """adds outgoing charged lepton energy CDF values to output file

    Args:
        nu_type (str): type of neutrino particle; can be neutrino or anti_neutrino
        ch_lepton (str): type of charged lepton; can be tau or muon
        idepth (int): depth of water layer, in km
        cross_section_model (str): neutrino cross-section model
        pn_model (str): photonuclear energy loss model
        prop_type (str): type of energy loss mechanism; can be stochastic or continuous
        stats (int): statistics or number of neutrinos injected
        bins (ndarray, optional): bins for computing CDF values. Defaults to np.logspace(-5,0,51)
        arg (str, optional): additional arguments at the end of the file name. Defaults to None.

    Returns:
        None
    """
    out_file = output_file(nu_type,ch_lepton,idepth,cross_section_model,pn_model,prop_type,stats,arg)

    cdf_meta = OrderedDict({'Description':'Outgoing %s energy CDF' % ch_lepton,
                            'z':'z=E_tau/E_nu',
                            'x':'Outgoing %s energy cdf values for x degrees' % ch_lepton})
    with h5py.File(out_file, 'a') as hf:
        energies = sorted([float(i) for i in hf['CLep_out_energies'].keys()])

        for energy in energies: # these are log_10(GeV)
            energy_str = str(energy)
            angles = sorted([float(i) for i in hf['CLep_out_energies'][energy_str].keys()])
            cdf_angles = []
            cdf_angles.append(bins)
            for angle in angles:
                e_out = get_clep_out(nu_type, ch_lepton, 10**energy, angle, idepth, cross_section_model, pn_model, prop_type, stats, arg=arg)
                if np.array_equal(bins,np.logspace(-5,0,51)):count, bins_count = np.histogram(e_out/10**energy, bins) # because for nuSpaceSim, z=E_tau(or E_mu)/E_nu
                else:count, bins_count = np.histogram(e_out, bins) # for user defined bins
                pdf = count / sum(count)
                cdf = np.insert(np.cumsum(pdf),0,0) # pad at the beginning with 0 because np.histogram 'eats' the first index/value
                cdf_angles.append(cdf)
            cdf_table = Table(cdf_angles, names=('z',*angles), meta=cdf_meta)
            cdf_table.write(out_file, path='CLep_out_cdf/%s' % (energy_str), append=True, overwrite=True)

    print("CDF tables created!")
    return None

def get_cdf(nu_type, ch_lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats, out=False, arg=None):
    """get charged lepton outgoing energy CDF values

    Args:
        nu_type (str): type of neutrino particle; can be neutrino or anti_neutrino
        ch_lepton (str): type of charged lepton; can be tau or muon
        energy (float): ingoing neutrino (or anti-neutrino) energy, in GeV
        idepth (int): depth of water layer, in km
        cross_section_model (str): neutrino cross-section model
        pn_model (str): photonuclear energy loss model
        prop_type (str): type of energy loss mechanism; can be stochastic or continuous
        stats (int): statistics or number of neutrinos injected
        out (bool, optional): saves the data as an ascii ecsv file if set to True; returns the array
        value if set to False. Defaults to False.
        arg (str): type of neutrino for part_type=nu or the material for part_type=tau and part_type=muon

    Returns:
        None/tuple: None if out=True otherwise tuple containing:
            z_vals (ndarray): 1D array of shape (n,) containing bins from the output file -> CLep_out_cdf table
            angles (ndarray): 1D array of shape (n,) containing earth emergence angles in the output file
            cdf_arr (ndarray): 1D array of shape (n,) containing CDF values of outgoing charged lepton energies,
            binned according to the output file -> CLep_out_cdf table
    """
    log_energy = np.log10(energy)
    energy_str = str(log_energy)

    in_file = output_file(nu_type,ch_lepton,idepth,cross_section_model,pn_model,prop_type,stats,arg)

    cdf_table = Table.read(in_file, 'CLep_out_cdf/%s' % energy_str)

    cols = cdf_table.columns
    z_vals = cdf_table['z'].data
    nn = cols.pop('z')
    angles = np.asarray([float(i) for i in cols])
    cdf_arr = np.asarray([cdf_table[str(i)].data for i in angles])

    if out:
        fnm = "cdf_%s_%s_%s_%skm_%s_%s_%s_%s.ecsv" % (nu_type, ch_lepton, energy_str, idepth, cross_section_model, pn_model, prop_type, sci_str(stats))
        ascii.write(cdf_table, fnm, format='ecsv', fast_writer=True, overwrite=True)
        return print('Outgoing %s energy CDF data saved to file %s' % (ch_lepton,fnm))
    return z_vals, angles, cdf_arr

def interp_pexit(nu_type, ch_lepton, energy, angle, idepth, cross_section_model, pn_model, prop_type, stats, method='linear', arg=None):
    """interpolates exit probability value at a given energy & angle

    Args:
        nu_type (str): type of neutrino particle; can be neutrino or anti_neutrino
        ch_lepton (str): type of charged lepton; can be tau or muon
        energy (float): energy to be interpolated at, in GeV
        angle (float): earth emergence angle to be interpolated at, in degrees
        idepth (int): depth of water layer, in km
        cross_section_model (str): neutrino cross-section model
        pn_model (str): photonuclear energy loss model
        prop_type (str): type of energy loss mechanism; can be stochastic or continuous
        stats (int): statistics or number of neutrinos injected
        method (str, optional): method for 2D interpolation; can be linear, nearest or splinef2d. Defaults to linear
        arg (str, optional): additional arguments at the end of the file name. Defaults to None

    Returns:
        float: interpolated p_exit value at (energy,angle)
    """
    in_file = output_file(nu_type,ch_lepton,idepth,cross_section_model,pn_model,prop_type,stats,arg=arg)

    with h5py.File(in_file, 'r') as hf:
        energies = 10**np.asarray(sorted([float(i) for i in hf['CLep_out_energies'].keys()])) # these energies are in GeV
        angles = np.asarray(sorted([float(i) for i in hf['CLep_out_energies'][str(np.log10(energies[0]))].keys()])) # get Earth emergence angles

    p_exit = np.asarray([get_pexit(nu_type, ch_lepton, i, idepth, cross_section_model, pn_model, prop_type, stats)[1] for i in energies]).reshape(len(energies),len(angles)) # p_exit[i,j] = [energy,angle]

    points = (energies, angles)
    point = (energy, angle)
    interp_val = float(interpn(points, p_exit, point, method=method))

    return interp_val

def interp_cdf(nu_type, ch_lepton, energy, angle, idepth, cross_section_model, pn_model, prop_type, stats, z=None, arg=None):
    """interpolates CDF values at given energy, angle and z (bin) value

    Args:
        nu_type (str): type of neutrino particle; can be neutrino or anti_neutrino
        ch_lepton (str): type of charged lepton; can be tau or muon
        energy (float): energy to be interpolated at, in GeV
        angle (float): earth emergence angle to be interpolated at, in degrees
        idepth (int): depth of water layer, in km
        cross_section_model (str): neutrino cross-section model
        pn_model (str): photonuclear energy loss model
        prop_type (str): type of energy loss mechanism; can be stochastic or continuous
        stats (int): statistics or number of neutrinos injected
        z (float, optional): z value to be interpolated at; z=E_tau (or E_muon)/E_nu
        If None, it will be the array in the output file -> CLep_out_cdf
        Defaults to None
        arg (str, optional): additional arguments at the end of the file name. Defaults to None

    Returns:
        float/ndarray: interpolated cdf value (float) at z if z is provided, or
        interpolated cdf array (ndarray) of size (len(z_vals)) at z array if z is None
    """
    in_file = output_file(nu_type,ch_lepton,idepth,cross_section_model,pn_model,prop_type,stats,arg=arg)

    with h5py.File(in_file, 'r') as hf:
        energies = 10**np.asarray(sorted([float(i) for i in hf['CLep_out_energies'].keys()])) # these energies are in GeV
        # angles = np.asarray(sorted([float(i) for i in hf['CLep_out_energies'][str(np.log10(energies[0]))].keys()])) # get Earth emergence angles
        z_vals = get_cdf(nu_type, ch_lepton, energies[0], idepth, cross_section_model, pn_model, prop_type, stats)[0] # get whatever z_vals are in the output file
        angles = get_cdf(nu_type, ch_lepton, energies[0], idepth, cross_section_model, pn_model, prop_type, stats)[1] # get Earth emergence angles

    cdf_vals = np.asarray([get_cdf(nu_type, ch_lepton, i, idepth, cross_section_model, pn_model, prop_type, stats)[2] for i in energies]).reshape(len(energies),len(angles),len(z_vals)) # cdf_vals[i,j,k] = [energy,angle,cdf_val]

    points = (energies, angles, z_vals) # 3D 'coordinates' or grid

    if z is None: # if the z value to be interpolated at is not provided
        out_arr = []
        for z in z_vals:
            point = (energy, angle, z)
            out_arr.append(interpn(points, cdf_vals, point))
        interp_arr = np.asarray(out_arr).flatten()
    else: # if the z value is provided by the user
        point = (energy, angle, z)
        interp_arr = float(interpn(points, cdf_vals, point)) # not an array, just a value at z

    return interp_arr

def sort_htc_files(nu_type, ch_lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats, cdf_bins=None):
    """processes files created when the code is run with HTC mode on

    Args:
        nu_type (str): type of neutrino particle; can be neutrino or anti_neutrino
        ch_lepton (str): type of charged lepton; can be tau or muon
        energy (float): ingoing neutrino (or anti-neutrino) energy, in GeV
        idepth (int): depth of water layer, in km
        cross_section_model (str): neutrino cross-section model
        pn_model (str): photonuclear energy loss model
        prop_type (str): type of energy loss mechanism; can be stochastic or continuous
        stats (int): statistics or number of neutrinos injected

    Returns:
        None
    """

    files_path = "osg_out/" # my output folder from OSG is osg_out/energy/*

    make_array = lambda x : x if isinstance(x, Iterable) else np.array([x]) # to avoid errors with single or no columns

    log_energy = int(np.log10(energy))

    eout_files = sorted(glob.glob(files_path + str(log_energy) + "/" + "eout_*"))
    eout_files = sorted(eout_files, key=len)

    assert len(eout_files) == len(sorted(glob.glob(files_path + str(log_energy) + "/" + "pexit_*"))) # make sure len(eout_files) == len(pexit_files)

    p_angle_lst = []
    p_noregen_lst = []
    p_regen_lst = []

    for i in range(len(eout_files)):
        fnm = eout_files[i].replace(".dat","")
        angle = float(fnm.split("_")[-1])
        e_out = make_array(np.genfromtxt(eout_files[i]))
        e_out = patch_for_astropy(e_out)

        p_angle, p_noregen, p_regen = np.genfromtxt(files_path + str(log_energy) + "/" + "pexit_%.2E_%.2f.dat" % (energy,angle), usecols=(1,2,3), unpack=True)

        p_angle_lst.append(p_angle)
        p_noregen_lst.append(p_noregen)
        p_regen_lst.append(p_regen)

        clep_meta = OrderedDict({'Description':'Outgoing %s energies' % ch_lepton,
                                'lep_energy':'Outgoing %s energy, in log_10(E) GeV'})
        clep_table = Table([e_out], names=('lep_energy',), meta=clep_meta)

        add_clep_out(nu_type, ch_lepton, energy, angle, idepth, cross_section_model, pn_model, prop_type, stats, clep_table)
        print("CLep_out processed successfully")

    pexit_angle = patch_for_astropy(np.asarray(p_angle_lst))
    pexit_noregen = patch_for_astropy(np.asarray(p_noregen_lst))
    pexit_regen = patch_for_astropy(np.asarray(p_regen_lst))

    pexit_meta = OrderedDict({'Description':'Exit probability for %s' % ch_lepton,
                              'angle':'Earth emergence angle, in degrees',
                              'no_regen':'Exit probability without including any %s regeneration' % ch_lepton,
                              'regen':'Exit probability including %s regeneration' % ch_lepton})

    pexit_table = Table([pexit_angle, pexit_noregen, pexit_regen], names=('angle','no_regen','regen'), meta=pexit_meta)

    add_pexit(nu_type, ch_lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats, pexit_table) # adds p_exit results to output file
    print("P_exit processed successfully")

    if cdf_bins is None:add_cdf(nu_type, ch_lepton, idepth, cross_section_model, pn_model, prop_type, stats) # adds the binned cdf values for all neutrino energies and angles in an output file, to the output file.
    else:add_cdf(nu_type, ch_lepton, idepth, cross_section_model, pn_model, prop_type, stats, bins=cdf_bins) # adds the binned cdf values for all neutrino energies and angles in an output file, to the output file.
    return None

# =============================================================================
# Test
# =============================================================================
if __name__ == "__main__":
    nu_type = 'neutrino'
    ch_lepton = 'tau'
    energy = 1e7
    idepth = 4
    cross_section_model = 'ct18nlo'
    pn_model = 'allm'
    prop_type = 'stochastic'
    stats = 1e7
    cdf_only = 'no'
    cdf_bins = np.logspace(-5,0,51) # nuSpaceSim binning
    pass
    # sort_htc_files(nu_type, ch_lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats, cdf_only)