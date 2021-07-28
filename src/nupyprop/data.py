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
import os

# pwd = os.getcwd()

E_nu = np.logspace(3,12,91,base=10).astype(np.float64)
E_lep = np.logspace(0,12,121,base=10).astype(np.float64)

nu_models = ['allm', 'bdhm', 'ct18nlo', 'nct15']
pn_models = ['brem', 'pair', 'pn_allm', 'pn_bb']

ref = importlib_resources.files('nupyprop.datafiles') / 'lookup_tables.h5' # path for lookup_tables
custom_models = importlib_resources.files('nupyprop.models') / '*.ecsv' # path for custom model files
nu_xc_ctw = importlib_resources.files('nupyprop.models') / 'xc_neutrino_ctw.ecsv'

os.path.exists(nu_xc_ctw)

with importlib_resources.as_file(nu_xc_ctw) as custom_model:
    xc_table = Table.read(custom_model, format='ascii.ecsv')

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
    pn_model = pn_model.replace("pn_","")
    fnm = "output_%s_%s_%s_%s_%s_%s_%s.h5" % (nu_type,lepton,idepth_str,cross_section_model,pn_model,prop_type,stats_str)
    return fnm

def sci_str(exp_value):
    dec = Decimal(exp_value)
    str_val = ('{:.' + str(len(dec.normalize().as_tuple().digits) - 1) + 'e}').format(dec).replace('+', '')
    return str_val

# def chk_file(nu_type, lepton, energy, angle, idepth, cross_section_model, pn_model, prop_type, stats):
#     '''

#     Parameters
#     ----------
#     nu_type : str
#         Type of neutrino particle. Can be nu (neutrino) or anu (anti-neutrino).
#     lepton : str
#         Type of lepton. Can be tau or muon.
#     energy : float
#         Neutrino energy, in GeV.
#     angle : int
#         Earth emergence angle (beta), in degrees.
#     idepth : int
#         Depth of water layer in km.
#     cross_section_model : str
#         Neutrino cross-section model.
#     pn_model : str
#         Photonuclear energy loss model.
#     prop_type : str
#         Type of energy loss mechanism. Can be stochastic or continuous.
#     stats : float
#         Statistics or number of neutrinos injected.

#     Returns
#     -------
#     int
#         0 is to stop execution (save the old ouput file); 1 is for replacing the old output file; 2 is to append to/overwrite the old file (only do this if you know what you're doing or else you'll end up with mixed results!).
#         Option no. 2 can be used if you need to 'add' more results for the same set of parameters or in case of abrupt code termination.

#     '''
#     fnm = output_file(nu_type,lepton,idepth,cross_section_model,pn_model,prop_type,stats)
#     if os.path.exists(fnm):
#         # try:

#         choice = input('There already exists an output file with these set of parameters (%s). Press \'d\' for deleting the output old file and creating a new output file, \'s\' for keeping the old output file or \'o\' for overwriting the old output file: ' % fnm)

#         if choice not in {"d", "s", "o"}:
#             print("Invalid option. Please enter \'d\', \'s\' or \'o\'")
#             return chk_file(nu_type, lepton, idepth, cross_section_model, pn_model, prop_type, stats)
#         elif choice == 's':
#             return 0
#         elif choice == 'd':
#             os.remove(fnm)
#             return 1

#     else: # so basically choice = 'o'
#         return 2 # output file non existant or overwrite enabled
#     # return out_val


def add_trajs(type_traj, idepth, traj_table):
    '''

    Parameters
    ----------
    type_traj : str
        Type of trajectory. Can be col (for column depth) or water (for water depth).
    idepth : int
        Depth of water layer in km.
    traj_table : `~astropy.table.Table`
        Table containing the trajectories. E
    data_table : `~astropy.table.Table` or list of `~astropy.table.Table`
        Table containing the observed spectrum. If multiple tables are passed
        as a string, they will be concatenated in the order given. Each table
        needs at least these columns, with the appropriate associated units
        (with the physical type indicated in brackets below) as either a
        `~astropy.units.Unit` instance or parseable string:

        - ``energy``: Observed photon energy [``energy``]
        - ``flux``: Observed fluxes [``flux`` or ``differential flux``]
        - ``flux_error``: 68% CL gaussian uncertainty of the flux [``flux`` or
          ``differential flux``]. It can also be provided as ``flux_error_lo``
          and ``flux_error_hi`` (see below).

    Returns
    -------
    None
        Adds trajectory lookup tables to lookup_tables.h5.

    '''
    if type_traj == 'col':branch = 'Column_Trajectories' # sub-sub branch inside the Earth/traj_idepth branch
    elif type_traj == 'water':branch = 'Water_Trajectories'

    with importlib_resources.as_file(ref) as lookup_tables:
        traj_table.write(lookup_tables, path='Earth/%s/%skm' % (branch,str(idepth)), append=True, overwrite=True)
    return print("%s lookup table successfully created for idepth = %s" % (branch,str(idepth)))

def get_trajs(type_traj, angle, idepth, out=False):
    '''

    Parameters
    ----------
    type_traj : str
        Type of trajectory. Can be col (for column depth) or water (for water depth).
    angle : float
        Earth emergence angle in degrees.
    idepth : int
        Depth of water layer in km.
    out : boolean, optional
        Set this to True to write output to file. The default is False.

    Returns
    -------
    tuple
        tuple of 1D arrays - (xalong,cdalong) [for col] and (chord,water) [for water].
        xalong - distance in water, in km.
        cdalong - column depth at xalong, in g/cm^2.
        chord - chord length, in km.
        water - final water layer distance, in km.

    '''
    if type_traj == 'col':
        with importlib_resources.as_file(ref) as lookup_tables:
            traj_table = Table.read(lookup_tables,path='Earth/Column_Trajectories/%skm' % str(idepth))

        sliced_table = traj_table[traj_table['beta']==angle]
        xalong = np.asfortranarray(sliced_table['xalong'].T)
        cdalong = np.asfortranarray(sliced_table['cdalong'].T)

        if out:
            fnm = "%s_%sdeg_%skm.ecsv" % (type_traj,angle,idepth)
            ascii.write(traj_table, fnm, format='ecsv', fast_writer=True, overwrite=True)
            return print('Column trajectory data saved to file %s' % fnm)
        return xalong, cdalong

    elif type_traj == 'water':
        with importlib_resources.as_file(ref) as lookup_tables:
            traj_table = Table.read(lookup_tables,path='Earth/Water_Trajectories/%skm' % str(idepth))

        chord = float(traj_table['chord'][traj_table['beta']==angle])
        water = float(traj_table['water'][traj_table['beta']==angle])

        if out:
            fnm = "%s_%sdeg_%skm.ecsv" % (type_traj,angle,idepth)
            ascii.write(traj_table, fnm, format='ecsv', fast_writer=True, overwrite=True)
            return print('Water trajectory data saved to file %s' % fnm)
        return chord, water
    return "Error in get_trajs in Data"

def add_xc(part_type, xc_table, **kwargs):
    '''

    Parameters
    ----------
    part_type : str
        Neutrino or lepton? Can be nu or tau or muon.
    xc_table : dict or ndarray
        1D containing neutrino cross-section values or 1D array containing lepton cross-section values, in cm^2.
    model : str
        Neutrino cross section model.
    **kwargs
        nu_type: Type of neutrino particle. Can be neutrino or anti-neutrino.
        material: Material of propagation, for leptons.

    Returns
    -------
    None
        Creates neutrino/lepton cross-section lookup entries in lookup_tables.h5.

    '''
    if part_type=='nu':
        nu_type = kwargs['nu_type']
        with importlib_resources.as_file(ref) as lookup_tables:
            xc_table.write(lookup_tables, path='Neutrinos/%s/xc' % nu_type, append=True, overwrite=True)
        return print("%s_sigma CC & NC lookup tables successfully created" % part_type)
    else: # lepton energy loss XC
        material = kwargs['material']
        with importlib_resources.as_file(ref) as lookup_tables:
            xc_table.write(lookup_tables, path='Charged_Leptons/%s/%s/xc' % (part_type,material), append=True, overwrite=True)
        return print("%s_sigma lookup table successfully created in %s" % (part_type, material))
    return None

def get_xc(part_type, model, out=False, **kwargs):
    '''

    Parameters
    ----------
    part_type : str
        Neutrino or lepton? Can be nu or tau or muon.
    model : str
        Neutrino cross-section/lepton photonuclear energy loss model.
    out : boolean, optional
        Set this to True to write output to file. The default is False.
    **kwargs
        nu_type: Type of neutrino particle. Can be neutrino or anti-neutrino, for neutrinos.
        material: Material of propagation, for leptons.

    Returns
    -------
    ndarray - 2D (neutrino) or 1D (lepton) cross-section array.
    if part_type = nu; cscc/csnc = neutrino-nucleon charged/neutral current cross sections, in cm^2.
    if part_type = tau/muon; out_arr = lepton-nucleon cross section, in cm^2.

    '''
    if part_type=='nu':
        nu_type = kwargs['nu_type']
        if nu_type=='anti-neutrino':nu_type='anti_neutrino'

        if model in nu_models:
            with importlib_resources.as_file(ref) as lookup_tables:
                xc_table = Table.read(lookup_tables,path='Neutrinos/%s/xc' % nu_type)
        else:
            fnm = "xc_%s_%s.ecsv" % (nu_type,model)
            if not os.path.exists(fnm):
                return print("Error! %s file does not exist" % fnm)

            with importlib_resources.as_file(fnm) as lookup_table:
                xc_table = Table.read(lookup_table, format='ascii.ecsv')

        cscc = xc_table['sigma_cc_%s' % model]
        csnc = xc_table['sigma_nc_%s' % model]
        out_arr = np.asarray([cscc,csnc])

        if out:
            fnm = "xc_%s_%s.ecsv" % (nu_type,model)
            out_table = Table([xc_table['energy'], cscc, csnc], meta=xc_table.meta)
            ascii.write(out_table, fnm, format='ecsv', fast_writer=False, overwrite=True)
            return print('%s cross-section data saved to file %s' % (nu_type,fnm))

            return np.asfortranarray(out_arr.T)

    else: # energy loss; part_type == 'tau' or 'muon'
        material = kwargs['material']
        # if
        try:
            with importlib_resources.as_file(ref) as lookup_tables:
                xc_table = Table.read(lookup_tables,path='Charged_Leptons/%s/%s/xc' % (part_type,material))
            out_arr = np.asarray(xc_table['sigma_%s' % model])

            if out:
                fnm = "xc_%s_%s_%s.ecsv" % (part_type,model,material)
                out_table = Table([xc_table['energy'], out_arr], meta=xc_table.meta)
                ascii.write(out_table, fnm, format='ecsv', fast_writer=False, overwrite=True)
                return print('%s cross-section data saved to file %s' % (part_type,fnm))

            return np.asfortranarray(out_arr.T)

        except KeyError:
            model = str(input(("Error finding cross-section values for %s model, please enter a valid model name: " % model)))
            return None
    return None

def add_ixc(part_type, ixc_table, **kwargs):
    '''

    Parameters
    ----------
    part_type : str
        Neutrino or lepton? Can be nu or tau or muon.
    ixc_dict : dict
        Integrated cross-section CDF value dictionary for neutrinos/anti-neutrinos or leptons.
    model : str
        Neutrino cross-section/lepton photonuclear energy loss model.
    **kwargs
        nu_type: Type of neutrino particle. Can be neutrino or anti-neutrino, for neutrinos.
        material: Material of propagation, for leptons.

    Returns
    -------
    None
        Creates neutrino/lepton integrated cross-section lookup entries in lookup_tables.h5.

    '''
    if part_type == 'nu':
        nu_type = kwargs['nu_type']

        with importlib_resources.as_file(ref) as lookup_tables:
            ixc_table.write(lookup_tables, path='Neutrinos/%s/ixc' % nu_type, append=True, overwrite=True)

        return print("%s_sigma CDF CC & NC lookup tables successfully created for" % nu_type)

    else: # energy_loss; part_type == 'muon' or 'tau'
        material = kwargs['material']

        with importlib_resources.as_file(ref) as lookup_tables:
            ixc_table.write(lookup_tables, path='Charged_Leptons/%s/%s/ixc' % (part_type,material), append=True, overwrite=True)

        return print("%s_sigma CDF lookup table successfully created in %s" % (part_type, material))
    return None

def get_ixc(part_type, model, out=False, **kwargs):
    '''

    Parameters
    ----------
    part_type : str
        Neutrino or lepton? Can be nu or tau or muon.
    model : str
        Neutrino cross-section/lepton photonuclear energy loss model.
    out : boolean, optional
        Set this to True to write output to file. The default is False.
    **kwargs
        nu_type: Type of neutrino particle. Can be neutrino or anti-neutrino.
        material: Material of propagation, for leptons.

    Returns
    -------
    ndarray - 2D (neutrino) or 1D (lepton) integrated cross-section CDF array.
    if part_type = nu; cscc/csnc = neutrino-nucleon charged/neutral current integrated cross section CDF values.
    if part_type = tau/muon; out_arr = lepton-nucleon integrated cross section CDF values.

    '''
    # v2 = -np.linspace(0.1,3,num=30)
    # yvals = 10**v2 # The integrated cross-section values should go from y = 0, 10^(-0.1),..., 10^(-3). This is a convention we chose to adopt.
    # yvals = np.insert(yvals,0,0)

    if part_type == 'nu':
        nu_type = kwargs['nu_type']
        if nu_type=='anti-neutrino':nu_type='anti_neutrino'
        try:
            with importlib_resources.as_file(ref) as lookup_tables:
                ixc_table = Table.read(lookup_tables,path='Neutrinos/%s/ixc' % nu_type)

            ixc_cc = np.asarray([ixc_table['cc_cdf_%s' % model][ixc_table['energy']==i] for i in E_nu])
            ixc_nc = np.asarray([ixc_table['nc_cdf_%s' % model][ixc_table['energy']==i] for i in E_nu])

            out_arr = np.asarray([ixc_cc, ixc_nc])

            if out:
                fnm = "ixc_%s_%s.ecsv" % (nu_type,model)
                out_table = Table([ixc_table['energy'], ixc_table['y'], ixc_table['cc_cdf_%s' % model], ixc_table['nc_cdf_%s' % model]], meta=ixc_table.meta)
                ascii.write(out_table, fnm, format='ecsv', fast_writer=False, overwrite=True)
                return print('%s cross-section CDF data saved to file %s' % (nu_type,fnm))

            return np.asfortranarray(out_arr.T)


        except KeyError or TypeError:
            model = str(input(("Error finding integrated cross-section values for %s model, please enter a valid model name." % str(model))))
            return None
    else: # energy loss; ixc_type == 'tau' or 'muon'
        try:
            material = kwargs['material']

            with importlib_resources.as_file(ref) as lookup_tables:
                ixc_table = Table.read(lookup_tables,path='Charged_Leptons/%s/%s/ixc' % (part_type,material))

            out_arr = np.asarray([ixc_table['cdf_%s' % model][ixc_table['energy']==i] for i in E_lep])

            if out:
                fnm = "ixc_%s_%s_%s.ecsv" % (part_type,model,material)
                out_table = Table([ixc_table['energy'], ixc_table['y'], ixc_table['cdf_%s' % model]], meta=ixc_table.meta)
                ascii.write(out_table, fnm, format='ecsv', fast_writer=True, overwrite=True)
                return print('%s cross-section CDF data saved to file %s' % (part_type,fnm))

            return np.asfortranarray(out_arr)

        except KeyError or TypeError:
            model = str(input("Error finding energy loss cross-section values for %s, please enter a valid model name: " % str(model)))
            return None
    return None

def add_alpha(lepton, material, alpha_table):
    '''

    Parameters
    ----------
    alpha : dict
        Ionization energy loss dictionary.
    lepton : str
        Lepton. Can be tau or muon.
    material : str
        Material of propagation.

    Returns
    -------
    None
        Creates lepton ionization energy loss lookup entries in lookup_tables.h5.

    '''

    with importlib_resources.as_file(ref) as lookup_tables:
        alpha_table.write(lookup_tables, path='Charged_Leptons/%s/%s/alpha' % (lepton,material), append=True, overwrite=True)
    return print("%s_alpha lookup table successfully created for %s" % (lepton,material))

def get_alpha(lepton, material, out=False):
    '''

    Parameters
    ----------
    lepton : str
        Lepton. Can be tau or muon.
    material : str
        Material of propagation.
    out : boolean, optional
        Set this to True to write output to file. The default is False.

    Returns
    -------
    ndarray
        1D array of lepton ionization energy loss, in (GeV*cm^2)/g.

    '''
    with importlib_resources.as_file(ref) as lookup_tables:
        alpha_table = Table.read(lookup_tables,path='Charged_Leptons/%s/%s/alpha' % (lepton,material))
    out_arr = alpha_table['alpha']

    if out:
        fnm = "alpha_%s_%s.ecsv" % (lepton,material)
        ascii.write(alpha_table, fnm, format='ecsv', fast_writer=True, overwrite=True)
        return print('Alpha data saved to file %s' % fnm)

    return np.asfortranarray(out_arr.T)

def add_beta(lepton, material, beta_table):
    '''

    Parameters
    ----------
    beta_arr : ndarray
        1D array containing energy loss parameter (beta), in cm^2/g.
    lepton : str
        Lepton. Can be tau or muon.
    material : str
        Material of propagation.
    # model : str
    #     Lepton energy loss model/process.

    Returns
    -------
    None
        Creates lepton (non-ionization) energy loss lookup entries in lookup_tables.h5.

    '''
    with importlib_resources.as_file(ref) as lookup_tables:
        beta_table.write(lookup_tables, path='Charged_Leptons/%s/%s/beta' % (lepton,material), append=True, overwrite=True)
    return print("%s_beta in %s lookup table successfully created" % (lepton,material))

def get_beta(lepton, material, model, beta_type, out=False):
    '''

    Parameters
    ----------
    lepton : str
        Lepton. Can be tau or muon.
    material : str
        Material of propagation.
    model : str
        Lepton energy loss model/process.
    beta_type : str
        Can be cut (for stochastic energy loss) or full (for continuous energy loss).
    out : boolean, optional
        Set this to True to write output to file. The default is False.

    Returns
    -------
    ndarray
        1D array containing lepton energy loss model/process beta value, in cm^2/g.

    '''
    with importlib_resources.as_file(ref) as lookup_tables:
        beta_table = Table.read(lookup_tables,path='Charged_Leptons/%s/%s/beta' % (lepton,material))
    out_arr = beta_table['beta_%s_%s' % (model,beta_type)]

    if out:
        fnm = "beta_%s_%s_%s_%s.ecsv" % (lepton,model,material,beta_type)
        out_table = Table([beta_table['energy'], out_arr], meta=beta_table.meta)
        ascii.write(out_table, fnm, format='ecsv', fast_writer=True, overwrite=True)
        return print('Beta data saved to file %s' % fnm)

    return np.asfortranarray(out_arr.T)

def combine_lep(data_type, lepton, material, pn_model, **kwargs):
    '''

    Parameters
    ----------
    data_type : str
        Can be xc (lepton cross-section type), beta (lepton energy loss) or ixc (lepton integrated cross-section CDFs).
    lepton : str
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
        xc_brem = get_xc(lepton, 'brem', material=material)
        xc_pair = get_xc(lepton, 'pair', material=material)
        xc_pn = get_xc(lepton, pn_model, material=material)
        xc_arr = np.asarray([xc_brem, xc_pair, xc_pn])
        xc = np.asfortranarray(xc_arr.T)
        return xc

    elif data_type == 'beta':
        beta_type = kwargs['beta_type']
        beta_brem = get_beta(lepton, material, 'brem', beta_type)
        beta_pair = get_beta(lepton, material, 'pair', beta_type)
        beta_pn = get_beta(lepton, material, pn_model, beta_type)
        beta_arr = np.asarray([beta_brem, beta_pair, beta_pn])
        beta = np.asfortranarray(beta_arr.T)
        return beta

    elif data_type == 'ixc':
        ixc_brem = get_ixc(lepton, 'brem', material=material)
        ixc_pair = get_ixc(lepton, 'pair', material=material)
        ixc_pn = get_ixc(lepton, pn_model, material=material)
        ixc_arr = np.asarray([ixc_brem, ixc_pair, ixc_pn])
        ixc = np.asfortranarray(ixc_arr.T)
        return ixc

    return None

def add_pexit(nu_type, lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats, pexit_table):
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

    pexit_table.write(output_file(nu_type,lepton,idepth,cross_section_model,pn_model,prop_type,stats), path='Exit_Probability/%s' % energy_str, append=True, overwrite=True)

    return None

def get_pexit(nu_type, lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats, out=False):
    '''

    Parameters
    ----------
    nu_type : str
        Type of neutrino particle. Can be neutrino or anti-neutrino.
    lepton : str
        Type of lepton. Can be tau or muon.
    energy : float
        Neutrino energy, in GeV.
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
    out : boolean, optional
        Set this to True to write output to file. The default is False.

    Returns
    -------
    ndarray
        2D array containing no_regen and with_regen exit probabilities.

    '''
    log_energy = np.log10(energy)
    energy_str = str(log_energy)

    pexit_table = Table.read(output_file(nu_type,lepton,idepth,cross_section_model,pn_model,prop_type,stats),path='Exit_Probability/%s' % energy_str)

    no_regen = pexit_table['no_regen']
    regen = pexit_table['regen']
    out_arr = np.asarray([no_regen, regen])

    if out:
        fnm = "pexit_%s_%s_%s_%skm_%s_%s_%s_%s.ecsv" % (nu_type, lepton, energy_str, idepth, cross_section_model, pn_model, prop_type, sci_str(stats))
        ascii.write(pexit_table, fnm, format='ecsv', fast_writer=True, overwrite=True)
        return print('Exit probability data saved to file %s' % fnm)

    return out_arr

def add_lep_out(nu_type, lepton, energy, angle, idepth, cross_section_model, pn_model, prop_type, stats, lep_table):
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

    lep_table.write(output_file(nu_type,lepton,idepth,cross_section_model,pn_model,prop_type,stats), path='Lep_out_energies/%s/%d' % (energy_str,angle), append=True, overwrite=True)
    return None

def get_lep_out(nu_type, lepton, energy, angle, idepth, cross_section_model, pn_model, prop_type, stats, out=False):
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
    out : boolean, optional
        Set this to True to write output to file. The default is False.

    Returns
    -------
    out_lep : ndarray
        1D array containing lepton out energies, in GeV.

    '''
    log_energy = np.log10(energy)
    energy_str = str(log_energy)

    e_out = Table.read(output_file(nu_type,lepton,idepth,cross_section_model,pn_model,prop_type,stats),'Lep_out_energies/%s/%s' % (energy_str,angle))
    out_lep = 10**(np.asarray(e_out['lep_energy'])) # changed 13/7/21

    if out:
        fnm = "lep_out_%s_%s_%s_%sdeg_%skm_%s_%s_%s_%s.ecsv" % (nu_type, lepton, energy_str, angle, idepth, cross_section_model, pn_model, prop_type, sci_str(stats))
        ascii.write(10**e_out['lep_energy'], fnm, format='ecsv', fast_writer=True, overwrite=True)
        return print('Lepton outgoing energy data saved to file %s' % fnm)

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

def add_cdf(nu_type, lepton, energy, angle, idepth, cross_section_model, pn_model, prop_type, stats, lep_table):
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

    bins = np.logspace(-5,0,51) # Default binning for use with nuSpaceSim. Change if different binning required.
    lep_out = 10**lep_table['lep_energy'] # because lep_out energies are in log10(GeV)
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
    # z_cdf_df = pd.DataFrame({'z':bins, 'CDF':np.around(z_cdf,decimals=8)})

    cdf_meta = OrderedDict({'Description':'Outgoing %s energy CDF' % lepton,
                            'z':'z=E_tau/E_nu',
                            'cdf':'Outgoing %s energy cdf value' % lepton})

    cdf_table = Table([bins,z_cdf], names=('z','cdf'), meta=cdf_meta)

    cdf_table.write(output_file(nu_type,lepton,idepth,cross_section_model,pn_model,prop_type,stats), path='Lep_out_cdf/%s/%d' % (energy_str,angle), append=True, overwrite=True)

    return None

def get_cdf(nu_type, lepton, energy, angle, idepth, cross_section_model, pn_model, prop_type, stats, out=False):
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
    out : boolean, optional
        Set this to True to write output to file. The default is False.

    Returns
    -------
    cdf : dict
        Dictionary of lepton outgoing value CDFs binned in z=E_nu/E_lep.

    '''

    log_energy = np.log10(energy)
    energy_str = str(log_energy)

    cdf_table = Table.read(output_file(nu_type,lepton,idepth,cross_section_model,pn_model,prop_type,stats),'Lep_out_cdf/%s/%s' % (energy_str,angle))

    cdf = cdf_table['cdf']

    if out:
        fnm = "cdf_%s_%s_%s_%sdeg_%skm_%s_%s_%s_%s.ecsv" % (nu_type, lepton, energy_str, angle, idepth, cross_section_model, pn_model, prop_type, sci_str(stats))
        # np.savetxt(fnm, np.transpose([df.z, df.cdf]), header="z=E_lep/E_nu" + "\t" + "cdf")
        ascii.write(cdf_table, fnm, format='ecsv', fast_writer=True, overwrite=True)
        return print('Lepton outgoing energy CDF data saved to file %s' % fnm)
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
# part_type, xc_obj, model, **kwargs

# def get_add_xc(part_type, model, **kwargs):
#     if part_type == 'nu':
#         nu_type = kwargs['nu_type']
#         xc_arr = get_xc(part_type, model, nu_type=nu_type).T
#         xc_cc = xc_arr[0]
#         xc_nc = xc_arr[1]
#         xc_meta = OrderedDict({'Description':'%s-nucleon cross-section values for %s' % (nu_type,model),
#                                'energy':'Neutrino energy, in GeV',
#                                'sigma_cc':'Charged current cross-section, in cm^2',
#                                'sigma_nc':'Neutral current cross-section, in cm^2'})
#         xc_table = Table([E_nu, xc_cc, xc_nc], names=('energy','sigma_cc','sigma_nc'), meta=xc_meta)
#         add_xc(part_type,xc_table,model,nu_type=nu_type)
#     else:
#         material = kwargs['material']
#         xc_arr = get_xc(part_type, model, material=material).T
#         xc_meta = OrderedDict({'Description':'%s-nucleon cross-section values for %s in %s' % (part_type,model,material),
#                                'energy':'%s energy, in GeV' % part_type,
#                                'sigma_%s' % model:'cross-section for %s * N_A/A, in cm^2/g' % model})
#         xc_table = Table([E_lep, xc_arr], names=('energy','sigma_%s' % model), meta=xc_meta)
#         add_xc(part_type,xc_table,model,material=material)

#     return None

# =============================================================================
# Test
# =============================================================================
if __name__ == "__main__":
    # arr = get_beta('tau', 'rock', 'total', 'allm', True)
    # add_pexit_manual(1e9, np.arange(1,36), 1e8)
    pass
    # for nu_type in nu_types:
    #     for model in nu_xc:
    #         nu_ixc = get_ixc('nu', model, out=True, nu_type=nu_type)
    #         nu_xc = get_xc('nu', model, out=True, nu_type=nu_type)

    # for lepton in lep:
    #     for material in materials:
    #         get_alpha(lepton, material, out=True)
    #         for process in lep_processes:
    #             lep_xc = get_ixc(lepton, process, out=True, material=material)
    #             lep_ixc = get_ixc(lepton, process, out=True, material=material)
    #             for beta_type in beta_types:
    #                 lep_beta = get_beta(lepton, material, process, beta_type, out=True)