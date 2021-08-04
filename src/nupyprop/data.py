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
import glob
from collections.abc import Iterable

# pwd = os.getcwd()

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

    def __init__(self, fnm, message="This is either an incorrectly formatted model file or model file not found"):
        self.fnm = fnm
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.fnm} -> {self.message}'

def get_custom_path(data_type, part_type, model, *args): # get custom model file posix path
    if part_type == 'nu':
        nu_type = args[0]
        fnm = data_type + '_%s_%s.ecsv' % (nu_type,model)
        file = importlib_resources.files('nupyprop.models') / fnm
        if not os.path.exists(file):
            raise ModelError(fnm)
        return file

    elif part_type == 'tau' or part_type == 'muon':
        material = args[0]
        fnm = data_type + '_%s_pn_%s_%s.ecsv' % (part_type,model,material)
        file = importlib_resources.files('nupyprop.models') / fnm
        if not os.path.exists(file):
            raise ModelError(fnm)
        return file

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

def add_xc(part_type, xc_table, *arg):
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
        nu_type = arg[0]
        with importlib_resources.as_file(ref) as lookup_tables:
            xc_table.write(lookup_tables, path='Neutrinos/%s/xc' % nu_type, append=True, overwrite=True)
        return print("%s_sigma CC & NC lookup tables successfully created" % part_type)
    else: # lepton energy loss XC
        material = arg[0]
        with importlib_resources.as_file(ref) as lookup_tables:
            xc_table.write(lookup_tables, path='Charged_Leptons/%s/%s/xc' % (part_type,material), append=True, overwrite=True)
        return print("%s_sigma lookup table successfully created in %s" % (part_type, material))
    return None

def get_xc(part_type, model, *arg, out=False):
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
        nu_type = arg[0]
        if nu_type=='anti-neutrino':nu_type='anti_neutrino'

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
        material = arg[0]
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

def add_ixc(part_type, ixc_table, *arg):
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
        nu_type = arg[0]
        if nu_type=='anti-neutrino':nu_type='anti_neutrino'

        with importlib_resources.as_file(ref) as lookup_tables:
            ixc_table.write(lookup_tables, path='Neutrinos/%s/ixc' % nu_type, append=True, overwrite=True)

        return print("%s_sigma CDF CC & NC lookup tables successfully created for" % nu_type)

    else: # energy_loss; part_type == 'muon' or 'tau'
        material = arg[0]

        with importlib_resources.as_file(ref) as lookup_tables:
            ixc_table.write(lookup_tables, path='Charged_Leptons/%s/%s/ixc' % (part_type,material), append=True, overwrite=True)

        return print("%s_sigma CDF lookup table successfully created in %s" % (part_type, material))
    return None

def get_ixc(part_type, model, *arg, out=False):
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
        nu_type = arg[0]
        if nu_type=='anti-neutrino':nu_type='anti_neutrino'

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
        material = arg[0]

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
    alpha_arr = alpha_table['alpha']

    if out:
        fnm = "alpha_%s_%s.ecsv" % (lepton,material)
        ascii.write(alpha_table, fnm, format='ecsv', fast_writer=True, overwrite=True)
        return print('Alpha data saved to file %s' % fnm)

    return np.asfortranarray(alpha_arr.T)

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

def get_beta(lepton, model, material, *arg, out=False):
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
    beta_type = arg[0]

    if beta_type=='cut':beta_type_str = 'y_max=1e-3'
    elif beta_type=='total':beta_type_str = 'y_max'

    with importlib_resources.as_file(ref) as lookup_tables:
        beta_table = Table.read(lookup_tables,path='Charged_Leptons/%s/%s/beta' % (lepton,material))

    beta_brem = beta_table['beta_brem_%s' % beta_type] # brem is not a custom model
    beta_pair = beta_table['beta_pair_%s' % beta_type] # pair is not a custom model

    if model in pn_models: # default PN model selection
        beta_pn = beta_table['beta_pn_%s_%s' % (model,beta_type)]

    else: # custom PN model
        with importlib_resources.as_file(get_custom_path('beta',lepton,model,material)) as lookup_table:
            pn_table = Table.read(lookup_table,format='ascii.ecsv')
        beta_pn = pn_table['beta_pn_%s_%s' % (model,beta_type)]

    beta_arr = np.asarray([beta_brem, beta_pair, beta_pn])

    if out:
        fnm = "beta_%s_pn_%s_%s_%s.ecsv" % (lepton,model,material,beta_type)
        beta_meta = OrderedDict({'Description':'Model/parameterization dependent photonuclear energy loss lookup table for %s in %s' % (lepton.capitalize(),material),
                                  'energy':'%s energy, in GeV' % lepton.capitalize(),
                                  'beta_pn_%s_%s' % (model,beta_type):'Photonuclear %s energy loss model beta values integrated from y_min to %s, in cm^2/g' % (str.upper(model),beta_type_str)})
        out_table = Table([beta_table['energy'], beta_pn], names=('energy','beta_pn_%s_%s' % (model,beta_type)), meta=beta_meta)
        ascii.write(out_table, fnm, format='ecsv', fast_writer=False, overwrite=True)
        return print('Beta data saved to file %s' % fnm)

    return np.asfortranarray(beta_arr.T)


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
        ascii.write(cdf_table, fnm, format='ecsv', fast_writer=True, overwrite=True)
        return print('Lepton outgoing energy CDF data saved to file %s' % fnm)
    return cdf

def sort_htc_files(nu_type, lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats, cdf_only='no'):

    files_path = "osg_out/" # my output folder from OSG is osg_out/energy/*

    make_array = lambda x : x if isinstance(x, Iterable) else np.array([x]) # to avoid errors with single or no columns

    log_energy = int(np.log10(energy))

    eout_files = sorted(glob.glob(files_path + str(log_energy) + "/" + "eout_*"))

    assert len(eout_files) == len(sorted(glob.glob(files_path + str(log_energy) + "/" + "pexit_*")))

    for i in range(len(eout_files)): # make sure len(eout_files) ==
        fnm = eout_files[i].replace(".dat","")
        angle = float(fnm.split("_")[-1])
        e_out = make_array(np.genfromtxt(eout_files[i]))
        p_angle, p_noregen, p_regen = np.genfromtxt(files_path + str(log_energy) + "/" + "pexit_%.2E_%.2f.dat" % (energy,angle), usecols=(1,2,3), unpack=True)

        lep_meta = OrderedDict({'Description':'Outgoing %s energies' % lepton,
                                'lep_energy':'Outgoing %s energy, in log_10(E) GeV'})
        lep_table = Table([e_out], names=('lep_energy',), meta=lep_meta)
        add_cdf(nu_type, lepton, energy, angle, idepth, cross_section_model, pn_model, prop_type, stats, lep_table) # adds the binned cdf values for each energy and angle to output file

        if cdf_only == 'no': # adds lep_out energies to output file
            add_lep_out(nu_type, lepton, energy, angle, idepth, cross_section_model, pn_model, prop_type, stats, lep_table)

        pexit_meta = OrderedDict({'Description':'Exit probability for %s' % lepton,
                                  'angle':'Earth emergence angle, in degrees',
                                  'no_regen':'Exit probability without including any %s regeneration' % lepton,
                                  'regen':'Exit probability including %s regeneration' % lepton})

        pexit_table = Table([make_array(p_angle), make_array(p_noregen), make_array(p_regen)], names=('angle','no_regen','regen'), meta=pexit_meta)

        add_pexit(nu_type, lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats, pexit_table) # adds p_exit results to output file
    return None

# =============================================================================
# Test
# =============================================================================
if __name__ == "__main__":
    nu_type = 'neutrino'
    lepton = 'tau'
    energy = 1e7
    idepth = 4
    cross_section_model = 'ct18nlo'
    pn_model = 'allm'
    prop_type = 'stochastic'
    stats = 1e8
    cdf_only = 'yes'
    pass
    # sort_htc_files(nu_type, lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats, cdf_only)