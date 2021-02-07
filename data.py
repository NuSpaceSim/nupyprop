#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 16:44:13 2020

@author: sam
"""

import numpy as np
import pandas as pd
from pandas import HDFStore
# import collections
import time
from scipy import interpolate

E_nu = np.logspace(3,12,91,base=10).astype(np.float64)
E_lep = np.logspace(0,12,121,base=10).astype(np.float64)
# v2 = -np.linspace(0.1,3,num=30).astype(np.float64)
# v2 = np.insert(v2,0,0) # padding 0 at the beginning to match index with ixc entry nos. as those start with 1


def add_trajs(type_traj, idepth, traj_array):
    hdf = HDFStore('lookup_tables.h5','a')
    if type_traj == 'col':branch = 'Column_Trajectories' # sub-sub branch inside the Earth/traj_idepth branch
    elif type_traj == 'water':branch = 'Water_Trajectories'
    hdf.put('Earth/traj_%s/%s' % (str(idepth),branch), traj_array)
    hdf.close()
    return print("%s lookup table successfully created for idepth = %s" % (branch,str(idepth)))

def get_trajs(type_traj, beta, idepth): # returns {xalong:cdalong} for beta if type=col or returns chord, water for beta if type=water
    if type_traj == 'col':
        dataset = pd.read_hdf('lookup_tables.h5','Earth/traj_%s/Column_Trajectories' % str(idepth))
        dataset_sliced = dataset[dataset['beta']==beta]
        xalong = np.array(dataset_sliced.xalong)
        cdalong = np.array(dataset_sliced.cdalong)
        return xalong, cdalong

    elif type_traj == 'water':
        dataset = pd.read_hdf('lookup_tables.h5','Earth/traj_%s/Water_Trajectories' % str(idepth))
        chord = float(dataset.chord[dataset['beta']==beta])
        water = float(dataset.water[dataset['beta']==beta])
        return chord, water
    return "Error in get_trajs in Data"

def add_xc(xc_type, xc_dict, **kwargs):
    '''
    Parameters
    ----------
    model : string
        Name of the model you want to add the cross-section values of, to the hdf file
    xc_dict : dictionary
        If ixc_type='nu', this dictionary should contain 4 keys: 'nucc', 'nunc', 'anucc' & 'anunc' and their corresponding keys , and if ixc_type='tau' or 'muon', this dictionary should contain 4 keys: 'brem', 'pair', 'pn_bb' & 'pn' and their corresponding keys

    Returns
    -------
    None.

    '''
    hdf = HDFStore('lookup_tables.h5','a')
    if xc_type=='nu':
        model = kwargs['model']
        particle_type = ['nu','anu']
        for particle in particle_type:
            cc = xc_dict[particle]['cc']
            nc = xc_dict[particle]['nc']
            dframe = pd.DataFrame({'energy':E_nu, 'sigma_cc':cc, 'sigma_nc':nc})
            if 'a' in particle: # anti-neutrino group
                hdf.put('Neutrino_Cross_Sections/anti_neutrino/xc/%s' % model, dframe)
            else: # neutrino group
                hdf.put('Neutrino_Cross_Sections/neutrino/xc/%s' % model, dframe)
        hdf.close()
        return print("%s_sigma CC & NC lookup tables successfully created for %s model" % (xc_type, model))
    else: # energy loss XC
        material = kwargs['material']
        # data_lst = xc_dict[model]
        brem = xc_dict['brem']
        pair = xc_dict['pair']
        pn_bb = xc_dict['pn_bb']
        pn = xc_dict['pn']
        dframe = pd.DataFrame({'energy':np.logspace(0,12,121,base=10), 'sigma_brem':brem, 'sigma_pair':pair, 'sigma_pn_bb':pn_bb, 'sigma_pn':pn})
        hdf.put('Energy_Loss/%s/%s/xc' % (xc_type,material),dframe)
        hdf.close()
        return print("%s_sigma lookup table successfully created for brem, pair, pn_bb & pn in %s" % (xc_type, material))
    return None

def get_xc(xc_type, model, **kwargs):
    if xc_type=='nu':
        particle = kwargs['particle']
        if particle=='anti-neutrino':particle='anti_neutrino'
        try:
            dataset_xc = pd.read_hdf('lookup_tables.h5','Neutrino_Cross_Sections/%s/xc/%s' % (particle,model))
            cscc = dataset_xc.sigma_cc
            csnc = dataset_xc.sigma_nc
            return np.asarray([cscc,csnc])
        except KeyError:
            model = str(input(("Error finding cross-section values for %s model, please enter a valid model name: " % model)))
            return None
    else: # energy loss; xc_type == 'tau' or 'muon'
        try:
            material = kwargs['material']
            dataset_xc = pd.read_hdf('lookup_tables.h5','Energy_Loss/%s/%s/xc' % (xc_type,material))
            cs_brem = dataset_xc.sigma_brem
            cs_pair = dataset_xc.sigma_pair
            if model == 'BB' or model == 'bb': # Bezrukov Bugaev PN energy loss
                cs_pn = dataset_xc.sigma_pn_bb
            else: # ALLM/custom model for PN energy loss
                cs_pn = dataset_xc.sigma_pn

            return np.asarray([cs_brem,cs_pair,cs_pn])
        except KeyError:
            model = str(input(("Error finding alpha & cross-section values for PN %s model, please enter a valid model name: " % model)))
            return None
    return None

def add_ixc(ixc_type, ixc_dict, **kwargs):
    '''
    Parameters
    ----------
    model : string
        Name of the model you want to add the integrated cross-section values of, to the hdf file
    cross_dict : dictionary
        If ixc_type='nu', this dictionary should contain 4 keys: 'nucc', 'nunc', 'anucc' & 'anunc' and their corresponding keys as dictionaries, and if ixc_type='tau' or 'muon', this dictionary should contain 4 keys: 'brem', 'pair', 'pn_bb' & 'pn' and their corresponding keys as dictionaries

    Returns
    -------
    None.

    '''
    hdf = HDFStore('lookup_tables.h5','a')
    model = kwargs['model']
    if ixc_type == 'nu':
        particle_current = ['anucc','anunc','nucc','nunc']
        for particle in particle_current:
            if 'a' in particle: # anti-neutrino group
                hdf.put('Neutrino_Cross_Sections/anti_neutrino/ixc/%s_%s' % (str(model),particle[3:]),ixc_dict[particle])
            else: # neutrino group
                hdf.put('Neutrino_Cross_Sections/neutrino/ixc/%s_%s' % (str(model),particle[2:]),ixc_dict[particle])
        hdf.close()
        return print("%s_sigma CDF CC & NC lookup tables successfully created for %s model" % (ixc_type, model))

    else: # energy_loss; ixc_type == 'muon' or 'tau'
        material = kwargs['material']
        hdf.put('Energy_Loss/%s/%s/ixc/%s' % (ixc_type,material,model),ixc_dict[model])
        hdf.close()
        return print("%s_sigma CDF lookup table successfully created for %s model in %s" % (ixc_type, model, material))
    return None

def get_ixc(ixc_type, model, **kwargs):
    if ixc_type == 'nu':
        particle = kwargs['particle']
        if particle=='anti-neutrino':particle='anti_neutrino'
        try:
            dataset_ixc_cc = pd.read_hdf('lookup_tables.h5','Neutrino_Cross_Sections/%s/ixc/%s_cc' % (particle,model))
            dataset_ixc_nc = pd.read_hdf('lookup_tables.h5','Neutrino_Cross_Sections/%s/ixc/%s_nc' % (particle,model))
            icc = dataset_ixc_cc.to_dict()
            inc = dataset_ixc_nc.to_dict()

            ixc = {'cc':icc,'nc':inc}
            return ixc
        except KeyError or TypeError:
            model = str(input(("Error finding integrated cross-section values for %s model, please enter a valid model name." % str(model))))
            return None
    else: # energy loss; ixc_type == 'tau' or 'muon'
        try:
            material = kwargs['material']
            dataset_ixc_brem = pd.read_hdf('lookup_tables.h5','Energy_Loss/%s/%s/ixc/brem' % (ixc_type,material))
            dataset_ixc_pair = pd.read_hdf('lookup_tables.h5','Energy_Loss/%s/%s/ixc/pair' % (ixc_type,material))
            if model == 'BB' or model == 'bb': # Bezrukov Bugaev PN energy loss
                dataset_ixc_pn = pd.read_hdf('lookup_tables.h5','Energy_Loss/%s/%s/ixc/pn_bb' % (ixc_type,material))
            else: # ALLM/custom model for PN energy loss
                dataset_ixc_pn = pd.read_hdf('lookup_tables.h5','Energy_Loss/%s/%s/ixc/pn' % (ixc_type,material))

            ixc_brem = dataset_ixc_brem.to_dict()
            ixc_pair = dataset_ixc_pair.to_dict()
            ixc_pn = dataset_ixc_pn.to_dict()
            ixc = {'brem':ixc_brem, 'pair':ixc_pair, 'pn':ixc_pn}
            return ixc
        except KeyError or TypeError:
            model = str(input("Error finding energy loss cross-section values for %s model, please enter a valid model name: " % str(model)))
            return None
    return None

def get_nu_ixc(model, **kwargs):
    particle = kwargs['particle']
    if particle=='anti-neutrino':particle='anti_neutrino'
    try:
        dataset_ixc_cc = pd.read_hdf('lookup_tables.h5','Neutrino_Cross_Sections/%s/ixc/%s_cc' % (particle,model))
        dataset_ixc_nc = pd.read_hdf('lookup_tables.h5','Neutrino_Cross_Sections/%s/ixc/%s_nc' % (particle,model))
        return {'cc':dataset_ixc_cc, 'nc':dataset_ixc_nc}
    except KeyError or TypeError:
        model = str(input(("Error finding integrated cross-section values for %s model, please enter a valid model name." % str(model))))
        return None

def get_lep_ixc(model, material):
    dataset_ixc_brem = pd.read_hdf('lookup_tables.h5','Energy_Loss/%s/%s/ixc/brem' % ('tau',material))
    dataset_ixc_pair = pd.read_hdf('lookup_tables.h5','Energy_Loss/%s/%s/ixc/pair' % ('tau',material))
    if model == 'BB' or model == 'bb': # Bezrukov Bugaev PN energy loss
        dataset_ixc_pn = pd.read_hdf('lookup_tables.h5','Energy_Loss/%s/%s/ixc/pn_bb' % ('tau',material))
    else: # ALLM/custom model for PN energy loss
        dataset_ixc_pn = pd.read_hdf('lookup_tables.h5','Energy_Loss/%s/%s/ixc/pn' % ('tau',material))
        return {'brem':dataset_ixc_brem, 'pair':dataset_ixc_pair, 'pn':dataset_ixc_pn}
        # return dataset_ixc_brem, dataset_ixc_pair, dataset_ixc_pn
    return None

def add_alpha(alpha, **kwargs):
    particle = kwargs['particle']
    material = kwargs['material']
    hdf = HDFStore('lookup_tables.h5','a')
    alpha_df = pd.DataFrame({'energy':np.logspace(0,12,121,base=10),'alpha':alpha})
    hdf.put('Energy_Loss/%s/%s/alpha' % (particle,material),alpha_df)
    hdf.close()
    return print("%s_alpha lookup table successfully created in %s" % (particle,material))

def get_alpha(particle, material):
    alpha_df = pd.read_hdf('lookup_tables.h5','Energy_Loss/%s/%s/alpha' % (particle,material))
    alpha_arr = alpha_df.alpha
    return np.asarray(alpha_arr)

def add_beta(beta_dict, **kwargs):
    particle = kwargs['particle']
    material = kwargs['material']
    beta_type = kwargs['beta_type'] # needs to be 'continuous' or 'total'

    brem = beta_dict['brem']
    pair = beta_dict['pair']
    pn_bb = beta_dict['pn_bb']
    pn = beta_dict['pn']

    hdf = HDFStore('lookup_tables.h5','a')
    beta_df = pd.DataFrame({'energy':E_lep,'beta_brem':brem,'beta_pair':pair,'beta_pn_bb':pn_bb,'beta_pn':pn})
    hdf.put('Energy_Loss/%s/%s/beta_%s' % (particle,material,beta_type),beta_df)
    hdf.close()
    return print("%s_beta_%s lookup table successfully created in %s" % (particle,beta_type,material))

def get_beta(particle, material, beta_type, pn_model):

    beta_df = pd.read_hdf('lookup_tables.h5','Energy_Loss/%s/%s/beta_%s' % (particle,material,beta_type))

    brem = beta_df.beta_brem
    pair = beta_df.beta_pair

    if pn_model=='bb' or pn_model=='BB':
        pn = beta_df.beta_pn_bb
    else:
        pn = beta_df.beta_pn

    return np.array([brem, pair, pn])


def add_pexit(energy_val, prob_dict):
    energy_str = str("%.0e" % energy_val).replace("+",'')
    angle = prob_dict['angle']
    no_regen = prob_dict['no_regen']
    regen = prob_dict['regen']

    hdf = HDFStore('output.h5','a')
    prob_df = pd.DataFrame({'angle':angle, 'no_regen':no_regen,'regen':regen})
    prob_df.set_index("angle", inplace = True)

    hdf.put('Exit_Probability/%s' % energy_str,prob_df)
    hdf.close()
    return None

def get_pexit(energy_val, p_type='regen'):
    energy_str = str("%.0e" % energy_val).replace("+",'')
    p_exit = pd.read_hdf('output.h5','Exit_Probability/%s' % energy_str)
    no_regen = dict(zip(p_exit.index, p_exit.no_regen))
    regen = dict(zip(p_exit.index, p_exit.regen))

    if p_type == 'no_regen':
        return no_regen
    elif p_type == 'regen':
        return regen
    return "Error in get_prob in data"

def add_lep_out(energy_val, angle_val, lep_dict):
    energy_str = str("%.0e" % energy_val).replace("+",'')
    lep_energies = lep_dict["lep_energy"]

    hdf = HDFStore('output.h5','a')
    lep_df = pd.DataFrame({'lep_energy':lep_energies})

    hdf.put('Lep_out_energies/%s/%d' % (energy_str,angle_val),lep_df)
    hdf.close()
    return None

def get_lep_out(energy_val, angle_val):
    energy_str = str("%.0e" % energy_val).replace("+",'')
    e_out = pd.read_hdf('output.h5','Lep_out_energies/%s/%s' % (energy_str,angle_val))
    no_cdf = np.asarray(e_out.lep_energy)
    return no_cdf

# def make_2d_lerp(ixc):
#     v2 = -np.linspace(0.1,3,num=30).astype(np.float64)
#     v2 = np.insert(v2,0,0)
#     x = np.tile(ixc.columns.to_numpy(), (31, 1))
#     y = ixc.to_numpy()
#     z = p.tile(v2, (y.shape[1], 1)).T
#     print(x.shape, y.shape, z.shape)
#     lerp = interpolate.interp2d(x, y, z)
#     return lerp

def make_lerps(ixc):
    v2 = np.linspace(0,-3,num=31).astype(np.float64)
    lerps = [interpolate.interp1d(ixc[e].to_numpy(), v2) for e in ixc.columns]
    return lerps

def get_nu_iixc(model, **kwargs):
    ixc_nu = get_nu_ixc(model, **kwargs)
    ixc_nu['cc'] = make_lerps(ixc_nu['cc'])
    ixc_nu['nc'] = make_lerps(ixc_nu['nc'])
    return ixc_nu

def get_lep_iixc(model, material):
    ixc_lep = get_lep_ixc(model, material)
    ixc_lep['brem'] = make_lerps(ixc_lep['brem'])
    ixc_lep['pair'] = make_lerps(ixc_lep['pair'])
    ixc_lep['pn'] = make_lerps(ixc_lep['pn'])
    return ixc_lep


# =============================================================================
# Test
# =============================================================================
if __name__ == "__main__":
    arr = get_beta('tau', 'rock', 'total', 'allm', True)
