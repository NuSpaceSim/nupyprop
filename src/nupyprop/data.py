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
# import collections
import time
import os

try:
    import importlib.resources as  importlib_resources
except:
    import importlib_resources

# data_dir = '/home/sam/nupyprop_test/output'
# data_dir = '/home/sam/nupyprop_test'

E_nu = np.logspace(3,12,91,base=10).astype(np.float64)
E_lep = np.logspace(0,12,121,base=10).astype(np.float64)

def sci_str(exp_value):
    dec = Decimal(exp_value)
    str_val = ('{:.' + str(len(dec.normalize().as_tuple().digits) - 1) + 'e}').format(dec).replace('+', '')
    return str_val

def chk_file(nu_type, lepton,idepth,nu_cs,lep_pn,loss_type,stats):
    # print('inside chk_file')
    # if nu_type == 'neutrino':nu_type = 'nu'
    # if nu_type == 'antii-neutrino':nu_type = 'anu'
    idepth_str = str(idepth) + 'km'
    stats_str = sci_str(stats)
    if os.path.exists("output_%s_%s_%s_%s_%s_%s_%s.h5" % (nu_type,lepton,idepth_str,nu_cs,lep_pn,loss_type,stats_str)):
        y_n = input('There already exists a output file with these set of parameters (output_%s_%s_%s_%s_%s_%s_%s.h5). Do you wish to overwrite it? Press \'y\' for yes and \'n\' for no: ' % (nu_type,lepton,idepth_str,nu_cs,lep_pn,loss_type,stats_str))
        # if y_n not in {"y", "n"}:print("please enter valid input")
        if y_n not in {"y", "n"}:
            print("Invalid option. Please press \'y\' for yes and \'n\' for no")
            return chk_file(nu_type, lepton,idepth,nu_cs,lep_pn,loss_type,stats)
        elif y_n == 'y':return 1
        elif y_n == 'n':return 0
        # else:
        #     print('Invalid option. Please try again')
            # return chk_file(lepton,idepth,nu_cs,lep_pn,loss_type,stats)
    else:
        return 1 # output file non existant or overwrite enabled
    # return out_val


def add_trajs(type_traj, idepth, traj_array):
    ref = importlib_resources.files('nupyprop.datafiles') / 'lookup_tables.h5'
    with importlib_resources.as_file(ref) as lookup_tables:
        hdf = HDFStore(lookup_tables,'a')
        if type_traj == 'col':branch = 'Column_Trajectories' # sub-sub branch inside the Earth/traj_idepth branch
        elif type_traj == 'water':branch = 'Water_Trajectories'
        hdf.put('Earth/traj_%s/%s' % (str(idepth),branch), traj_array, format='t', data_columns=True)
        hdf.close()
        return print("%s lookup table successfully created for idepth = %s" % (branch,str(idepth)))

def get_trajs(type_traj, beta, idepth): # returns {xalong:cdalong} for beta if type=col or returns chord, water for beta if type=water
    if type_traj == 'col':
        ref = importlib_resources.files('nupyprop.datafiles') / 'lookup_tables.h5'
        with importlib_resources.as_file(ref) as lookup_tables:
            dataset = pd.read_hdf(lookup_tables,'Earth/traj_%s/Column_Trajectories' % str(idepth))
            dataset_sliced = dataset[dataset['beta']==beta]
            xalong = np.asfortranarray(dataset_sliced.xalong.T)
            cdalong = np.asfortranarray(dataset_sliced.cdalong.T)
            return xalong, cdalong

    elif type_traj == 'water':
        ref = importlib_resources.files('nupyprop.datafiles') / 'lookup_tables.h5'
        with importlib_resources.as_file(ref) as lookup_tables:
            dataset = pd.read_hdf(lookup_tables,'Earth/traj_%s/Water_Trajectories' % str(idepth))
            chord = float(dataset.chord[dataset['beta']==beta])
            water = float(dataset.water[dataset['beta']==beta])
            return chord, water
    return "Error in get_trajs in Data"

def add_xc(part_type, xc_obj, model, **kwargs):
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
    ref = importlib_resources.files('nupyprop.datafiles') / 'lookup_tables.h5'
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

def get_xc(part_type, model, **kwargs):
    if part_type=='nu':
        particle = kwargs['particle']
        if particle=='anti-neutrino':particle='anti_neutrino'
        try:
            ref = importlib_resources.files('nupyprop.datafiles') / 'lookup_tables.h5'
            with importlib_resources.as_file(ref) as lookup_tables:
                dataset_xc = pd.read_hdf(lookup_tables,'Neutrino_Cross_Sections/%s/xc/%s' % (particle,model))
                cscc = dataset_xc.sigma_cc
                csnc = dataset_xc.sigma_nc
                out_arr = np.asarray([cscc,csnc])
                return np.asfortranarray(out_arr.T)
        except KeyError:
            model = str(input(("Error finding cross-section values for %s model, please enter a valid model name: " % model)))
            return None
    else: # energy loss; part_type == 'tau' or 'muon'
        try:
            material = kwargs['material']

            ref = importlib_resources.files('nupyprop.datafiles') / 'lookup_tables.h5'
            with importlib_resources.as_file(ref) as lookup_tables:
                dataset_xc = pd.read_hdf(lookup_tables,'Energy_Loss/%s/%s/xc/%s' % (part_type,material,model))

            out_arr = np.asarray(dataset_xc['sigma_%s' % model])
            return np.asfortranarray(out_arr.T)

        except KeyError:
            model = str(input(("Error finding cross-section values for %s model, please enter a valid model name: " % model)))
            return None
    return None

def add_ixc(part_type, ixc_dict, model, **kwargs):
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
    ref = importlib_resources.files('nupyprop.datafiles') / 'lookup_tables.h5'
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

def get_ixc(part_type, model, **kwargs):
    if part_type == 'nu':
        particle = kwargs['particle']
        if particle=='anti-neutrino':particle='anti_neutrino'
        try:
            ref = importlib_resources.files('nupyprop.datafiles') / 'lookup_tables.h5'
            with importlib_resources.as_file(ref) as lookup_tables:
                dataset_ixc_cc = pd.read_hdf(lookup_tables,'Neutrino_Cross_Sections/%s/ixc/%s_cc' % (particle,model))
                dataset_ixc_nc = pd.read_hdf(lookup_tables,'Neutrino_Cross_Sections/%s/ixc/%s_nc' % (particle,model))

                ixc_cc, ixc_nc = [], []

                for energy in range(len(E_nu)):
                    ixc_cc.append(np.asarray(dataset_ixc_cc[E_nu[energy]]))
                    ixc_nc.append(np.asarray(dataset_ixc_nc[E_nu[energy]]))

                ixc_cc = np.asarray(ixc_cc)
                ixc_nc = np.asarray(ixc_cc)

                out_arr = np.asarray([ixc_cc, ixc_nc])
                return np.asfortranarray(out_arr.T)

        except KeyError or TypeError:
            model = str(input(("Error finding integrated cross-section values for %s model, please enter a valid model name." % str(model))))
            return None
    else: # energy loss; ixc_type == 'tau' or 'muon'
        try:
            material = kwargs['material']
            ref = importlib_resources.files('nupyprop.datafiles') / 'lookup_tables.h5'
            with importlib_resources.as_file(ref) as lookup_tables:
                dataset_ixc = pd.read_hdf(lookup_tables,'Energy_Loss/%s/%s/ixc/%s' % (part_type,material,model))


            ixc = []
            for energy in range(len(E_lep)):
                ixc.append(np.asarray(dataset_ixc[E_lep[energy]]))

            ixc = np.asarray(ixc)

            out_arr = ixc
            return np.asfortranarray(out_arr)

        except KeyError or TypeError:
            model = str(input("Error finding energy loss cross-section values for %s model, please enter a valid model name: " % str(model)))
            return None
    return None

def add_alpha(alpha, particle, material):
    ref = importlib_resources.files('nupyprop.datafiles') / 'lookup_tables.h5'
    with importlib_resources.as_file(ref) as lookup_tables:
        hdf = HDFStore(lookup_tables,'a')
        alpha_df = pd.DataFrame({'energy':E_lep,'alpha':alpha})
        hdf.put('Energy_Loss/%s/%s/alpha' % (particle,material),alpha_df, format='t', data_columns=True)
        hdf.close()
        return print("%s_alpha lookup table successfully created in %s" % (particle,material))

def get_alpha(particle, material):
    ref = importlib_resources.files('nupyprop.datafiles') / 'lookup_tables.h5'
    with importlib_resources.as_file(ref) as lookup_tables:
        alpha_df = pd.read_hdf(lookup_tables,'Energy_Loss/%s/%s/alpha' % (particle,material))
        alpha_arr = alpha_df.alpha
        out_arr = np.asarray(alpha_arr)
        return np.asfortranarray(out_arr.T)

def add_beta(beta_arr, particle, material, model, beta_type):
    ref = importlib_resources.files('nupyprop.datafiles') / 'lookup_tables.h5'
    with importlib_resources.as_file(ref) as lookup_tables:
        hdf = HDFStore(lookup_tables,'a')
        beta_df = pd.DataFrame({'energy':E_lep, 'beta_%s' % model:beta_arr})
        hdf.put('Energy_Loss/%s/%s/beta_%s/%s' % (particle,material,beta_type,model), beta_df, format='t', data_columns=True)
        hdf.close()
        return print("%s_beta_%s lookup table successfully created in %s" % (particle,beta_type,material))

def get_beta(particle, material, model, beta_type):
    ref = importlib_resources.files('nupyprop.datafiles') / 'lookup_tables.h5'
    with importlib_resources.as_file(ref) as lookup_tables:
        beta_df = pd.read_hdf(lookup_tables,'Energy_Loss/%s/%s/beta_%s/%s' % (particle,material,beta_type,model))
        beta_arr = beta_df['beta_%s' % model]
        out_arr = np.asarray(beta_arr)
        return np.asfortranarray(out_arr.T)

def combine_lep(data_type, particle, material, pn_model, **kwargs):

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

def add_pexit(nu_type, lepton, energy_val, prob_dict, idepth, nu_cs, lep_pn, loss_type, stats):

    log_energy = np.log10(energy_val)
    # energy_str = str("%.0e" % energy_val).replace("+",'')
    energy_str = str(log_energy)
    idepth_str = str(idepth) + 'km'
    stats_str = sci_str(stats)

    angle = prob_dict['angle']
    no_regen = prob_dict['no_regen']
    regen = prob_dict['regen']

    hdf = HDFStore('output_%s_%s_%s_%s_%s_%s_%s.h5' % (nu_type,lepton,idepth_str,nu_cs,lep_pn,loss_type,stats_str),'a')
    prob_df = pd.DataFrame({'angle':angle, 'no_regen':no_regen,'regen':regen})
    prob_df.set_index("angle", inplace = True)

    hdf.put('Exit_Probability/%s' % energy_str,prob_df, format='t', data_columns=True)
    hdf.close()
    return None

def get_pexit(nu_type, lepton, energy_val, p_type, idepth, nu_cs, lep_pn, loss_type, stats):
    # os.chdir(data_dir)
    log_energy = np.log10(energy_val)
    # energy_str = str("%.0e" % energy_val).replace("+",'')
    energy_str = str(log_energy)
    idepth_str = str(idepth) + 'km'
    stats_str = sci_str(stats)

    p_exit = pd.read_hdf('output_%s_%s_%s_%s_%s_%s_%s.h5' % (nu_type,lepton,idepth_str,nu_cs,lep_pn,loss_type,stats_str),'Exit_Probability/%s' % energy_str)
    no_regen = dict(zip(p_exit.index, p_exit.no_regen))
    regen = dict(zip(p_exit.index, p_exit.regen))

    if p_type == 'no_regen':
        return no_regen
    elif p_type == 'regen':
        return regen
    return "Error in get_prob in data"

def add_lep_out(nu_type, lepton, energy_val, angle_val, lep_dict, idepth, nu_cs, lep_pn, loss_type, stats):
    log_energy = np.log10(energy_val)
    # energy_str = str("%.0e" % energy_val).replace("+",'')
    energy_str = str(log_energy)
    idepth_str = str(idepth) + 'km'
    stats_str = sci_str(stats)

    lep_energies = lep_dict["lep_energy"]

    hdf = HDFStore('output_%s_%s_%s_%s_%s_%s_%s.h5' % (nu_type,lepton,idepth_str,nu_cs,lep_pn,loss_type,stats_str),'a')
    try:
        if len(lep_energies) > 1:
            lep_df = pd.DataFrame({'lep_energy':lep_energies})
        elif len(lep_energies) == 1:
            lep_df = pd.DataFrame({'lep_energy':lep_energies}, index=[0])
        else:
            lep_df = pd.DataFrame({'lep_energy':np.array([0])}, index=[0])
    except TypeError or ValueError:
        lep_df = pd.DataFrame({'lep_energy':np.array([0])}, index=[0])

    hdf.put('Lep_out_energies/%s/%d' % (energy_str,angle_val),lep_df, format='t', data_columns=True)
    hdf.close()
    return None

def get_lep_out(nu_type, lepton, energy_val, angle_val, idepth, nu_cs, lep_pn, loss_type, stats):
    # os.chdir(data_dir)
    log_energy = np.log10(energy_val)
    # energy_str = str("%.0e" % energy_val).replace("+",'')
    energy_str = str(log_energy)
    idepth_str = str(idepth) + 'km'
    stats_str = sci_str(stats)

    e_out = pd.read_hdf('output_%s_%s_%s_%s_%s_%s_%s.h5' % (nu_type,lepton,idepth_str,nu_cs,lep_pn,loss_type,stats_str),'Lep_out_energies/%s/%s' % (energy_str,angle_val))
    out_lep = np.asarray(e_out.lep_energy)
    return out_lep

def add_pexit_manual(nu_type, energy_val, angles, idepth, nu_cs, lep_pn, loss_type, stats): # manual will only work for regen (so basically, for muons)
    log_energy = np.log10(energy_val)
    lepton = 'muon'
    # energy_str = str("%.0e" % energy_val).replace("+",'')
    energy_str = str(log_energy)
    idepth_str = str(idepth) + 'km'
    stats_str = sci_str(stats)

    angle, no_regen, regen = [],[],[]
    for angle_val in angles:
        e_out = pd.read_hdf('output_%s_%s_%s_%s_%s_%s_%s.h5' % (nu_type,lepton,idepth_str,nu_cs,lep_pn,loss_type,stats_str),'Lep_out_energies/%s/%s' % (energy_str,angle_val))
        e_out_len =  len(e_out.lep_energy)
        if e_out_len == 1:
            if np.asarray(e_out.lep_energy)[0] == 0:
                p_exit = 0.0
            else:
                p_exit = e_out_len/stats
        else:
            p_exit = e_out_len/stats
        angle.append(angle_val)
        no_regen.append(p_exit)
        regen.append(p_exit)

    hdf = HDFStore('output_%s_%s_%s_%s_%s_%s_%s.h5' % (nu_type,lepton,idepth_str,nu_cs,lep_pn,loss_type,stats_str),'a')
    prob_df = pd.DataFrame({'angle':angle, 'no_regen':no_regen,'regen':regen})
    prob_df.set_index("angle", inplace = True)

    hdf.put('Exit_Probability/%s' % energy_str,prob_df, format='t', data_columns=True)
    hdf.close()
    return None

def add_cdf(energies, angles, nu_type, lepton, idepth, nu_cs, lep_pn, loss_type, stats):
    # os.chdir(data_dir)

    idepth_str = str(idepth) + 'km'
    stats_str = sci_str(stats)

    # energies = np.logspace(7,11,17)
    # angles = np.arange(1,36)

    for energy in energies:
        # print('energy = ', energy)
        log_energy = np.log10(energy)
        # energy_str = str("%.0e" % energy_val).replace("+",'')
        energy_str = str(log_energy)
        for angle in angles:
            # print("angle = ", angle)
            lep_out = get_lep_out(nu_type, lepton, energy, angle, idepth, nu_cs, lep_pn, loss_type, stats)
            hdf = HDFStore('output_%s_%s_%s_%s_%s_%s_%s.h5' % (nu_type,lepton,idepth_str,nu_cs,lep_pn,loss_type,stats_str),'a')
            bins = np.logspace(-5,0,51)
            z = lep_out/energy
            binned_z = np.digitize(z, bins)
            bin_counts = np.bincount(binned_z)
            if len(bin_counts) < len(bins):
                zeros_z = np.zeros(len(bins) - len(bin_counts))
                binned_z_fixed = np.concatenate((bin_counts, zeros_z))
            else:
                binned_z_fixed = bin_counts
            z_cumsum = np.cumsum(binned_z_fixed)
            # z_cdf = binned_z_fixed/max(binned_z_fixed)
            z_cdf = z_cumsum/z_cumsum[-1]
            z_cdf[0] = 0
            # print("len_z = ",len(bins))
            # print("len_z_cdf = ",len(z_cdf))
            z_cdf_df = pd.DataFrame({'z':bins, 'cdf':z_cdf})
            # z_cdf_df.set_index("z", inplace = True)
            hdf.put('Lep_out_cdf/%s/%d' % (energy_str,angle),z_cdf_df, format='t', data_columns=True)
            hdf.close()
    return None

# =============================================================================
# Test
# =============================================================================
if __name__ == "__main__":
    # arr = get_beta('tau', 'rock', 'total', 'allm', True)
    # add_pexit_manual(1e9, np.arange(1,36), 1e8)
    pass
    # pexit = get_pexit('tau', 1e7, p_type='regen', loss_type='stochastic')
    # add_cdf('tau', 'full')
