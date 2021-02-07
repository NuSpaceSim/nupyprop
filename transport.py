#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 19:31:14 2020

@author: sam
"""
import data as Data
import geometry as Geometry
import energy_loss as Energy_loss
import my_interpolation as Interpolation

import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd
import scipy.constants as scc
import time
from numba import njit
# import timeit
# from numba.typed import Dict
# import numba as nb

E_nu = Data.E_nu
E_lep = Data.E_lep
Re = Geometry.Re
N_A = scc.Avogadro

rho_water = 1.02 # g/cm^3
rho_rock = 2.65 # g/cm^3

v2 = -np.linspace(0.1,3,num=30).astype(np.float64)
v2 = np.insert(v2,0,0) # padding 0 at the beginning to match index with ixc entry nos. as those start with 1


# @njit(nogil=True)
def my_rand(): # because we need a random no. in the range of (0,1) and not [0,1)
    random_no = random.random()
    while random_no == 0:
        random_no = random.random()
    return random_no
    # return 0.33 # for debugging only!

# @njit(nogil=True)
def int_length_nu(energy, nu_xc): # interaction lengths for neutrino & leptons; in cm
    sig_cc, sig_nc = Interpolation.int_xc_nu(energy, nu_xc) # initialize CC & NC xc interpolations
    x_int = 1/(((sig_cc + sig_nc))*fac_nu) # check the * or / fac_nu!
    return x_int

# @njit(nogil=True)
def int_length_lep(energy, xc_arr, rho): # interaction lengths for neutrino & leptons; in cm
    sig_cc = 0 # placeholder for CC lepton interactions
    sig_nc = 0 # placeholder for NC lepton interactions
    gamma = energy/m_le
    sig_brem, sig_pair, sig_pn = Interpolation.int_xc_lep(energy, xc_arr) # initialize brem, pair & pn xc interpolations
    x_int = 1/((sig_brem + sig_pair + sig_pn + (1/(gamma*c_tau*rho)) + sig_cc + sig_nc))
    return x_int

# @njit(nogil=True)
def interaction_type_nu(energy, nu_xc):
    sig_cc, sig_nc = Interpolation.int_xc_nu(energy, nu_xc)
    tot_frac = 1/int_length_nu(energy, nu_xc)

    cc_frac = sig_cc/tot_frac

    x = np.random.uniform(0.0, 1.0, len(energy))
    # if x<=cc_frac:
    #     return 'cc' # interaction type
    # else:
    #     return 'nc' # interaction type
    return np.where(x <= cc_frac, 'cc', 'nc')

# interaction_type_strings = ('decay', 'brem', 'pair', 'pn', 'lep_nc')
# @njit(nogil=True)
def interaction_type_lep(energy, xc_arr, rho):
    sig_brem,sig_pair,sig_pn = Interpolation.int_xc_lep(energy, xc_arr)
    gamma = energy/m_le
    # tot_frac = 1/int_length_lep(energy, xc_arr, rho)
    tot_frac = (sig_brem + sig_pair + sig_pn + (1/(gamma*c_tau*rho)))
    # decay_frac = (1/(gamma*c_tau*rho))/tot_frac
    # brem_frac = sig_brem/tot_frac
    # pair_frac = sig_pair/tot_frac
    # pn_frac = sig_pn/tot_frac

    # x = my_rand()
    # x = 0.33 # for debugging only!
    x = np.random.uniform(0.0, 1.0, len(energy))
    part_type = np.empty(len(energy), dtype="object")

    targ_frac = (1/(gamma*c_tau*rho))/tot_frac
    mask = x<targ_frac
    part_type[mask] = 'decay'
    filled_mask = np.copy(mask)

    targ_frac += sig_brem / tot_frac
    mask = ~filled_mask & (x < targ_frac)
    part_type[mask] = 'brem'
    filled_mask |= mask

    targ_frac += sig_pair / tot_frac
    mask = ~filled_mask & (x < targ_frac)
    part_type[mask] = 'pair'
    filled_mask |= mask

    targ_frac += sig_pn / tot_frac
    mask = ~filled_mask & (x < targ_frac)
    part_type[mask] = 'pn'
    filled_mask |= mask
    # if x<targ_frac:
    #     return 'pn' # ip = 3
    part_type[~filled_mask] = 'lep_nc'

    return part_type

# @njit(nogil=True)
def find_y(energy, ixc_dict, ip):
    '''
    Piecewise interpolation over masked CDF tables for given energies and
    ip types.
    '''
    dlv = np.empty_like(energy)
    for j in ixc_dict.keys():
        lerps = ixc_dict[j]
        E__ = E_nu if j == 'cc' or j == 'nc' else E_lep
        ipmask = ip == j
        idxs = np.searchsorted(E__, energy[ipmask])
        djlv = np.empty(np.count_nonzero(ipmask))

        for i, lerp in enumerate(lerps):
            msk = idxs == i
            sz = np.count_nonzero(msk)
            dy = np.random.uniform(0.0, 1.0, sz)
            lv = lerp(dy)
            djlv[msk] = lv

        dlv[ipmask] = djlv

    y = 10**dlv
    y[y > 1.0] = 1.0
    return y


# @njit(nogil=True)
def full_sub_mask(full_mask, sub_mask):
    fsmask = np.copy(full_mask)
    newmsk = np.copy(full_mask)
    fsmask[fsmask] = sub_mask
    newmsk[newmsk] = ~sub_mask
    return fsmask, newmsk


# =============================================================================
#                 Nu propagation
# =============================================================================
# @njit(nogil=True)
def propagate_nu(e_init, nu_xc, nu_ixc, depth_max):

    part_type = np.full(len(depth_max), 'nu', dtype="object") # starting off as a neutrino
    col_depth_total = depth_max*1e5 # in cm.w.e
    e_nu = np.full_like(depth_max, e_init)
    e_fin = np.full_like(depth_max, e_init)
    x_0 = np.zeros_like(depth_max) # starting depth in cm.w.e
    d_travel = np.copy(depth_max) # added this in case there is a problem and needed for breaking out when E<1e3

    mask = np.full_like(depth_max, True, dtype=np.bool)

    # 10 continue
    while np.any(e_nu[mask]>1e3):

        int_len = int_length_nu(e_nu[mask], nu_xc) # get the interaction length
        dy = np.random.uniform(0.0, 1.0, len(int_len))
        x = -int_len*np.log(dy) # logarithmic sampling
        x_f = x_0[mask] + x

        cdmsk = (x_f > col_depth_total[mask])
        mask[mask] = ~cdmsk
        x_0[mask] = x_f[~cdmsk] # keep going and update x_0

        int_type = interaction_type_nu(e_nu[mask], nu_xc) # what type of interaction is occurring?
        ccmsk, _ = full_sub_mask(mask, (part_type[mask] == 'nu') & (int_type == 'cc'))
        part_type[ccmsk] = 'lep'

        y = find_y(e_nu[mask],nu_ixc,int_type) # how much energy does the neutrino/lepton have after an interaction
        y[y>0.999] = 0.999

        e_fin[mask] = e_nu[mask]*(1-y) # this is the energy for the next interaction

        lep_mask = mask & (part_type == 'lep')
        d_travel[lep_mask] = x_0[lep_mask] * 1e-5
        mask &= ~(part_type == 'lep')

        e_nu[mask] = e_fin[mask] # now we have a nu with energy e_fin, which has gone x_f
        d_travel[mask] = x_0[mask]*1e-5

    return part_type, d_travel, e_fin

# =============================================================================
#                 Tau propagation in water
# =============================================================================
# @njit(nogil=True)
def propagate_lep_water(e_init, xc_water, lep_ixc, alpha_water, beta_water, d_in, type_loss):

    # if type_loss == 'stochastic':

    e_min = 1e3 # minimum tau energy, in GeV

    part_id = np.full(len(d_in), 'not_decayed', dtype="object") # start with tau that's obviously not decayed
    d_fin = np.zeros_like(d_in)

    # cd_left = d_in*1e5 # how much to go, in cm.w.e
    e_lep = np.copy(e_init)
    e_fin = np.copy(e_init) # in case the first interaction is too far
    x_0 = np.zeros_like(e_lep) # haven't gone anywhere yet
    mask = np.full(len(e_lep), True)
    cnt = 0

    # 10 continue
    while np.any(e_lep[mask] > e_min):
        cnt+=1
        # find interaction distance in g/cm^2
        int_len = int_length_lep(e_lep[mask], xc_water, rho_water)
        dy = np.random.uniform(0.0, 1.0, len(int_len))
        x = -int_len*np.log(dy) # t is rho*L and has units g/cm^2
        x_f = x_0[mask] + x # going 1D; how far have we traveled here
        d_fin[mask] = x_f/1e5 # make sure it is not past the old number, in km.w.e
        # print('1', mask)

        # if x_f >= cd_left: # already past maximum depth but still a tau; go to 30
        #     # 30 continue
        #     d_rem = cd_left - x_0
        #     alpha = Interpolation.int_alpha(e_lep, alpha_water) # changed 12/9/2020
        #     beta = Interpolation.int_beta(e_lep, beta_water) # changed 12/9/2020
        #     e_fin = e_lep - (e_lep*beta + alpha)*d_rem
        #     d_fin = d_in

        #     if e_fin > e_init: # sanity check
        #         e_fin = e_init

        #     if e_fin <= e_min: # tau has decayed
        #         d_fin = d_in # just in case; added 12/9/2020
        #         e_fin = e_min
        #         part_id = 'no_count'
        #     # 50 continue
        #     return part_id, d_fin, e_fin
        dmsk = (x_f >= (d_in[mask] * 1e5))
        full_dmsk, new_msk = full_sub_mask(mask, dmsk)

        d_rem = (1e5 * d_in[full_dmsk]) - x_0[full_dmsk]
        alpha = Interpolation.int_alpha(e_lep[full_dmsk], alpha_water)
        beta = Interpolation.int_beta(e_lep[full_dmsk], beta_water)
        e_fin[full_dmsk] = e_lep[full_dmsk] - (e_lep[full_dmsk]*beta + alpha)*d_rem
        d_fin[full_dmsk] = d_in[full_dmsk]

        efhimsk, _ = full_sub_mask(full_dmsk, e_fin[full_dmsk] > e_init[full_dmsk])
        e_fin[efhimsk] = e_init[efhimsk]
        eflomsk, _ = full_sub_mask(full_dmsk, e_fin[full_dmsk] <= e_min)
        e_fin[eflomsk] = e_min
        part_id[eflomsk] = 'no_count'

        mask = new_msk
        x_0[mask] = x_f[~dmsk] # update x_0 and keep going
        x = x[~dmsk]
        x_f = x_f[~dmsk]
        # print(mask)

        alpha = Interpolation.int_alpha(e_lep[mask], alpha_water)
        beta = Interpolation.int_beta(e_lep[mask], beta_water)
        e_int = e_lep[mask] - (e_lep[mask]*beta + alpha)*x # find some intermediate energy to get reasonable values of energy between initial and final energy, a la MUSIC
        # if e_int <= e_min: e_int = e_min
        # e_avg = 10**((np.log10(e_lep)+np.log10(e_int))/2) # does this work?; changed 12/9/2020
        e_int[e_int < e_min] = e_min
        e_avg = 10**((np.log10(e_lep[mask])+np.log10(e_int))/2) # does this work?

        alpha = Interpolation.int_alpha(e_avg, alpha_water)
        beta = Interpolation.int_beta(e_avg, beta_water)

        e_int = Energy_loss.em_cont_part(e_lep[mask], alpha, beta, x) # get the continuous energy loss
        # print(mask)

        # if e_int <= e_min: # below minimum energy; go to 20
        #     # 20 continue
        #     d_fin = d_in # changed 12/9/2020
        #     e_fin = e_min
        #     part_id = 'no_count' # don't count this
        #     return part_id, d_fin, e_fin
        low_ei_msk = e_int <= e_min
        full_low_msk = np.copy(mask)
        full_low_msk[full_low_msk] = low_ei_msk
        e_fin[full_low_msk] = e_min
        d_fin[full_low_msk] = d_in[full_low_msk]
        part_id[full_low_msk] = 'no_count'

        mask[mask] = ~low_ei_msk
        e_int = e_int[~low_ei_msk]
        x_f = x_f[~low_ei_msk]
        # print(mask)

        # print(e_int)
        int_type = interaction_type_lep(e_int, xc_water, rho_water)
        # print(int_type)

        # if int_type == 'decay': # tau has decayed
        #     part_id = 'decayed'
        #     e_fin = e_int
        #     d_fin = x_f/1e5 # commented out 12/9/2020
        #     # go to 50
        #     # 50 continue
        #     return part_id, d_fin, e_fin # basically what d_fin does this return?
        decay_msk = int_type == 'decay'
        full_decay_msk = np.copy(mask)
        full_decay_msk[full_decay_msk] = decay_msk
        part_id[full_decay_msk] = 'decayed'
        # print(type(int_type), len(e_fin[full_decay_msk]), len(e_init))
        e_fin[full_decay_msk] = e_int[decay_msk]
        d_fin[full_decay_msk] = x_f[decay_msk]

        mask[mask] = ~decay_msk
        int_type = int_type[~decay_msk]
        # print(decay_msk, ~decay_msk, e_int)
        e_int = e_int[~decay_msk]

        # tau didn't decay. Now how much energy does it have after interaction?
        y = find_y(e_int, lep_ixc, int_type)

        # outgoing tau energy is old e_lep*(1-y)
        # e_lep = e_int*(1-y) # this is the energy for the next interaction
        # e_fin = e_lep
        e_lep[mask] = e_int*(1-y) # this is the energy for the next interaction
        e_fin[mask] = e_lep[mask]
        # go to 10

    # # Outside the while loop, e_lep has to be < e_min
    # if e_lep <= e_min: # only continuous energy loss; go to 20
    #     d_fin = d_in
    #     e_fin = e_min
    #     part_id = 'no_count' # don't count this
    #     return part_id, d_fin, e_fin
    outer_mask = np.copy(mask)
    outer_mask[mask] = (e_lep[mask] <= e_min)
    part_id[outer_mask] = 'no_count'
    d_fin[outer_mask] = d_in[outer_mask]
    e_fin[outer_mask] = e_min

    return part_id, d_fin, e_fin

# =============================================================================
#                 Tau propagation in rock
# =============================================================================
# @njit(nogil=True)
def propagate_lep_rock(angle, e_init, xc_rock, lep_ixc, alpha_rock, beta_rock, d_entry, d_in, type_loss, cd2distd, densityatx):
    e_min = 1e3 # minimum tau energy, in GeV

    # Prep return arrays
    part_id = np.full(len(d_entry), 'not_decayed', dtype="object")
    d_fin = np.zeros_like(d_entry)
    e_fin = np.zeros_like(d_entry)

    col_depth = d_entry*1e5 # how far in
    # d_max = d_in*1e5 # how much to go, in cm.w.e
    e_lep = e_init
    x_0 = np.zeros_like(e_lep)
    mask = np.full(len(e_lep), True)
    cnt = 0
    # return part_id, d_fin, e_fin

    while np.any(e_lep[mask] > e_min):
        cnt +=1

        # if e_lep <= e_min: # continuous energy loss only
        #     # go to 20
        #     d_fin = d_max/1e5 # in km.w.e
        #     e_fin = e_min
        #     part_id = 'decayed'
        #     return part_id, d_fin, e_fin
        low_e_msk = e_lep[mask] <= e_min
        full_low_msk = np.copy(mask)
        full_low_msk[full_low_msk] = low_e_msk
        d_fin[full_low_msk] = d_in[full_low_msk]
        e_fin[full_low_msk] = e_min
        part_id[full_low_msk] = 'decayed'

        mask[mask] = ~low_e_msk

        # 10 continue
        x_interp = cd2distd(col_depth[mask]) # find how far we are along the chord for given beta
        _, rho = densityatx(x_interp) # find the density at x

        int_len = int_length_lep(e_lep[mask], xc_rock, rho)

        # dy = my_rand()
        dy = np.random.uniform(0.0, 1.0, len(int_len))
        # dy = 0.33 # for debugging only!
        x = -int_len*np.log(dy)
        col_depth[mask] += x # update along trajectory, from the start of the chord
        x_f = x_0[mask] + x # going 1D
        d_fin[mask] = x_f/1e5 # make sure it is not past the old number, in km.w.e

        # if x_f > d_max: # already past max depth
        #     # go to 30
        #     d_rem = d_max - x_0
        #     alpha = Interpolation.int_alpha(e_lep, alpha_rock)
        #     beta = Interpolation.int_beta(e_lep, beta_rock)
        #     e_fin = e_lep - (e_lep*beta + alpha)*d_rem
        #     d_fin = d_max/1e5
        #     if e_fin > e_init: e_fin=e_init
        #     return part_id, d_fin, e_fin
        dmsk = (x_f > (d_in[mask] * 1e5))
        full_dmsk = np.copy(mask)
        full_dmsk[full_dmsk] = dmsk
        # print(len(e_fin), len(e_lep), len(d_in), len(x_0), len(full_dmsk), len(dmsk))
        d_rem = (1e5 * d_in[full_dmsk]) - x_0[full_dmsk]
        alpha = Interpolation.int_alpha(e_lep[full_dmsk], alpha_rock)
        beta = Interpolation.int_beta(e_lep[full_dmsk], beta_rock)
        d_fin[full_dmsk] = d_in[full_dmsk]
        e_fin[full_dmsk] = e_lep[full_dmsk] - (e_lep[full_dmsk]*beta + alpha)*d_rem
        e_fin_init_msk = e_fin[full_dmsk] > e_init[full_dmsk]
        e_fin[e_fin_init_msk] = e_init[e_fin_init_msk]

        mask[mask] = ~dmsk
        x_0[mask] = x_f[~dmsk] # update x_0 and keep going
        x = x[~dmsk]
        rho = rho[~dmsk]

        alpha = Interpolation.int_alpha(e_lep[mask], alpha_rock)
        beta = Interpolation.int_beta(e_lep[mask], beta_rock)
        e_int = e_lep[mask] - (e_lep[mask]*beta + alpha)*x # find some intermediate energy to get reasonable values of energy between initial and final energy

        # if e_int <= e_min: e_int=e_min
        e_int[e_int < e_min] = e_min

        e_avg = 10**((np.log10(e_lep[mask])+np.log10(e_int))/2) # does this work?

        alpha = Interpolation.int_alpha(e_avg, alpha_rock)
        beta = Interpolation.int_beta(e_avg, beta_rock)

        e_int = Energy_loss.em_cont_part(e_lep[mask], alpha, beta, x) # get the continuous energy

        # if e_int <= e_min: # is it below minimum energy now?
        #     e_fin = e_int
        #     # go to 20
        #     d_fin = d_max/1e5 # in km.w.e
        #     e_fin = e_min
        #     part_id = 'decayed'
        #     return part_id, d_fin, e_fin
        low_ei_msk = e_int <= e_min
        full_low_msk = np.copy(mask)
        full_low_msk[full_low_msk] = low_ei_msk
        e_fin[full_low_msk] = e_min
        d_fin[full_low_msk] = d_in[full_low_msk]
        part_id[full_low_msk] = 'decayed'

        mask[mask] = ~low_ei_msk
        e_int = e_int[~low_ei_msk]
        rho = rho[~low_ei_msk]

        # print(e_int, rho)
        int_type = interaction_type_lep(e_int, xc_rock, rho)
        # print(int_type)

        # if int_type == 'decay': # tau has decayed
        #     part_id = 'decayed'
        #     e_fin = e_int
        #     # go to 50
        #     # 50 continue
        #     return part_id, d_fin, e_fin
        decay_msk = int_type == 'decay'
        part_id[mask][decay_msk] = 'decayed'
        e_fin[mask][decay_msk] = e_int[decay_msk]

        mask[mask] = ~decay_msk
        int_type = int_type[~decay_msk]
        e_int = e_int[~decay_msk]

        # tau didn't decay. Now how much energy does it have after interaction?
        y = find_y(e_int, lep_ixc, int_type)

        # outgoing tau energy is old e_lep*(1-y)
        e_lep[mask] = e_int*(1-y) # this is the energy for the next interaction
        e_fin[mask] = e_lep[mask]

        # if e_lep > e_min: # go back and propagate again
        # go to 10

    # # Outside the while loop, e_lep has to be < e_min
    # if e_lep <= e_min: # only continuous energy loss; go to 20
    #     d_fin = d_max/1e5
    #     e_fin = e_min
    #     part_id = 'decayed' # don't count this
    #     return part_id, d_fin, e_fin
    outer_mask = np.copy(mask)
    outer_mask[mask] = (e_lep[mask] <= e_min)
    part_id[outer_mask] = 'decayed'
    d_fin[outer_mask] = d_in[outer_mask]
    e_fin[outer_mask] = e_min

    return part_id, d_fin, e_fin


# @njit(nogil=True)
def tau_thru_layers(angle, depth, d_water, depth_traj, e_lep_in, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, cd2distd, densityatx): # note: angle is now in angle

    d_fin = np.copy(depth_traj)
    col_depth = depth_traj*1e5 # g/cm^2
    e_fin = np.copy(e_lep_in)
    part_type = np.full(len(depth_traj), 'not_decayed', dtype="object") # tau going in

    # if e_lep_in < 1e3: # just in case
    #     part_type = 'decayed'
    #     d_fin = depth
    #     return part_type, d_fin, e_fin
    # e_lep_in_mask = e_lep_in < 1e3

    # if angle <= 1.5 or depth-depth_traj < d_water:
    #     rho = rho_water # water
    # else:
    #     x = cd2distd(col_depth)
    #     r, rho = densityatx(x)
    #     if rho <= 0: # round off error happening here; went too far
    #         print("col_depth = ", col_depth)
    #         print("x = ", x)
    #         print('rho is 0')
    #     # print('r, rho = ', r, rho)
    #     if rho <= 1.5 and r < 6365.0:
    #         print('rho too small!')
    rho_mask_1 = ~((depth-depth_traj < d_water) | (angle <= 1.5))

    x = cd2distd(col_depth[rho_mask_1])
    _, rho_prime = densityatx(x)

    rho = np.empty_like(col_depth)
    rho[rho_mask_1] = rho_prime
    rho[~rho_mask_1] = rho_water

    # if rho > 1.5: # we aren't in water yet
    #     d_in = depth - depth_traj - d_water # propagate this far in rock

    #     part_type, d_f, e_fin = propagate_lep_rock(angle, e_lep_in, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, depth_traj, d_in, 'stochastic', cd2distd, densityatx)

    #     if part_type == 'not_decayed': # still a tau
    #         e_lep_in = e_fin
    #         d_in = d_water
    #         depth_traj = depth_traj + d_f # now propagate through final layer of water
    #         d_fin = depth_traj
    #         part_type, d_f, e_fin = propagate_lep_water(e_lep_in, xc_water, lep_ixc_water, alpha_water, beta_water, d_in, 'stochastic')

    #     else: # neutrino
    #         return part_type, d_fin, e_fin
    rho_mask_2 = (rho > 1.5)
    d_in = (depth - depth_traj - d_water)[rho_mask_2] # propagate this far in rock

    p_t, d_f, e_f = propagate_lep_rock(angle, e_lep_in[rho_mask_2], xc_rock, lep_ixc_rock, alpha_rock, beta_rock, depth_traj[rho_mask_2], d_in, 'stochastic', cd2distd, densityatx)
    part_type[rho_mask_2] = p_t
    e_fin[rho_mask_2] = e_f

    nd_msk = p_t == 'not_decayed'
    full_nd_msk = np.copy(rho_mask_2)
    full_nd_msk[full_nd_msk] = nd_msk
    d_in[nd_msk] = d_water
    depth_traj[full_nd_msk] += d_f[nd_msk] # now propagate through final layer of water
    d_fin[full_nd_msk] = depth_traj[full_nd_msk]

    p_t, _, e_f = propagate_lep_water(e_f[nd_msk], xc_water, lep_ixc_water, alpha_water, beta_water, d_in[nd_msk], 'stochastic')
    part_type[full_nd_msk] = p_t
    e_fin[full_nd_msk] = e_f

    #     # 202 continue
    # else:
    #     d_in = depth - depth_traj
    #     part_type, d_f, e_fin = propagate_lep_water(e_lep_in, xc_water, lep_ixc_water, alpha_water, beta_water, d_in, 'stochastic')
    #     depth_traj += d_f
    #     return part_type, depth_traj, e_fin
    d_in = (depth - depth_traj)[~rho_mask_2] # propagate this far in rock
    p_t, d_f, e_f = propagate_lep_water(e_lep_in[~rho_mask_2], xc_water, lep_ixc_water, alpha_water, beta_water, d_in, 'stochastic')
    part_type[~rho_mask_2] = p_t
    d_fin[~rho_mask_2] = depth_traj[~rho_mask_2] + d_f
    e_fin[~rho_mask_2] = e_f

    return part_type, d_fin, e_fin

# @njit(nogil=True)
def distnu(r, ithird): # ithird = 1 => 1/3 or ithird = 2 => dist; rename to decay_distnu
    # fnu=lambda y: y/3 * (5 - 3* y**2 + y**3) - y/3 * (1 - 3 * y**2 + 2 * y**3) # polarized
    fnu=lambda y: (4*y - y**4)/3.
    if ithird !=1:
        # fm = 1 # max value of distribution
        # ff = r*fm
        # now solve the equation f(y) = ff
        y0 = np.zeros_like(r)
        y1 = np.ones_like(r)
        yy = np.zeros_like(r)
        msk = np.full(len(r), True)
        # 1
        while np.any(msk):
            # abs(y1-y0) > 0.1e-2:
            # mask[mask] = 
            y = (y0[msk]+y1[msk])/2
            x = fnu(y)
            ltmsk = np.copy(msk)
            gtmsk = np.copy(msk)
            ltmsk[ltmsk] = (x < r)
            gtmsk[gtmsk] = ~(x < r)

            y0[ltmsk] = y[x<r]
            y1[gtmsk] = y[~(x<r)]
            # # if fnu(y) < r:
            # #     y0 = y
            # else:
            #     y1 = y
            yy[msk] = (y0[msk] + y1[msk])/2.
            msk[msk] = np.abs(y1[msk] - y0[msk]) > 0.1e-2

        # return (y0+y1)/2 # solved
        return yy

    else:
        return np.full_like(r, 1/3)

# @njit(nogil=True)
def regen(angle, e_lep, depth, d_water, d_lep, nu_xc, nu_ixc, ithird, xc_water, xc_rock, ixc_water, ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, cd2distd, densityatx):

    # find the neutrino energy from the tau decay with tau energy e_lep
    r = np.random.uniform(0.0, 1.0, len(e_lep))
    frac = distnu(r,ithird)
    e_nu = frac * e_lep
    d_lep = np.copy(d_lep)

    d_left = depth-d_lep # this is how far the neutrino can go
    e_fin = e_nu # in case we need to exit
    part_type = np.full(len(e_lep), 'exit', dtype="object") # in case we need to exit; change later to string (HLS = 2)

    # if d_left <= 0: # past the point of interations allowed
    #     d_exit = depth
    #     # go to 60
    #     # 60 continue
    #     return part_type, d_exit, e_fin

    d_exit = np.copy(d_lep) # we are starting this far into the Earth with a neutrino
    # int_part = 'nu' # starting with a neutrino with energy e_nu; change later to string

    # tnu follows NC to the end, or gives results if CC interactions
    int_part, dtr, etau2 = propagate_nu(e_nu, nu_xc, nu_ixc, d_left) # does the neutrino interact?

    # if int_part != 'lep': # neutrinos at the end;
    #     d_exit = depth
    #     part_type = 'decayed' # (HLS = 0); changed 22/12/2020
    #     e_fin = etau2 # final neutrino energy
    #     # go to 60; all done
    #     # 60 continue
    #     return part_type, d_exit, e_fin
    int_mask = int_part != 'lep'
    part_type[int_mask] = 'decayed'
    d_exit[int_mask] = depth[int_mask]
    e_fin[int_mask] = etau2[int_mask]

    mask = ~int_mask

    # otherwise we have a tau
    # d_lep = d_lep + dtr
    # d_left = d_left - dtr
    d_lep[mask] += dtr[mask]
    d_left[mask] -= dtr[mask]

    # if d_left <= 0:
    #     d_exit = depth
    #     e_fin = etau2
    #     part_type = 'decayed' # went too far to make a tau, so don't count (HLS = 0); changed 22/12/2020
    #     # go to 60; no, still a neutrino
    #     # 60 continue
    #     return part_type, d_exit, e_fin
    dlmask = d_left[mask] <= 0
    full_dlmask = np.copy(mask)
    full_dlmask[full_dlmask] = dlmask
    part_type[full_dlmask] = 'decayed'
    d_exit[full_dlmask] = depth[full_dlmask]
    e_fin[full_dlmask] = etau2[full_dlmask]

    mask[mask] = ~dlmask

    # we have a tau with room to travel for tauthrulayers

    p_t, d_e, e_f = tau_thru_layers(angle, depth[mask], d_water, d_lep[mask], etau2[mask], xc_water, xc_rock, ixc_water, ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, cd2distd, densityatx)

    part_type[mask] = p_t
    d_exit[mask] = d_e
    e_fin[mask] = e_f

    # 60 continue

    return part_type, d_exit, e_fin

# following functions are not useful for the main program and are only designed for debugging!
'''
@njit(nogil=True)
def regen_water(angle, e_lep, depth, d_water, d_lep, nu_xc, nu_ixc, ithird, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong): # note: angle is now in angle

    # find the neutrino energy from the tau decay with tau energy e_lep
    r = my_rand()
    frac = distnu(r,ithird)
    e_nu = frac * e_lep

    d_left = depth-d_lep # this is how far the neutrino can go
    e_fin = e_nu # in case we need to exit
    part_type = 'exit' # in case we need to exit; change later to string (HLS = 2)

    if d_left <= 0: # past the point of interactions allowed
        d_exit = depth
        # go to 60
        # 60 continue
        return part_type, d_exit, e_fin

    d_exit = d_lep # we are starting this far into the Earth with a neutrino
    int_part = 'nu' # starting with a neutrino with energy e_nu; change later to string

    # tnu follows NC to the end, or gives results if CC interactions
    int_part, dtr, etau2 = propagate_nu(e_nu, nu_xc, nu_ixc, d_left) # does the neutrino interact?

    if int_part != 'lep': # neutrinos at the end;
        d_exit = depth
        part_type = 'decayed' # (HLS = 0)
        e_fin = etau2 # final neutrino energy
        # go to 60; all done
        # 60 continue
        return part_type, d_exit, e_fin

    # otherwise we have a tau
    d_lep = d_lep + dtr
    d_left = d_left - dtr

    if d_left <= 0:
        d_exit = depth
        e_fin = etau2
        part_type = 'decayed' # went too far to make a tau, so don't count (HLS = 0)
        # go to 60; no, still a neutrino
        # 60 continue
        return part_type, d_exit, e_fin

    # we have a tau with room to travel

    part_type, d_exit, e_fin = propagate_lep_water(etau2, xc_water, lep_ixc_water, alpha_water, beta_water, d_left, 'stochastic')

    # 60 continue
    return part_type, d_exit, e_fin

@njit(nogil=True)
def regen_rock(angle, e_lep, depth, d_water, d_lep, nu_xc, nu_ixc, ithird, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, cd2distd, densityatx): # note: angle is now in angle

    # find the neutrino energy from the tau decay with tau energy e_lep
    r = my_rand()
    frac = distnu(r,ithird)
    e_nu = frac * e_lep

    d_left = depth-d_lep # this is how far the neutrino can go
    e_fin = e_nu # in case we need to exit
    part_type = 'exit' # in case we need to exit; change later to string (HLS = 2)

    if d_left <= 0: # past the point of interations allowed
        d_exit = depth
        # go to 60
        # 60 continue
        return part_type, d_exit, e_fin

    d_exit = d_lep # we are starting this far into the Earth with a neutrino
    int_part = 'nu' # starting with a neutrino with energy e_nu; change later to string

    # tnu follows NC to the end, or gives results if CC interactions
    int_part, dtr, etau2 = propagate_nu(e_nu, nu_xc, nu_ixc, d_left) # does the neutrino interact?

    if int_part != 'lep': # neutrinos at the end;
        d_exit = depth
        part_type = 'decayed' # (HLS = 0)
        e_fin = etau2 # final neutrino energy
        # go to 60; all done
        # 60 continue
        return part_type, d_exit, e_fin

    # otherwise we have a tau
    d_lep = d_lep + dtr
    d_left = d_left - dtr

    if d_left <= 0:
        d_exit = depth
        e_fin = etau2
        part_type = 'decayed' # went too far to make a tau, so don't count (HLS = 0)
        # go to 60; no, still a neutrino
        # 60 continue
        return part_type, d_exit, e_fin

    # we have a tau with room to travel 
    part_type, d_exit, e_fin = propagate_lep_rock(0, etau2, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, 0.0, d_left, 'stochastic', cd2distd, densityatx)

    # 60 continue
    return part_type, d_exit, e_fin

@njit(nogil=True)
def pexit_w(energy, dm, xc_water, lep_ixc_water, alpha_water, beta_water):

    # dm = dm[0:20]
    imax = int(1)
    # f = open('results_%d' % imax,'w+')

    # propagate_lep_water(e_init, xc_water, lep_ixc, alpha_water_sig, beta_water, d_in, type_loss)

    for j in range(len(dm)):
        surv = 0.0
        ic = 0
        for k in prange(imax):
            p_id, df, e_fin = propagate_lep_water(energy[j], xc_water, lep_ixc_water, alpha_water, beta_water, dm[j], 'stochastic')
            print(energy[j], dm[j], df, e_fin)
            if p_id == 'not_decayed' and e_fin>50:
                # print(e_fin)
                ic+=1
        surv = float(ic)/float(imax)
        # print(energy[j], dm[j], surv, df)
        # f.write(str("%.4f" % energy[j]) + "\t" + str("%.4f" % dm[j]) + "\t" + str("%.4f" % surv) + "\n")
    # f.close()
    return None

@njit(nogil=True)
def pexit_r(angle, energy, dm, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, xalong, cdalong):

    # dm = dm[:5]
    dm = dm[3:5]
    # energy = np.genfromtxt('psurv-ALLM-v9r0-ebin.17', usecols = 0)
    # dm = np.genfromtxt('psurv-ALLM-v9r0-ebin.17', usecols = 1)
    imax = int(1e4)
    # f = open('results_%d_rock' % imax,'w+')

    for j in range(len(dm)):
        surv = 0.0
        ic = 0
        for k in prange(imax):
            # p_id, df, e_fin = propagate_lep_water(energy[j], xc_water, lep_ixc_water, alpha_water, beta_water, dm[j], 'stochastic')
            # propagate_lep_rock(e_init, xc_rock, lep_ixc, alpha_rock, beta_rock, d_entry, d_in, type_loss)
            p_id, df, e_fin = propagate_lep_rock(0, energy[j], xc_rock, lep_ixc_rock, alpha_rock, beta_rock, 0.0, dm[j], 'stochastic', cd2distd, densityatx)
            print(dm[j], e_fin, df)
            if p_id == 'not_decayed' and e_fin>50:
                ic+=1
        # surv = float(ic)/float(imax)
        # print(energy[j], dm[j], surv, df)
        # f.write(str("%.4f" % energy[j]) + "\t" + str("%.4f" % dm[j]) + "\t" + str("%.4f" % surv) + "\n")
    # f.close()
    return None

@njit(nogil=True)
def pexit_n(angle, energy, dm, nu_xc, nu_ixc, xc_water, lep_ixc_water, alpha_water, beta_water, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, xalong, cdalong):
    ebin = np.zeros(91)

    imax = int(1e5)

    for j in range(len(dm)):
        surv = 1.0
        ic = 0
        for k in prange(imax):
            ip, dtravel, ef = propagate_nu(1e10, nu_xc, nu_ixc, dm[j])
            # print(float(1e10), dm[j], dtravel, ip)

            if ip == 'lep':
                etau0 = ef
                dleft = dm[j]-dtravel
                # p_id, df, etauf = propagate_lep_water(etau0, xc_water, lep_ixc_water, alpha_water, beta_water, dleft, 'stochastic')
                p_id, df, etauf = propagate_lep_rock(0, etau0, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, 0.0, dleft, 'stochastic', cd2distd, densityatx)
                if etauf>50 and p_id=='not_decayed':
                    ic+=1
                jef = (np.abs(E_nu-etauf)).argmin()
                if p_id=='not_decayed':
                    ebin[jef] += 1

        surv = float(ic)/float(imax)
        print(float(1e10), dm[j], etauf, df, surv)

    return None

@njit(nogil=True,fastmath=True)
def p_exit_regen(prop_type, angle, energy, dm, nu_xc, nu_ixc, xc_water, lep_ixc_water, alpha_water, beta_water, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, xalong, cdalong):

    ebin = np.zeros(92)
    ebin_regen = np.zeros(92)

    imax = int(1e5)
    # imax = int(1)
    # dm = dm[1:]
    for j in range(len(dm)):
        surv = 1.0
        ic = 0
        ic2 = 0
        # if surv<=0.1:imax=int(1e5)
        # if surv<1e-3:
        #     surv=0.0
        #     continue
        for k in prange(imax):
            # print("propagating neutrinos...")
            ip, dtravel, ef = propagate_nu(1e10, nu_xc, nu_ixc, dm[j])
            # print(float(1e10), dm[j], dtravel, ip)
            etau0 = ef
            dleft = dm[j]-dtravel
            # print("dm = ",dm[j])
            # print("dtravel = ", dtravel)
            # print("dleft = ", dleft)
            # print(float(1e10),dm[j],ef,dtravel,ip)


            if ip == 'lep':

                if prop_type == 'water':
                    # print("propagating taus in water...")
                    p_id, df, etauf = propagate_lep_water(etau0, xc_water, lep_ixc_water, alpha_water, beta_water, dleft, 'stochastic')
                    # print(etau0,dleft,etauf,df,p_id)

                else:
                    # print("propagating taus in rock...")
                    # p_id, df, etauf = propagate_lep_rock(0, etau0, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, 0.0, dleft, xalong, cdalong, 'stochastic')
                    # print("dleft = ", dleft)
                    # print(etau0,dleft)
                    p_id, df, etauf = propagate_lep_rock(angle, etau0, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, 0.0, dleft, 'stochastic', cd2distd, densityatx)

                    # print(etau0,dleft,etauf,df,p_id)


                if etauf>50 and p_id=='not_decayed':
                    jef = (np.abs(E_nu-etauf)).argmin()
                    ic += 1
                    ebin[jef] += 1
                    ic2 += 1
                    ebin_regen[jef] +=1
                # jef = (np.abs(E_nu-etauf)).argmin()

                # if p_id=='not_decayed':
                #     ebin[jef] += 1
                #     ebin_regen[jef] +=1
            else:
                if prop_type == 'water':
                    # print("regenerating taus in water...")
                    # p_id, df, etauf = regen_water(0, etau0, dm[j], 0.0, dtravel, nu_xc, nu_ixc, 0, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong) # ithird = 0
                    p_id, df, etauf = regen(angle, etau0, dm[j], 0.0, dtravel, nu_xc, nu_ixc, 0, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong)
                    # print(etau0,dm[j],0.0,dtravel,etauf,df,p_id)
                else:
                    # print("regenerating taus in rock...")
                    # p_id, df, etauf = regen_rock(0, etau0, dm[j], 0.0, dtravel, nu_xc, nu_ixc, 0, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong) # ithird = 0
                    p_id, df, etauf = regen(angle, etau0, dm[j], 0.0, dtravel, nu_xc, nu_ixc, 0, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong)
                    # print(etau0,dm[j],0.0,dtravel,etauf,df,p_id)

                if etauf>50 and p_id=='not_decayed':
                    jef = (np.abs(E_nu-etauf)).argmin()
                    ic2 +=1
                    ebin_regen[jef] +=1

        surv = float(ic)/float(imax)
        surv2 = float(ic2)/float(imax)
        print(float(1e10), dm[j], surv, surv2)
        # print(float(1e10), dm[j], etauf, df, surv)
        # dmoe = dm[j]*(1e5/etau0)
    for i in range(1,92):
        emid = 10**((float(i+30)-0.5)/10.0)
        # print(dm[j],etau0,emid,ebin[i],ebin_regen[i])
        # print(ebin)
    return None

def pexit_plot(material):
    if material == 'water':
        hls_df = pd.read_table('psurv-ALLM-v9w0-ebin.16', header=None,delimiter=r"\s+")
        hls_df.columns = ["energy","d_in","surv"]
    elif material == 'rock':
        hls_df = pd.read_table('psurv-ALLM-v9r0-ebin.17', header=None,delimiter=r"\s+")
        hls_df.columns = ["energy","d_in","surv"]
    elif material == 'nu_water':
        hls_df = pd.read_table('nupsurv-ALLM-v9w0-ebin.16', header=None,delimiter=r"\s+")
        hls_df.columns = ["energy","d_in","surv"]
    elif material == 'nu_rock':
        hls_df = pd.read_table('nupsurv-ALLM-v9r0-ebin.17', header=None,delimiter=r"\s+")
        hls_df.columns = ["energy","d_in","surv"]
    elif material == 'regen_water':
        hls_df = pd.read_table('nupsurv-ALLM-v9w0-ebin.16', header=None,delimiter=r"\s+")
        hls_df.columns = ["energy","d_in","no_regen","regen_1"]
        # print(hls_df.columns)
    elif material == 'regen_rock':
        hls_df = pd.read_table('nupsurv-ALLM-v9r0-ebin.17', header=None,delimiter=r"\s+")
        hls_df.columns = ["energy","d_in","no_regen","regen_1"]


    # hls_df.columns = ["energy","d_in","surv"]
    # imax = 10000

    if material == 'rock' or material == 'water':

        hls_dict = {1e7:(np.array(hls_df['d_in'][0:26]), np.array(hls_df['surv'][0:26])), 1e8:(np.array(hls_df['d_in'][26:52]),np.array(hls_df['surv'][26:52])), 1e9:(np.array(hls_df['d_in'][52:78]),np.array(hls_df['surv'][52:78])), 1e10:(np.array(hls_df['d_in'][78:104]),np.array(hls_df['surv'][78:104]))}

        my_df = pd.read_table('results_%s' % material, header=None,delimiter=r"\s+")
        my_df.columns = ["energy","d_in","surv"]

        my_dict = {1e7:(np.array(my_df['d_in'][0:26]), np.array(my_df['surv'][0:26])), 1e8:(np.array(my_df['d_in'][26:52]),np.array(my_df['surv'][26:52])), 1e9:(np.array(my_df['d_in'][52:78]),np.array(my_df['surv'][52:78])), 1e10:(np.array(my_df['d_in'][78:]),np.array(my_df['surv'][78:]))}

        plt.figure(1)

        plt.plot(my_dict[1e7][0],my_dict[1e7][1],ls=':',color='r',label=r'$10^{7}$ GeV - NuPyProp')
        plt.plot(hls_dict[1e7][0],hls_dict[1e7][1],color='r',label=r'$10^{7}$ GeV - HLS')

        plt.plot(my_dict[1e8][0],my_dict[1e8][1],ls=':',color='g',label=r'$10^{8}$ GeV - NuPyProp')
        plt.plot(hls_dict[1e8][0],hls_dict[1e8][1],color='g',label=r'$10^{8}$ GeV - HLS')

        plt.plot(my_dict[1e9][0],my_dict[1e9][1],ls=':',color='b',label=r'$10^{9}$ GeV - NuPyProp')
        plt.plot(hls_dict[1e9][0],hls_dict[1e9][1],color='b',label=r'$10^{9}$ GeV - HLS')

        plt.plot(my_dict[1e10][0],my_dict[1e10][1],ls=':',color='k',label=r'$10^{10}$ GeV - NuPyProp')
        plt.plot(hls_dict[1e10][0],hls_dict[1e10][1],color='k',label=r'$10^{10}$ GeV - HLS')

        # plt.xscale('log')

        plt.legend(loc='best')
        plt.xlabel(r"d [km.w.e]")
        plt.ylabel(r"$P_{surv}$")
        # plt.xlim(100, 1e9)
        plt.minorticks_on()
        plt.grid(which='minor', linestyle=':', linewidth='0.2', color='black')

    elif material == 'nu_water' or material == 'nu_rock':
        hls_dict = {1e10:(np.array(hls_df['d_in'][1:]),np.array(hls_df['surv'][1:]))}

        my_df = pd.read_table('results_%s' % material, header=None,delimiter=r"\s+")
        my_df.columns = ["energy","d_in","e_fin","d_fin","surv"]
        my_dict = {1e10:(np.array(my_df['d_in'][1:]),np.array(my_df['surv'][1:]))}

        plt.figure(1)
        plt.plot(my_dict[1e10][0],my_dict[1e10][1],ls=':',color='k',label=r'$10^{10}$ GeV - NuPyProp')
        plt.plot(hls_dict[1e10][0],hls_dict[1e10][1],color='k',label=r'$10^{10}$ GeV - HLS')

        plt.legend(loc='best')
        plt.xlabel(r"d [km.w.e]")
        plt.ylabel(r"$P_{surv}$")

    elif material == 'regen_water' or material == 'regen_rock':
        hls_dict = {1e10:(np.array(hls_df['d_in'][1:]),np.array(hls_df['no_regen'][1:]),np.array(hls_df['regen_1'][1:]))}

        my_df = pd.read_table('results_%s' % material, header=None,delimiter=r"\s+")
        my_df.columns = ["energy","d_in","no_regen","regen_1"]
        my_dict = {1e10:(np.array(my_df['d_in'][1:]),np.array(my_df['no_regen'][1:]),np.array(my_df['regen_1'][1:]))}

        plt.figure(1)
        plt.plot(my_dict[1e10][0],my_dict[1e10][1],ls=':',color='k',label=r'$10^{10}$ GeV - NuPyProp')
        plt.plot(hls_dict[1e10][0],hls_dict[1e10][1],color='k',label=r'$10^{10}$ GeV - HLS')
        plt.title("No Regen - %s" % material[6:])
        plt.legend(loc='best')
        plt.xlabel(r"d [km.w.e]")
        plt.ylabel(r"$P_{surv}$")

        plt.figure(2)
        plt.plot(my_dict[1e10][0],my_dict[1e10][2],ls=':',color='k',label=r'$10^{10}$ GeV - NuPyProp')
        plt.plot(hls_dict[1e10][0],hls_dict[1e10][2],color='k',label=r'$10^{10}$ GeV - HLS')
        plt.title("Regen_1 - %s" % material[6:])
        plt.legend(loc='best')
        plt.xlabel(r"d [km.w.e]")
        plt.ylabel(r"$P_{surv}$")
    return None

def tau_cdf(material):
    if material == 'water':
        hls_df = pd.read_table('nupsurv-ALLM-v9w0-ebin.66', header=None,delimiter=r"\s+")

    elif material == 'rock':
        hls_df = pd.read_table('nupsurv-ALLM-v9r0-ebin.67', header=None,delimiter=r"\s+")


    hls_df.columns = ["d_in","e_in","e_mid","bin_no_regen","bin_regen_1"]

    hls_no_regen = np.cumsum(hls_df['bin_no_regen'],axis=0)
    hls_no_regen = np.asarray(hls_no_regen)/np.asarray(hls_no_regen)[-1] # CDF

    hls_regen_1 = np.cumsum(hls_df['bin_regen_1'],axis=0)
    hls_regen_1 = np.asarray(hls_regen_1)/np.asarray(hls_regen_1)[-1] # CDF


    hls_dict = {1e10:(np.array(hls_df['e_mid']),np.array(hls_df['bin_no_regen']),np.array(hls_df['bin_regen_1']),hls_no_regen,hls_regen_1)}
    # print("HLS max(ebin) = ", max(np.array(hls_df['e_bin'])))

    my_df = pd.read_table('%s_cdf' % material, header=None,delimiter=r"\s+")
    my_df.columns = hls_df.columns

    my_no_regen = np.cumsum(my_df['bin_no_regen'],axis=0)
    my_no_regen = np.asarray(my_no_regen)/np.asarray(my_no_regen)[-1] # CDF

    my_regen_1 = np.cumsum(my_df['bin_regen_1'],axis=0)
    my_regen_1 = np.asarray(my_regen_1)/np.asarray(my_regen_1)[-1] # CDF

    my_dict = {1e10:(np.array(my_df['e_mid']),np.array(my_df['bin_no_regen']),np.array(my_df['bin_regen_1']),my_no_regen,my_regen_1)}
    # print("NPP max(ebin) = ", max(np.array(my_df['e_bin'])))

    plt.figure(1)
    plt.semilogx(my_dict[1e10][0],my_dict[1e10][3],ls=':',color='k',label=r'$10^{10}$ GeV - NuPyProp')
    plt.plot(hls_dict[1e10][0],hls_dict[1e10][3],color='k',label=r'$10^{10}$ GeV - HLS')

    # plt.xlim()

    plt.legend(loc='best')
    plt.xlabel(r"$E_{\tau}$")
    plt.ylabel("CDF")
    plt.title("No Regen - %s" % material)
    plt.show()

    plt.figure(2)
    plt.semilogx(my_dict[1e10][0],my_dict[1e10][4],ls=':',color='k',label=r'$10^{10}$ GeV - NuPyProp')
    plt.plot(hls_dict[1e10][0],hls_dict[1e10][4],color='k',label=r'$10^{10}$ GeV - HLS')

    # plt.xlim()

    plt.legend(loc='best')
    plt.xlabel(r"$E_{\tau}$")
    plt.ylabel("CDF")
    plt.title("Regen - %s" % material)
    plt.show()
    return None

'''
# =============================================================================
# Test
# =============================================================================
if __name__ == "__main__":
    idepth = Geometry.idepth = 4
    start_time = time.time()
    fac_nu = 1
    m_le = Energy_loss.m_tau
    lepton = 'tau'
    material = 'rock'
    cross_section_model = 'ncteq15'
    pn_model = 'allm'
    c_tau = 8.703e-3
    # %timeit -n 10000000

    # print(distnu(r=0.5,ithird=0))
    # x_val = np.linspace(0,1,num=50)
    # distval = [y/3 * (5 - 3* y**2 + y**3) - y/3 * (1 - 3 * y**2 + 2 * y**3) for y in x_val]
    # plt.plot(x_val,distval)
    # plt.show()

    # print(distnu(0,0))


    # rho =
    # print(find_interface())
    nu_xc = Data.get_xc(xc_type='nu',model='ncteq15',particle='neutrino')

    ixc_nu = Data.get_ixc('nu', 'ncteq15', particle='neutrino')
    nu_ixc = ixc_nb(ixc_nu)

    xc_water = Data.get_xc(xc_type=lepton,model=pn_model, material='water')

    alpha_water = Data.get_alpha(particle=lepton, material='water')

    beta_water = Data.get_beta(lepton, 'water', 'continuous', pn_model)

    ixc_water = Data.get_ixc(lepton, model=pn_model, material='water')
    lep_ixc_water = ixc_nb(ixc_water)

    xc_rock = Data.get_xc(xc_type=lepton,model=pn_model, material='rock')

    alpha_rock = Data.get_alpha(particle=lepton, material='rock')

    beta_rock = Data.get_beta(lepton, 'rock', 'continuous', pn_model)

    ixc_rock = Data.get_ixc(lepton, model=pn_model, material='rock')
    lep_ixc_rock = ixc_nb(ixc_rock)

    xalong, cdalong = Data.get_trajs('col', 10, 4)

    energy_nw, dm_nw = np.genfromtxt('nupsurv-ALLM-v9w0-ebin.16', usecols = (0,1), unpack=True)

    energy_nr, dm_nr = np.genfromtxt('nupsurv-ALLM-v9r0-ebin.17', usecols = (0,1), unpack=True)

    # propagate_lep_rock(10, 7706346933.017713, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, 0.0, 40, xalong, cdalong, 'stochastic')

    # pexit_n(0.0, energy_nr, dm_nr, nu_xc, nu_ixc, xc_water, lep_ixc_water, alpha_water, beta_water, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, xalong, cdalong)

    # p_exit_regen('water', 10, energy_nw, dm_nw, nu_xc, nu_ixc, xc_water, lep_ixc_water, alpha_water, beta_water, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, xalong, cdalong)

    # x_inter = Interpolation.cd2distd(xalong, cdalong, 31307.031565357076) # find how far we are along the chord for given beta
    # print("x_inter = ", x_inter)

    # r, rho = Geometry.densityatx(719814756.1985742, 10) # find the density at x
    # print("rho = ", rho)

    # p_exit_regen('rock', 10, energy_nr, dm_nr, nu_xc, nu_ixc, xc_water, lep_ixc_water, alpha_water, beta_water, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, xalong, cdalong)

    # propagate_lep_rock(10, 7706346933.017713, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, 0.0, 11.452173806385758, xalong, cdalong, 'stochastic')


    # p_id, df, etauf = regen_water(0, 1e10, 40, 0.0, 0.0, nu_xc, nu_ixc, 0, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong)
    # print(float(1e10),40,0.0,0.0,etauf,df,p_id)

    # for i in range(100):
    #     print(propagate_nu(1e10, nu_xc, nu_ixc, dm_nw[1]))

    # random.seed(223)

    # p_exit_regen('rock', 0, energy_nr, dm_nr, nu_xc, nu_ixc, xc_water, lep_ixc_water, alpha_water, beta_water, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, xalong, cdalong)

    # pexit_plot('nu_water')
    # pexit_plot('nu_rock')
    # pexit_plot('regen_water')
    # pexit_plot('regen_rock')
    # tau_cdf('water')
    # tau_cdf('rock')
    # energy_water, dm_water = np.genfromtxt('psurv-ALLM-v9w0-ebin.16', usecols = (0,1), unpack=True)

    # energy_rock, dm_rock = np.genfromtxt('psurv-ALLM-v9r0-ebin.17', usecols = (0,1), unpack=True)

    # pexit_w(energy_water, dm_water, xc_water, lep_ixc_water, alpha_water, beta_water)
    # pexit_r(0, energy_rock, dm_rock, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, xalong, cdalong)

    # for i in range(len(dm_rock)):
    #     p_id, df, e_fin = propagate_lep_rock(0, energy_rock[i], xc_rock, lep_ixc_rock, alpha_rock, beta_rock, 0.0, dm_rock[i], xalong, cdalong, 'stochastic')
    #     print(energy_rock[i], dm_rock[i], df, e_fin)

    # for i in range(len(dm_water)):
    #     p_id, df, e_fin = propagate_lep_water(energy_water[i], xc_water, lep_ixc_water, alpha_water, beta_water, dm_water[i], 'stochastic')
    #     print(energy_water[i], dm_water[i], df, e_fin)

    # dm_n = dm_nw
    # for i in range(len(dm_n)):
    #     ip, dtravel, ef = propagate_nu(1e11, nu_xc, nu_ixc, dm_n[i])
    #     # print(float(1e7), dm_n[i], ef, dtravel, ip)
    #     if ip == 'lep':
    #         etau0 = ef
    #         dleft = dm_n[i]-dtravel
    #         p_id, df, etauf = propagate_lep_water(etau0, xc_water, lep_ixc_water, alpha_water, beta_water, dleft, 'stochastic')
    #         print(etau0, dleft, etauf, df, p_id)

    # etau0 = 7381976469.0908289
    # dleft = 360.84061946057818
    # p_id, df, etauf = propagate_lep_water(etau0, xc_water, lep_ixc_water, alpha_water, beta_water, dleft, 'stochastic')
    # p_id, df, etauf = propagate_lep_rock(0, etau0, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, 0.0, dleft, xalong, cdalong, 'stochastic')

    # dm_n = dm_nr
    # for i in range(len(dm_n)):
    #     ip, dtravel, ef = propagate_nu(1e11, nu_xc, nu_ixc, dm_n[i])
    #     # print(float(1e10), dm_n[i], ef, dtravel, ip)
    #     if ip == 'lep':
    #         etau0 = ef
    #         dleft = dm_n[i]-dtravel
    #         p_id, df, etauf = propagate_lep_rock(0, etau0, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, 0.0, dleft, xalong, cdalong, 'stochastic')
    #         print(etau0, dleft, etauf, df, p_id)


    # print(etau0, dleft, etauf, df, p_id)

    end_time = time.time()
    print(f"It took {end_time-start_time:.2f} seconds to compute")
