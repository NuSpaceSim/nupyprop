#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 17:54:37 2020

@author: sam
"""
import data as Data
# import cross_section as Cross_section
# import energy_loss as Energy_loss

import numpy as np
import pandas as pd
import time
# from numba import njit
import scipy.constants as scc
from interpolation import interp

import warnings
warnings.filterwarnings('ignore')

pd.set_option("display.max_rows", None, "display.max_columns", None)


E_nu = Data.E_nu
E_lep = Data.E_lep
N_A = scc.Avogadro


# @njit(nogil=True)
def interpol(x_vals, x, y):
    return interp(x, y, x_vals)

# @njit(nogil=True)
def cd2distd(xalong, cdalong, col_depth):

    if col_depth < np.min(cdalong):
        return (col_depth/cdalong[0])*xalong[0]
    elif col_depth > np.max(cdalong):
        return np.max(xalong)
    else:
        return interpol(col_depth, cdalong, xalong)

def l_cd2distd(col_depth, xalong, cdalong, mn, mx, xx):
    d = np.empty_like(col_depth)
    d[col_depth < mn] = (col_depth[col_depth < mn] / cdalong[0]) * xalong[0]
    d[col_depth > mx] = xx
    mask = ~((col_depth < mn) | (col_depth > mx))
    d[mask] = interpol(col_depth[mask], cdalong, xalong)
    return d
    # # if col_depth < np.min(cdalong):
    # #     return (col_depth/cdalong[0])*xalong[0]
    # elif col_depth > np.max(cdalong):
    #     return np.max(xalong)
    # else:
    #     return interpol(col_depth, cdalong, xalong)

def f_cd2distd(xalong, cdalong):
    mn = np.min(cdalong)
    mx = np.max(cdalong)
    xx = np.max(xalong)
    return lambda col_depth: l_cd2distd(col_depth, xalong, cdalong, mn, mx, xx)

# @njit(nogil=True)
def int_xc_nu(energy, nu_xc): # interpolate xc's & multiply them by N_A (only for neutrinos)
    sig_cc = N_A*interpol(energy, E_nu, nu_xc[0]) # has to be multiplied by N_A because nu_int_length depends on N_A
    sig_nc = N_A*interpol(energy, E_nu, nu_xc[1]) # has to be multiplied by N_A because nu_int_length depends on N_A
    return sig_cc, sig_nc

# @njit(nogil=True)
def int_xc_lep(energy, xc_arr): # interpolate xc's

        sig_brem = interpol(energy, E_lep, xc_arr[0]) # doesn't need to be multiplied by N_A since they already are in the energy_loss module
        sig_pair = interpol(energy, E_lep, xc_arr[1]) # doesn't need to be multiplied by N_A since they already are in the energy_loss module
        sig_pn = interpol(energy, E_lep, xc_arr[2]) # doesn't need to be multiplied by N_A since they already are in the energy_loss module
        return sig_brem, sig_pair, sig_pn


# @njit(nogil=True)
def int_alpha(energy, alpha_sig):
    alpha = interpol(energy, E_lep, alpha_sig)
    return alpha

# @njit(nogil=True)
def int_beta(energy, beta_arr):
    brem = interpol(energy, E_lep, beta_arr[0])
    pair = interpol(energy, E_lep, beta_arr[1])
    pn = interpol(energy, E_lep, beta_arr[2])
    return brem+pair+pn

    # return "Error in int_beta function in interpolation"

# def int_ixc(energy, ixc_dict, ip): # interpolate ixc's; too slow so don't use for now - see find_y in transport
#     try:
#         df = pd.DataFrame.from_dict(ixc_dict[ip]) # ip=cc;nc;brem;pair;pn
#         df.insert(0,energy,np.nan)
#         df = df.sort_index(axis=1) # put the energy column where it belongs
#         cols = np.array(df.columns)
#         df_interp = df[cols].apply(lambda x: x.interpolate(method='index'), axis=1)
#         # df_interp =
#         return df_interp.to_dict()[float(energy)] # needs to be a dictionary for faster computing
#     except ValueError:
#         return ixc_dict[ip][float(energy)]
#     return "Error in int_ixc function in interpolation"

# def int_ixc_2(energy, ixc_dict, ip): # too slow so don't use for now - see find_y in transport
#     try:
#         df = pd.DataFrame.from_dict(ixc_dict[ip]) # ip=cc;nc;brem;pair;pn
#         df.insert(0,energy,np.nan)
#         df = df.sort_index(axis=1) # put the energy column where it belongs
#         cols = np.array(df.columns)
#         df_interp = df[cols].apply(lambda x: x.interpolate(method='index'), axis=1)
#         # df_interp =
#         return df_interp.to_dict()[energy] # needs to be a dictionary for faster computing
#     except ValueError:
#         return ixc_dict[ip][energy]
#     return None

# =============================================================================
# Test
# =============================================================================
if __name__ == "__main__":
    start_time = time.time()
    lepton = 'tau'
    material = 'rock'
    pn_model = 'allm'
    # print(cd2distd(xalong, cdalong, depth))
    # nu_xc = Data.get_xc(xc_type='nu',model='ncteq15',particle='neutrino')


    # lep_xc_water = Data.get_xc(xc_type='tau',model='allm', material='water')


    # alpha_water = Data.get_alpha(particle=lepton, material='water')


    beta_water = Data.get_beta(lepton, 'water', 'continuous', pn_model)
    print(int_beta(1e7, beta_water, log=False))



    # print(int_xc(123456, nu_xc, 'nu'))
    # for i in range(100):
    #     # xalong, cdalong = Data.get_trajs('col', 1.0, 4)
    #     energy = random.randrange(1e3,1e12)
    #     # int_xc_lep(energy, xc_water_brem, xc_water_pair, xc_water_pn)
    #     int_alpha(energy, alpha_water)
    #     int_beta(energy, beta_water)

    # xalong, cdalong = Data.get_trajs('col', 10, 4)
    # print(cd2distd(xalong, cdalong, 1e2))

    # x_inter = cd2distd(xalong, cdalong, 31307.031565357076) # find how far we are along the chord for given beta
    # print(x_inter)

    # print(cd2distd(xalong, cdalong, 1e3))
    # cd2 = jit(nopython=True)(cd2distd)
    # jit('float64(float64[:])', nopython=True)(cd2distd(xalong, cdalong, 1e3))
    # print(cd2distd(xalong, cdalong, 226827))
    # print(cd2distd(xalong, cdalong, min(cdalong)-1))
    # print(cd2distd(xalong, cdalong, max(cdalong)+1))
    end_time = time.time()
    print(f"It took {end_time-start_time:.2f} seconds to compute")
