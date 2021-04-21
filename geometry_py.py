#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 8 11:58:16 2020

@author: sam
"""
# import h5py
import data as Data

import numpy as np
# cimport numpy as np
from scipy import integrate
import multiprocessing as mp
from multiprocessing import Pool
# from tables import *
import pandas as pd
import sympy as sp
import warnings
# from numba import njit


pd.set_option("display.max_rows", None, "display.max_columns", None)
warnings.filterwarnings('ignore')


Re = 6371.0 # radius of the earth in km
Rlay = np.array([1221.5, 3480.0, 5701.0, 5771.0, 5971.0, 6151.0, 6346.6, 6356.0, 6368.0, 6371.0]) # PREM layers based on R_earth
rho_water = 1.02 # density of water in g/cm^3
beta_arr = np.asarray([float('{:.1f}'.format(i)) for i in np.concatenate((np.linspace(0.1,5,50), np.linspace(6,90,85)))])

# Rlay_mod = np.copy(Rlay)
# Rlay_mod[8] = 6368.0+(3.0-float(idepth)) # fixed Rlay, accounting for depth of water/ice.

# prem_density_functions = [
#         lambda y : 13.0885-8.8381*y**2,
#         lambda y : 12.5815-1.2638*y-3.6426*y**2-5.5281*y**3,
#         lambda y : 7.9565-6.4761*y+5.5283*y**2-3.0807*y**3,
#         lambda y : 5.3197-1.4836*y,
#         lambda y : 11.2494-8.0298*y,
#         lambda y : 7.1089-3.8045*y,
#         lambda y : 2.6910+0.6924*y,
#         lambda y : 2.900,
#         lambda y : 2.600,
#         lambda y : 1.020,
#         lambda y : 1.020,
#         lambda y : 0.0
#         ]

# @njit(nogil=True)
def sagitta_deg(beta_deg):
    '''

    Parameters
    ----------
    beta_deg : float
        Beta in degrees

    Returns
    -------
    sagitta : float
        The sagitta in km.

    '''
    tnadir = (90.0-beta_deg)*(np.pi/180.0)
    # sagitta = sagitta(tnadir)
    sagitta = Re*(1.0 - np.sin(tnadir))
    return sagitta

# @njit(nogil=True)
def pathlength(tnadir):
    '''

    Parameters
    ----------
    tnadir : float
        tnadir in rad.

    Returns
    -------
    pathlength : float
        Path length in km.

    '''
    pathlength = 2*Re*np.cos(tnadir)
    return pathlength

# @njit(nogil=True)
def trajlength(beta_deg):
    '''

    Parameters
    ----------
    beta_deg : float
        Beta in degrees.

    Returns
    -------
    trajlength : float
        Trajectory length in km.

    '''
    tnadir = (90.0-beta_deg)*(np.pi/180.0)
    traj_length = Re*np.cos(tnadir)*2
    return float(traj_length)



# @njit(nogil=True)
def PREMdensity(Rin, idepth):
    '''

    Parameters
    ----------
    Rin : float
        Distance from the Earth Center (at sagitta), in km.

    Returns
    -------
    rhoOUT : float
        Density in g/cm^3.

    '''
    # print("idepth = ", idepth)
    Rlay_2 = np.copy(Rlay)
    Rlay_2[8] = 6368.0+(3.0-float(idepth))
    # print("Rlay_2[8] = ", Rlay_2[8])

    x=Rin
    y=x/Re

    # print('x=', x)

    if (x<=Rlay_2[0]):
        # edens=13.0885-8.8381*y**2
        return 13.0885-8.8381*y**2
    elif (x<=Rlay_2[1]):
        # edens=12.5815-1.2638*y-3.6426*y**2-5.5281*y**3
        return 12.5815-1.2638*y-3.6426*y**2-5.5281*y**3
    elif (x<=Rlay_2[2]):
        # edens=7.9565-6.4761*y+5.5283*y**2-3.0807*y**3
        return 7.9565-6.4761*y+5.5283*y**2-3.0807*y**3
    elif (x<=Rlay_2[3]):
        # edens=5.3197-1.4836*y
        return 5.3197-1.4836*y
    elif (x<=Rlay_2[4]):
        # edens=11.2494-8.0298*y
        return 11.2494-8.0298*y
    elif (x<=Rlay_2[5]):
        # edens=7.1089-3.8045*y
        return 7.1089-3.8045*y
    elif (x<=Rlay_2[6]):
        # edens=2.6910+0.6924*y
        return 2.6910+0.6924*y
    elif (x<=Rlay_2[7]):
        # edens=2.900
        return 2.900
    elif (x<=Rlay_2[8]):
        # edens=2.600
        return 2.600
    elif (x<=Rlay_2[9]):
        # edens=1.020
        return 1.020
    elif (x<=Rlay_2[9]*1.001): # too close to call!
        return 1.020
    else:
        # edens=0.
        return 0.

def PREMgramVSang(z, idepth):
    '''

    Parameters
    ----------
    z : float or integer
        Angle in degrees of trajectory (relative to tangent to surface).

    Returns
    -------
    gramlen : float
        The "grammage", column density in g/cm2.

    '''
    Rlay_3 = np.copy(Rlay)
    Rlay_3[8] = 6368.0+(3.0-float(idepth))
    y = z*(np.pi/180)
    Chord = 2.0*Re*np.sin(y)
    Depth = Re-(0.5*np.sqrt(4.0*Re**2-Chord**2))
    Rdep = Re-Depth

    Rlen = np.zeros(10)
    RlenS = np.zeros(10)
    clen = np.zeros(10)
    glen = np.zeros(10)

    ilay = np.asarray([1 if Rdep<Rlay_3[i] else 0 for i in range(len(Rlay_3))])
    gramlen = 0

    first_non_zero = np.nonzero(ilay)[0][0]

    for i in range(first_non_zero, len(ilay)):
        if i == first_non_zero:
            ifirst = i
            Rlen[i] = Rlay_3[i]-Rdep
            clen[i] = 2*np.sqrt(Rlay_3[i]**2 - (Rlay_3[i] - Rlen[i])**2)
            Rin = Rlay_3[i] - (Rlen[i]/2.0)
            rho = PREMdensity(Rin, idepth)
            glen[i] = clen[i]*1e5*rho

        elif ilay[i] > 0:
            Rlen[i] = Rlay_3[i] - Rlay_3[i-1]
            RlenS[i] = Rlay_3[i] - Rlay_3[i-1]

            for j in range(ifirst,i):
                RlenS[i] = RlenS[i] + Rlen[j]

            clen[i] = 2.0*np.sqrt(Rlay_3[i]**2 - (Rlay_3[i] - RlenS[i])**2)
            # print(clen[i])

            for j in range(ifirst,i):
                clen[i] = clen[i] - clen[j]

            Rin = Rlay_3[i] - (Rlen[i]/2.0)
            rhoOUT = PREMdensity(Rin, idepth)
            glen[i] = clen[i]*1e5*rhoOUT

        gramlen += glen[i]

    return gramlen

def columndepth(beta_deg, idepth):
    '''

    Parameters
    ----------
    tnadir : float
        tnadir in rad.

    Returns
    -------
    columndepth : float
        Column density in g/cm^2.

    '''
    z = beta_deg
    if z<0.5:
        z1 = 1.0
        c = PREMgramVSang(z1, idepth)
        columndepth_val = c*(z/z1-1.0/6.0*(z/z1)**3)
    else:
        columndepth_val = PREMgramVSang(z, idepth)
    return columndepth_val

# @njit(nogil=True)
def densityatx(x, beta_deg, idepth):
    '''

    Parameters
    ----------
    x : float
        Position to find the density at, in km.
    beta_deg : float
        beta in degrees.

    Returns
    -------
    rho_at_x : float
        Density at the position x, in g/cm^3.

    '''
    tnadir = (90.0-beta_deg)*(np.pi/180.0)
    ell = Re*np.cos(tnadir)*2
    r2 = x**2 - (ell*x) + Re**2
    if beta_deg<5.0:
        r = Re*(1.0 + 0.5*(x**2-ell*x)/Re**2)
    else:
        r = np.sqrt(r2)
    rho_at_x = PREMdensity(r, idepth)

    print('r=',r)

    return r, rho_at_x

def cdtot(x_v, beta_deg, idepth):
    '''

    Parameters
    ----------
    x_v : lambda function
        Integration variable.
    beta_deg : float
        beta in degrees.

    Returns
    -------
    cdtot : float
        Integrated column density

    '''
    r, rho = densityatx(x_v, beta_deg, idepth)
    cdtot_val = rho*1e5
    return cdtot_val

def integrator(args):
    '''

    Parameters
    ----------
    args : list of lists
        Function arguments (used for multithreading later).

    Returns
    -------
    function
        Integrator function.

    '''
    fun = args[0]
    var = args[1]
    low_lim = args[2]
    up_lim = args[3]
    return integrate.quad(fun,low_lim,up_lim,args=(var))

def gen_col_trajs():

    p = Pool(mp.cpu_count()) # use all available cores

    betad_fix = np.asarray([np.ones(100)*i for i in beta_arr]).flatten() # used for storing in hdf5 and matching len of cdalong array which uses list comprehension for faster calcs.
    chord = np.asarray([trajlength(float(i)) for i in beta_arr])
    dx = np.asarray([i/100.0 for i in chord])
    xalong = np.asarray([i*float(j) for i in dx for j in range(1,101)])
    cdalong = np.asarray([p.map(integrator,[[cdtot, i, 0, j]])[0][0] for i,j in zip(betad_fix,xalong)])
    p.close()
    return betad_fix, xalong, cdalong

def gen_water_trajs(idepth):
    dw = idepth
    Rrock = Re-dw
    chord = np.asarray([trajlength(i) for i in beta_arr])
    water = np.asarray([0.5*(i - np.sqrt(i**2 - 4.0 * (dw**2 + 2.0*Rrock*dw))) if sagitta_deg(j)>dw else i for i,j in zip(chord, beta_arr)])
    return beta_arr, chord, water

def find_interface(idepth):
    # global idepth
    y = sp.Symbol('y')
    eqn = sp.Eq(Re*(1-sp.sin((90-y)*(np.pi/180))), idepth)
    return sp.solve(eqn)

def create_traj_table(idepth):
    try:
        pd.read_hdf('lookup_tables.h5','Earth/traj_%s/Column_Trajectories' % str(int(idepth)))[0:2]
        return print("idepth = %s already exists in the lookup table. Will initialize that data." % (int(idepth)))
    except (KeyError, FileNotFoundError) as e:
        beta_col, xalong, cdalong = gen_col_trajs()
        dataset_col = pd.DataFrame({'beta':beta_col, 'xalong':xalong, 'cdalong':cdalong}).sort_values(by=['beta','xalong','cdalong'])

        beta_water, chord, water = gen_water_trajs(idepth)
        dataset_water = pd.DataFrame({'beta':beta_water, 'chord':chord, 'water':water})

        Data.add_trajs('col', int(idepth), dataset_col)
        Data.add_trajs('water', int(idepth), dataset_water) # yikes! fixed!
        return None


# =============================================================================
# Test
# =============================================================================
if __name__ == "__main__":
    # import data as Data

    idepth = 6
    # print(find_interface(idepth)[0])
    print(PREMdensity(1203, idepth))
    # print(PREMdensity_2(1203))

    # PREMdensity(1203, idepth)
    # r, rho = densityatx(6752.231264859498,32.0, idepth)
    # print(r,rho)
    # r,rho = densityatx(719814756.1985742, 10, idepth)
    # print(r,rho)
    # print(columndepth(10, idepth))
    # thn = 1.3962622222222221
    # print(columndepth(thn, idepth))
    # print(PREMgramVSang(10, idepth))