#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 11:55:30 2021

@author: sam
"""

import nupyprop.data as Data

import numpy as np
from scipy import integrate
import multiprocessing as mp
from multiprocessing import Pool
import pandas as pd
import sympy as sp
import warnings
pd.set_option("display.max_rows", None, "display.max_columns", None)
warnings.filterwarnings('ignore')

# from propagate import geometry as Geometry
from .propagate import geometry as Geometry

import importlib_resources

ref = importlib_resources.files('nupyprop.datafiles') / 'lookup_tables.h5'

Re = 6371.0 # radius of the earth in km
Rlay = np.array([1221.5, 3480.0, 5701.0, 5771.0, 5971.0, 6151.0, 6346.6, 6356.0, 6368.0, 6371.0]) # PREM layers based on R_earth. If you're using another Earth model, be sure to change it here as well as propagate.f90, in PREMdensity subroutine.

rho_water = 1.02 # density of water in g/cm^3
beta_arr = np.asarray([float('{:.1f}'.format(i)) for i in np.concatenate((np.linspace(0.1,5,50), np.linspace(6,90,85)))])

def sagitta_deg(beta_deg):
    '''

    Parameters
    ----------
    beta_deg : float
        Earth emergence angle (beta), in degrees.

    Returns
    -------
    sagitta : float
        The sagitta (distance from the center of the arc to the center of its base), in km.

    '''
    tnadir = (90.0-beta_deg)*(np.pi/180.0)
    # sagitta = sagitta(tnadir)
    sagitta = Re*(1.0 - np.sin(tnadir))
    return sagitta

def pathlength(tnadir):
    '''

    Parameters
    ----------
    tnadir : float
        Nadir angle, in radians.

    Returns
    -------
    pathlength : float
        Path length, in km.

    '''
    pathlength = 2*Re*np.cos(tnadir)
    return pathlength

def trajlength(beta_deg):
    '''

    Parameters
    ----------
    beta_deg : float
        Earth emergence angle (beta), in degrees.

    Returns
    -------
    trajlength : float
        Trajectory length, in km.

    '''
    tnadir = (90.0-beta_deg)*(np.pi/180.0)
    traj_length = Re*np.cos(tnadir)*2
    return float(traj_length)

def PREMgramVSang(z, idepth):
    '''

    Parameters
    ----------
    z : float
        Angle of trajectory (relative to tangent to surface), in degrees.
    idepth : int
        Depth of water layer in km.

    Returns
    -------
    gramlen : float
        The "grammage" or column density, in g/cm^2.

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
            rho = Geometry.premdensity(Rin, idepth)
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
            rhoOUT = Geometry.premdensity(Rin, idepth)
            glen[i] = clen[i]*1e5*rhoOUT

        gramlen += glen[i]

    return gramlen

def columndepth(beta_deg, idepth):
    '''

    Parameters
    ----------
    beta_deg : float
        Earth emergence angle (beta), in degrees.
    idepth : int
        Depth of water layer in km.

    Returns
    -------
    columndepth_val : float
        Column density, in g/cm^2.

    '''
    z = beta_deg
    if z<0.5:
        z1 = 1.0
        c = PREMgramVSang(z1, idepth)
        columndepth_val = c*(z/z1-1.0/6.0*(z/z1)**3)
    else:
        columndepth_val = PREMgramVSang(z, idepth)
    return columndepth_val

def cdtot(x_v, beta_deg, idepth):
    '''

    Parameters
    ----------
    x_v : None
        Lambda function (integration variable).
    beta_deg : float
        Earth emergence angle (beta), in degrees.
    idepth : int
        Depth of water layer in km.

    Returns
    -------
    cdtot_val : float
        Integrated column density.

    '''
    r, rho = Geometry.densityatx(x_v, beta_deg, idepth)
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
    '''

    Returns
    -------
    betad_fix : ndarray
        1D array containing 13500 entries, with repeating entries from 0.1 deg to 90 deg.
    xalong : ndarray
        1D array containing distance in water, in km.
    cdalong : ndarray
        1D array containing column depth at xalong, in g/cm^2.

    '''
    p = Pool(mp.cpu_count()) # use all available cores

    betad_fix = np.asarray([np.ones(100)*i for i in beta_arr]).flatten() # used for storing in hdf5 and matching len of cdalong array which uses list comprehension for faster calcs.
    chord = np.asarray([trajlength(float(i)) for i in beta_arr])
    dx = np.asarray([i/100.0 for i in chord])
    xalong = np.asarray([i*float(j) for i in dx for j in range(1,101)])
    cdalong = np.asarray([p.map(integrator,[[cdtot, i, 0, j]])[0][0] for i,j in zip(betad_fix,xalong)])
    p.close()
    return betad_fix, xalong, cdalong

def gen_water_trajs(idepth):
    '''


    Parameters
    ----------
    idepth : int
        Depth of water layer in km.

    Returns
    -------
    beta_arr : ndarray
        1D array containing Earth emergence angles from 0.1 deg to 90 deg.
    chord : ndarray
        1D array containing chord length, in km.
    water : ndarray
        1D array containing final water layer distance, in km.

    '''
    dw = idepth
    Rrock = Re-dw
    chord = np.asarray([trajlength(i) for i in beta_arr])
    water = np.asarray([0.5*(i - np.sqrt(i**2 - 4.0 * (dw**2 + 2.0*Rrock*dw))) if sagitta_deg(j)>dw else i for i,j in zip(chord, beta_arr)])
    return beta_arr, chord, water

def find_interface(idepth):
    '''


    Parameters
    ----------
    idepth : int
        Depth of water layer in km.

    Returns
    -------
    float
        The (Earth emergence) angle at which the transition from water to rock occurs, in degrees.

    '''
    # global idepth
    y = sp.Symbol('y')
    eqn = sp.Eq(Re*(1-sp.sin((90-y)*(np.pi/180))), idepth)
    return sp.solve(eqn)

def create_traj_table(idepth):
    '''

    Parameters
    ----------
    idepth : int
        Depth of water layer in km.

    Returns
    -------
    None
        Adds trajectory table lookup entries in lookup_table.h5.

    '''
    try:
        with importlib_resources.as_file(ref) as lookup_tables:
            pd.read_hdf(lookup_tables,'Earth/traj_%s/Column_Trajectories' % str(int(idepth)))[0:2]
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
    print(Geometry.premdensity(1203, idepth))
    # print(premdensity_2(1203))

    # premdensity(1203, idepth)
    # r, rho = densityatx(6752.231264859498,32.0, idepth)
    # print(r,rho)
    # r,rho = densityatx(719814756.1985742, 10, idepth)
    # print(r,rho)
    # print(columndepth(10, idepth))
    # thn = 1.3962622222222221
    # print(columndepth(thn, idepth))
    # print(PREMgramVSang(10, idepth))
