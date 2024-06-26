#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 11:55:30 2021

@author: sam
    
Updated Wed April 10 14:02:27 2024
    
@author: luke
"""

import nupyprop.data as Data
import nupyprop.constants as const

from collections import OrderedDict
import numpy as np
from scipy import integrate
import sympy as sp
from astropy.table import Table
import warnings
warnings.filterwarnings('ignore')

#Earth Radius in km
Re = const.R_earth

#PREM layers based on R_Earth
Rlay = const.Rlay

#Water density in g/cm^3
rho_water = const.rho_water

#Earth Emergence angles in steps of 0.1 degrees
beta_arr = const.beta_arr

def premdensity(Rin , idepth):
    """
    Calculates the density at a radius inside the Earth according to the PREM model.
    hep-ph/9512364v1 eq. 25
        
    Args:
        Rin (float):
            Radius (in km) at which density is to be calculated.
        idepth (integer):
            Depth (in km) of water layer (sets the last layer).

    Returns:
        float:
            Density (in g/cm^3) at the specified radius.
    """
    #sets the last layer as water layer
    Rlay[9] = 6368.0 + (3.0 - int(idepth))
    y = Rin/Re
    edens = 0.0 
    
    expressions = (
        lambda y: 13.0885-8.8381*y**2,
        lambda y: 12.5815-1.2638*y-3.6426*y**2-5.5281*y**3,
        lambda y: 7.9565-6.4761*y+5.5283*y**2-3.0807*y**3,
        lambda y: 5.3197-1.4836*y,
        lambda y: 11.2494-8.0298*y,
        lambda y: 7.1089-3.8045*y,
        lambda y: 2.6910+0.6924*y,
        lambda y: 2.900,
        lambda y: 2.600,
        lambda y: 1.020,
        )
    
    for i in range(len(Rlay)):
        if Rin <= Rlay[i] or (i == len(Rlay) -1 and Rin <= Rlay[-1]*1.001):
            edens = expressions[i](y)
            break
    
    return edens

def densityatx(x, beta_deg, idepth):
    """Calculates the density at a distance x, for a given Earth emergence angle
        1905.13223v2 fig. 2 for chord length.

    Args:
        x (float): Distance along the chord of the trajectory in km
        beta_deg (float): Earth emergence angle in degrees
        idepth (integer): Depth of water layer in km

    Returns:
        float: Radial distance from the center of the Earth to x in km
        float: Density at x in g/cm^3
    """
    
    #2 R_E sin(beta)
    chord_length = 2*Re*np.sin(np.deg2rad(beta_deg))
    #Just the law of cosines
    r2 = (chord_length-x)**2 + Re**2 - 2*Re*(chord_length-x)*np.sin(np.deg2rad(beta_deg))
    
    if beta_deg < 5.0:
        r = Re*(1.0 + 0.5 *(x**2-chord_length*x)/Re**2)
    else:
        r = np.sqrt(r2)
    
    rho_at_x = premdensity(r,idepth)
    
    return r, rho_at_x

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
    tnadir = (90.0 - beta_deg)*(np.pi/180.0)
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


def PREMgramVSang(beta_deg, idepth):
    chord_length = 2*Re*np.sin(beta_deg*(np.pi/180.))

    col_depth = integrate.quad(lambda x: densityatx(x, beta_deg, idepth)[1], 0, chord_length, epsabs=1.49e-8, epsrel=1.49e-8, maxp1=1000, limit=1000)

    return col_depth[0]*1e5

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
        Column depth, in g/cm^2.

    '''
    if beta_deg<0.5:
        total_col_depth = PREMgramVSang(1.0, idepth)
        columndepth_val = total_col_depth*(beta_deg - (1.0/6.0)*(beta_deg**3))
    else:
        columndepth_val = PREMgramVSang(beta_deg, idepth)
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
        Integrated column depth.

    '''
    r, rho = densityatx(x_v, beta_deg, idepth)
    cdtot_val = rho*1e5
    return cdtot_val

def gen_col_trajs(idepth):
    '''
    Parameters
    ----------
    idepth : int
        Depth of water layer in km.

    Returns'beta':'Earth emergence angle, in degrees',
            'xalong':'Distance in water, in km',
            'cdalong':'Column depth at xalong, in g/cm^2'})
    -------
    betad_fix : ndarray
        1D array containing 13500 entries, with repeating (x100) entries from 0.1 deg to 90 deg.
    xalong : ndarray
        1D array containing distance in water, in km.
    cdalong : ndarray
        1D array containing column depth at xalong, in g/cm^2.

    Essentially, cdalong = integral(cdtot(x_v, beta, idepth)), where x_v limits go from 0 to cdalong, for each beta value.

    '''

    betad_fix = np.repeat(beta_arr,100) # used for storing in hdf5 and matching len of cdalong array which uses list comprehension for faster calcs.
    chord = np.asarray([trajlength(float(i)) for i in beta_arr])
    dx = chord/100.0 # step size
    xalong = np.asarray([i*float(j) for i in dx for j in range(1,101)])
    cdalong = np.asarray([integrate.quad(cdtot,0,j,args=(i,idepth))[0] for i,j in zip(betad_fix,xalong)])
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

    beta_col, xalong, cdalong = gen_col_trajs(idepth)
    col_meta = OrderedDict({'Description':'Column trajectories for water layer = %s km' % str(idepth),
                            'beta':'Earth emergence angle, in degrees',
                            'xalong':'Distance in water, in km',
                            'cdalong':'Column depth at xalong, in g/cm^2'})
    col_table = Table([beta_col, xalong, cdalong], names=('beta','xalong','cdalong'), meta=col_meta)

    beta_water, chord, water = gen_water_trajs(idepth)
    water_meta = OrderedDict({'Description':'Water trajectories for water layer = %s km' % str(idepth),
                              'beta':'Earth emergence angle, in degrees',
                              'chord':'Chord length, in km',
                              'water':'Final water layer distance, in km'})
    water_table = Table([beta_water, chord, water], names=('beta','chord','water'), meta=water_meta)

    Data.add_trajs('col', int(idepth), col_table)
    Data.add_trajs('water', int(idepth), water_table) # yikes! fixed!
    return None
