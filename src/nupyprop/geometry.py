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
from scipy.integrate import quad
import sympy as sp
from astropy.table import Table
import warnings
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.integrate import cumulative_trapezoid

warnings.filterwarnings('ignore')

#Earth Radius in km
Re = const.R_earth

#PREM layers based on R_Earth
Rlay = const.Rlay

#Water density in g/cm^3
rho_water = const.rho_water

#Earth Emergence angles in steps of 0.1 degrees
beta_arr = const.beta_arr

ak135_density_data = np.genfromtxt("ak135_density.txt", skip_header=1)
depth_ak135 = ak135_density_data[:,0]
density_ak135 = ak135_density_data[:,1]

def ak135density(Rin_array, idepth):
    """
    ak135 Earth density for multiple Earth radii at once.
    https://ses.cidp.edu.cn/Geophys-J-Int-1995-Kennett-108-24.pdf (Page 122)
    
    Args:
        Rin_array (numpy.ndarray): Array of radii (km) where density is to be computed.
        idepth (int): Depth of water layer (km).

    Returns:
        numpy.ndarray: Densities at the given radii (g/cm^3).
    """
    # Create local copies to avoid modifying global variables
    depth_local = depth_ak135
    density_local = density_ak135

    depth1_array = Re - Rin_array  # Convert radius to depth
    
    if idepth == 0:
        depth_local = depth_local[4:]
        density_local = density_local[4:]
    
    indices = np.searchsorted(depth_local, depth1_array) 

    # Assign values based on depth
    density_result = density_local[indices]

    # Apply idepth condition efficiently
    density_result[depth1_array <= idepth] = density_local[0]

    return density_result

def premdensity(Rin_array, idepth):
    """
    PREM Earth density model. Computes density at multiple radii of Earth simultaneously. hep-ph/9512364v1 eq. 25

    Args:
        Rin_array (float or ndarray): 
            Radius (in km) at which density is to be calculated.
        idepth (integer): 
            Depth (in km) of water layer (sets the last layer).

    Returns:
        ndarray: Density (in g/cm^3) at the specified radius (same shape as `Rin_array`).
    """
    # Ensure Rin_array is a NumPy array for vectorized operations
    Rin_array = np.asarray(Rin_array)

    # Update last layer for water depth
    Rlay[9] = 6368.0 + (3.0 - int(idepth))
    y = Rin_array / Re  # Normalize by Earth radius

    # Define density expressions as NumPy functions
    densities = np.array([
        13.0885 - 8.8381 * y**2,
        12.5815 - 1.2638 * y - 3.6426 * y**2 - 5.5281 * y**3,
        7.9565 - 6.4761 * y + 5.5283 * y**2 - 3.0807 * y**3,
        5.3197 - 1.4836 * y,
        11.2494 - 8.0298 * y,
        7.1089 - 3.8045 * y,
        2.6910 + 0.6924 * y,
        np.full_like(y, 2.900),
        np.full_like(y, 2.600),
        np.full_like(y, 1.020 if idepth > 0 else 2.600)
    ])

    # Use np.select() to apply the correct density function based on conditions
    conditions = [(Rin_array <= Rlay[i]) | ((i == len(Rlay) - 1) & (Rin_array <= Rlay[-1] * 1.001)) for i in range(len(Rlay))]
    edens = np.select(conditions, densities, default=0.0)  # Default to 0 if no condition matches

    return edens


def densityatx(x_array, beta_array, idepth, model_name):
    """Calculates the density at a distance x, for a given Earth emergence angle
        1905.13223v2 fig. 2 for chord length. Takes multiple x and beta values. 

    Args:
        x_array (numpy.ndarray): Distances along the trajectory in km.
        beta_array (numpy.ndarray): Corresponding Earth emergence angles in degrees.
        idepth (int): Depth of water layer in km.
        model_name (str): "PREM" or "ak135".

    Returns:
        tuple: (array of radii, array of densities).
    """
    x_array = np.atleast_1d(x_array)
    beta_array = np.atleast_1d(beta_array)  # Force array conversion
    
    # Compute trajectory length using vectorized operations
    chord_length = 2 * Re * np.sin(np.deg2rad(beta_array))

    # Compute radial distance using the law of cosines (vectorized)
    r2 = (chord_length - x_array) ** 2 + Re ** 2 - 2 * Re * (chord_length - x_array) * np.sin(np.deg2rad(beta_array))
    
    r = np.sqrt(r2)

    # Small-angle approximation for low beta values (avoid precision issues)
    small_angle_mask = beta_array < 5.0
    r[small_angle_mask] = Re * (1.0 + 0.5 * (x_array[small_angle_mask]**2 - chord_length[small_angle_mask] * x_array[small_angle_mask]) / Re**2)

    # Select density model and compute density
    if model_name == "PREM":
        rho_at_x = premdensity(r, idepth)  # Assuming premdensity is also vectorized
    elif model_name == "ak135":
        rho_at_x = ak135density(r, idepth)
    else:
        raise ValueError("Input model_name must be 'PREM' or 'ak135'.")

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


def PREMgramVSang(beta_deg, idepth, model_name):
    ''' 
    Computes column depth for given emergence angle using vectorized integration.

    Args:
        beta_deg (float or ndarray): Earth emergence angle(s) in degrees.
        idepth (int): Depth of water layer in km.
        model_name (str): Earth density model name (PREM or ak135).

    Returns:
        float or ndarray: Column depth in g/cm^2.
    '''
    beta_rad = np.radians(beta_deg)

    # Compute chord length (vectorized)
    chord_length = 2 * Re * np.sin(beta_rad)

    # Define the integrand function (ensuring x is always an array)
    def integrand(x, beta):
        x_array = np.atleast_1d(x)  # Convert scalar x to an array
        _, rho = densityatx(x_array, beta, idepth, model_name)
        return rho if np.isscalar(x) else rho[0]  # Return scalar if x is scalar

    # Perform integration for each beta value
    col_depths = np.array([
        integrate.quad(integrand, 0, cl, args=(b,), epsabs=1.49e-8, epsrel=1.49e-8)[0]
        for cl, b in zip(chord_length, beta_deg)
    ])

    return col_depths * 1e5 if col_depths.shape else col_depths[0] * 1e5 

def columndepth(beta_deg, idepth):
    '''

    Parameters
    ----------
    beta_deg : float
        Earth emergence angle (beta), in degrees.
    idepth : int
        Depth of water layer in km.
    model_name : string 
        Earth density model name (PREM or ak135)

    Returns
    -------
    columndepth_val : float
        Column depth, in g/cm^2.

    '''
    if beta_deg<0.5:
        total_col_depth = PREMgramVSang(1.0, idepth, model_name)
        columndepth_val = total_col_depth*(beta_deg - (1.0/6.0)*(beta_deg**3))
    else:
        columndepth_val = PREMgramVSang(beta_deg, idepth, model_name)
    return columndepth_val

def cdtot(x_v, beta_deg, idepth, model_name):
    '''

    Parameters
    ----------
    x_v : None
        Lambda function (integration variable).
    beta_deg : float
        Earth emergence angle (beta), in degrees.
    idepth : int
        Depth of water layer in km.
    model_name : string 
        Earth density model name (PREM or ak135)

    Returns
    -------
    cdtot_val : float
        Integrated column depth.

    '''
    r, rho = densityatx(x_v, beta_deg, idepth, model_name)
    cdtot_val = rho*1e5
    return cdtot_val

def gen_col_trajs(idepth, model_name):
    '''
    Parameters
    ----------
    idepth : int
        Depth of water layer in km.
    model_name : string
       Earth density model name (PREM or ak135)

    Returns
    -------
    betad_fix : ndarray
        1D array containing 13500 entries, with repeating (x100) entries from 0.1 deg to 90 deg.
    xalong : ndarray
        1D array containing distance in water, in km.
    cdalong : ndarray
        1D array containing column depth at xalong, in g/cm^2.

    '''
    betad_fix = np.repeat(beta_arr, 100)  # Repeat each beta 100 times for matching xalong length
    chord_length = 2 * Re * np.sin(np.deg2rad(beta_arr))  # Vectorized trajectory length
    dx = chord_length / 100.0  # Step size for integration

    # Generate xalong using vectorized operations
    xalong = np.repeat(dx, 100) * np.tile(np.arange(1, 101), len(beta_arr))

    # Compute density at all (xalong, beta) points at once
    r_values, rho_values = densityatx(xalong, betad_fix, idepth, model_name)

    # Compute column depth using cumulative integration
    cdalong =  np.asarray([quad(cdtot,0,j,args=(i,idepth,model_name))[0] for i,j in tqdm(zip(betad_fix, xalong), total=len(betad_fix), desc="Integrating cdalong")
    ])
    
    #print("new code = ", betad_fix, xalong, len(xalong), cdalong, len(cdalong))
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

    model_name = ["PREM", "ak135"]

    col_meta = OrderedDict({'Description':'Column trajectories for water layer = %s km' % str(idepth),
                            'beta':'Earth emergence angle, in degrees',
                            'xalong':'Distance in water, in km',
                            'cdalong':'Column depth at xalong, in g/cm^2'})
    
    '''beta, xalong_prem, cdalong_prem = gen_col_trajs(idepth, model_name[0])
    col_table = Table([beta, xalong_prem, cdalong_prem], names=('beta', 'xalong_prem','cdalong_prem'), meta=col_meta)

    Data.add_trajs('col', int(idepth), col_table)'''

    xalong = np.asarray([ Data.get_trajs('col', b, idepth, out=False)[0] for b in beta_arr ])
    cdalong = np.asarray([ Data.get_trajs('col', b, idepth, out=False)[1] for b in beta_arr ])
    
    beta, _, cdalong_ak135 = gen_col_trajs(idepth, model_name[1])
    #col_table_ak135 = Table([xalong_ak135, cdalong_ak135], names=('xalong_ak135','cdalong_ak135'), meta=col_meta)

    col_table = Table([beta, xalong.flatten(), cdalong.flatten(), cdalong_ak135], names=('beta', 'xalong','cdalong_prem','cdalong_ak135'), meta=col_meta)
    Data.add_trajs('col', int(idepth), col_table)
    
    '''beta_water, chord, water = gen_water_trajs(idepth)
    water_meta = OrderedDict({'Description':'Water trajectories for water layer = %s km' % str(idepth),
                              'beta':'Earth emergence angle, in degrees',
                              'chord':'Chord length, in km',
                              'water':'Final water layer distance, in km'})
    water_table = Table([beta_water, chord, water], names=('beta','chord','water'), meta=water_meta)

    Data.add_trajs('water', int(idepth), water_table) # yikes! fixed!'''
    return None

if __name__ == "__main__":
    #pass

    #Rin = np.arange(0, 6372, 1)
    #Rin = np.flip(Rin)

    '''den = ak135density(Rin,3)  
    den_prem = premdensity(Rin,3) 

    plt.plot(Rin, den, label="ak135")
    plt.plot(Rin, den_prem, label="PREM")
    plt.legend()
    plt.show()'''

    '''beta_deg = np.array([1, 10, 20, 30, 40, 50])
    cols_prem = PREMgramVSang(beta_deg, 10, "PREM")
    cols_ak135 = PREMgramVSang(beta_deg, 10, "ak135")
    plt.plot(beta_deg, cols_prem/1e5)
    plt.plot(beta_deg, cols_ak135/1e5)
    plt.legend()
    plt.show()'''

    create_traj_table(0)
    
    #gen_col_trajs(0, "PREM")
