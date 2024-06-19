#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 15:02:27 2024

@author: Diksha Garg
Comment: This code contains all the functions that are required in the propagation process of the particles. These functions will be called by the propagate.py code. 
"""

import numpy as np
import pandas as pd
import importlib_resources
import nupyprop.constants as const 

rho_rock = const.rho_rock
rho_iron = const.rho_iron

E_nu = const.E_nu
E_lep = const.E_lep
yvals = const.yvals

polarization_path = importlib_resources.files('nupyprop.datafiles') / 'polarization_data.txt'
pola_df = pd.read_csv(resource_path, delimiter='\s+', header='#') 

def cd2distd(xalong, cdalong, col_depth):
    '''
    Interpolate between column depth & distance in a medium.

    Args:
    xalong (float):
        Array of distance in the medium, in km.
    cdalong (float):
        Array of column depth at xalong, in g/cm^2.
    col_depth (float):
        Column depth to interpolate at, in g/cm^2.

    Returns:
        float:
            Interpolated distance in the medium, in km.
    '''

    if (col_depth < np.min(cdalong)):
       return (col_depth/cdalong[0])*xalong[0]
    elif (col_depth > np.max(cdalong)):
       return np.max(xalong)
    else:
       distance = np.interp(col_depth, cdalong, xalong)
       return distance 

def get_rho_frac(rho, frac, frac_pn):
    ''' 
    Calculates the correction/scaling fraction for material density between rock & iron.

    Args:
    rho (float):
        Density of material, in g/cm^3.
    
    Returns:
    frac (float):
        Scaling factor for density change between rock & iron (for Bremmstrahlung & pair production).
    frac_pn (float):
        Scaling factor for density change between rock & iron (for photonuclear).
   '''

    if (rho > rho_rock and rho < rho_iron):
       f_rock = (rho_iron - rho)/(rho_iron - rho_rock)
       frac = 1.97 - 0.97 * f_rock
       frac_pn = 0.91 + 0.09 * f_rock
    else: # for rho <= rho_water or rho>=iron (that shouldn't happen!)
       frac = 1
       frac_pn = 1

    return frac, frac_pn

def int_xc_nu(energy, nu_xc, fac_nu):
    '''
    Interpolate between neutrino energy & cross-section values.

    Args:
    energy (float):
        Energy value to interpolate at, in GeV.
    nu_xc (numpy.ndarray):
        2D array of neutrino cross-section values.
    fac_nu (float):
        Rescaling factor for SM neutrino cross-sections.
     
    Returns:
    sig_cc (float):
        Interpolated CC cross-section value, in cm^2.
    sig_nc (float):
        Interpolated NC cross-section value, in cm^2.
    '''

    # Perform interpolation
    sig_cc = np.interp(energy, E_nu, nu_xc[:, 0])
    sig_nc = np.interp(energy, E_nu, nu_xc[:, 1])

    # Apply rescaling factor
    sig_cc *= fac_nu
    sig_nc *= fac_nu

    return sig_cc, sig_nc

def int_xc_lep(energy, xc_arr, rho):
    '''
    Interpolate between lepton energy & cross-section values.

    Args:
    energy (float):
        Energy value to interpolate at, in GeV.
    xc_arr (numpy.ndarray):
        2D array of N_A/A*charged lepton cross-section values, in cm^2/g.
    rho (float):
        Density of the material, in g/cm^3.

    Returns:
    sig_brem (float):
        Interpolated cross-section value for Bremmstrahlung, in cm^2/g.
    sig_pair (float):
        Interpolated cross-section value for pair production, in cm^2/g.
    sig_pn (float):
        Interpolated cross-section value for photonuclear, in cm^2/g.
    '''

    # Get scaling factors based on material density
    frac, frac_pn = get_rho_frac(rho)

    # Note: The lookup tables already have N_A multiplied by lep_xc!
    sig_brem = np.interp(energy, E_lep, xc_arr[:, 0])
    sig_pair = np.interp(energy, E_lep, xc_arr[:, 1])
    sig_pn = np.interp(energy, E_lep, xc_arr[:, 2])

    # Apply scaling factors
    sig_brem *= frac
    sig_pair *= frac
    sig_pn *= frac_pn

    return sig_brem, sig_pair, sig_pn

def int_alpha(energy, alpha_sig):
    """
    Interpolate between charged lepton energy & ionization energy loss values.

    Args:
    energy (float):
        Energy value to interpolate at, in GeV.
    alpha_sig (numpy.ndarray):
        1D array of ionization energy loss, in (GeV*cm^2)/g.

    Returns:
    alpha (float):
        Interpolated ionization energy loss value, in (GeV*cm^2)/g.
    """
    # Perform interpolation
    alpha = np.interp(energy, E_lep, alpha_sig)
    return alpha

def int_beta(energy, beta_arr, rho):
    """
    Interpolate between charged lepton energy & beta (energy loss parameter) values.

    Args:
    energy (float):
        Energy value to interpolate at, in GeV.
    beta_arr (numpy.ndarray):
        2D array of beta values, in cm^2/g.
    rho (float):
        Density of the material, in g/cm^3.
    
    Returns:
    tot (float):
        Interpolated (& summed) value of beta, in cm^2/g.
    """
    # Get scaling factors based on material density
    frac, frac_pn = get_rho_frac(rho)

    # Perform interpolation for each beta component
    brem = np.interp(energy, E_lep, beta_arr[:, 0])
    pair = np.interp(energy, E_lep, beta_arr[:, 1])
    pn = np.interp(energy, E_lep, beta_arr[:, 2])

    # Calculate total beta value with scaling
    tot = (frac * brem) + (frac * pair) + (frac_pn * pn)
    return tot

def searchsorted(array, search_value):
    """
    Given an array and a value, returns the index of the element that is closest to, but less than, the given value.

    Uses a binary search algorithm provided by numpy.

    Args:
    array (numpy.ndarray): Input array.
    search_value (float): Value to search for.

    Returns:
    int: Index of the element closest to, but less than, the given value.
    """
    # Ensure the array is sorted for searchsorted to work correctly
    array = np.sort(array)
    
    # Use numpy's searchsorted to find the insertion point
    index = np.searchsorted(array, search_value, side='right') - 1
    
    # Ensure index is within the valid range
    if index < 0:
        return 0  # Return the first index if the search value is less than the smallest element
    return index

def idecay(energy, distance, m_le, c_tau):
    """
    Calculate decay probability of lepton.

    Args:
        energy (float): Charged lepton energy, in GeV.
        distance (float): Distance of charged lepton travel, in cm.
        m_le (float): Mass of charged lepton, in GeV.
        c_tau (float): Decay length of charged lepton, in cm.

    Returns:
        int: Decay status (0 means the charged lepton decayed, 1 means it did not decay).
    """
    # Calculate the Lorentz factor
    gamma_val = energy / m_le
    
    # Calculate the decay probability
    prob_decay = 1.0 - np.exp(-distance / (gamma_val * c_tau))
    
    # Generate a random number between 0 and 1
    dy = np.random.random()
    
    # Determine the decay status based on the probability
    if dy < prob_decay:
        decay = 0
    else:
        decay = 1
    
    return decay

def em_cont_part(E_init, alpha_val, beta_val, x, m_le):
    """
    Calculate the charged lepton electromagnetic energy loss (continuous part) a la MUSIC.

    Args:
        E_init (float): Initial charged lepton energy, in GeV.
        alpha_val (float): Ionization energy loss value, in (GeV*cm^2)/g.
        beta_val (float): Energy loss parameter (brem + pair + pn), in cm^2/g.
        x (float): Distance (column depth) of charged lepton travel, in g/cm^2.
        m_le (float): Mass of charged lepton, in GeV.

    Returns:
        float: Final charged lepton energy, in GeV.
    """
    if beta_val * x < 1e-6:
        E_fin = E_init * (1 - beta_val * x) - alpha_val * x
    else:
        E_fin = E_init * np.exp(-beta_val * x) - alpha_val / beta_val * (1 - np.exp(-beta_val * x))
    
    if E_fin < 0:
        E_fin = m_le

    return E_fin

def int_depth_nu(energy, nu_xc, fac_nu):
    """
    Calculate neutrino interaction depth.
    int_depth = M/(N_A*sigma_tot)
    
    Args:
        energy (float): Neutrino energy, in GeV.
        fac_nu (float): Rescaling factor for SM neutrino cross-sections.
        nu_xc (np.ndarray): 2D array containing neutrino CC & NC cross-section values, in cm^2.

    Returns:
        float: Neutrino interaction depth, in cm^2/g.
    """
      sig_cc, sig_nc = int_xc_nu(energy, nu_xc, fac_nu)
      sig_weak = sig_cc + sig_nc # weak interactions
      x_int = 1.0 / (N_A * sig_weak)

      return x_int

def int_depth_lep(energy, xc_arr, rho, m_le, c_tau):
    """
    Calculate charged lepton interaction depth.
    int_depth = M/((N_A/A)*sigma_tot + 1/(gamma*c*tau*rho)); here we need rho to convert decay distance to decay depth

    Args:
        energy (float): Charged lepton energy, in GeV.
        xc_arr (np.ndarray): 2D array containing N_A/A*charged lepton-nucleon cross-section values, in cm^2/g.
        rho (float): Density of material, in g/cm^3.
        m_le (float): Mass of charged lepton, in GeV.
        c_tau (float): Decay length of charged lepton, in cm.

    Returns:
        float: Charged lepton interaction length, in g/cm^2.
    """
    # Initialize CC and NC cross-sections
    sig_cc, sig_nc = 0.0, 0.0

    # Calculate decay length and decay depth inverse
    decay_length = (energy / m_le) * c_tau
    decay_depth_inv = 1.0 / (decay_length * rho)

    # Get interpolated cross-section values for brem, pair production, and photonuclear interactions
    sig_brem, sig_pair, sig_pn = int_xc_lep(energy, xc_arr, rho)

    # Calculate total EM and weak interactions
    sig_em = sig_brem + sig_pair + sig_pn
    sig_weak = sig_cc + sig_nc

    # Calculate the interaction depth
    x_int = 1.0 / (sig_em + sig_weak + decay_depth_inv)
    
    return x_int

def interaction_type_nu(energy, nu_xc, fac_nu):
    """
    Determine the type of neutrino-nucleon interaction.

    Args:
        energy (float): Neutrino energy, in GeV.
        nu_xc (np.ndarray): 2D array containing neutrino CC & NC cross-section values, in cm^2.
        fac_nu (float): Rescaling factor for SM neutrino cross-sections.

    Returns:
        int: Type of neutrino interaction. 0=CC; 1=NC.
    """
    # Interpolate CC & NC cross-section values
    sig_cc, sig_nc = int_xc_nu(energy, nu_xc, fac_nu)

    # Calculate the total and CC fraction of the cross-section
    tot_frac = sig_cc + sig_nc
    cc_frac = sig_cc / tot_frac

    # Generate a random number
    x = np.random.random()

    # Determine the interaction type based on the random number
    if x <= cc_frac:
        int_type = 0  # CC
    else:
        int_type = 1  # NC

    return int_type

def interaction_type_lep(energy, xc_arr, rho, m_le, c_tau):
    """
    Determine the type of charged lepton-nucleon interaction.
    
    Args:
        energy (float): Charged lepton energy, in GeV.
        xc_arr (np.ndarray): 2D array containing N_A/A*charged lepton-nucleon cross-section values, in cm^2/g.
        rho (float): Density of material, in g/cm^3.
        m_le (float): Mass of charged lepton, in GeV.
        c_tau (float): Decay length of charged lepton, in cm.

    Returns:
        int: Type of lepton interaction. 2=decay; 3=Bremmstrahlung; 4=pair-production; 5=photonuclear; 6=CC/NC (placeholder).
    """
    # Placeholders for CC and NC interactions (values can be read from a lookup table in the future)
    sig_cc, sig_nc = 0.0, 0.0
    
    # Interpolate cross-section values for Bremmstrahlung, pair-production, and photonuclear interactions
    sig_brem, sig_pair, sig_pn = int_xc_lep(energy, xc_arr, rho)
    
    # Calculate decay length and inverse decay depth
    decay_length = (energy / m_le) * c_tau
    decay_depth_inv = 1.0 / (decay_length * rho)
    
    # Calculate total interaction depth for the lepton
    sig_em = sig_brem + sig_pair + sig_pn
    sig_weak = sig_cc + sig_nc
    int_lep = 1.0 / (sig_em + sig_weak + decay_depth_inv)
    
    # Calculate fractions for different interaction types
    tot_frac = 1.0 / int_lep
    decay_frac = decay_depth_inv / tot_frac
    brem_frac = sig_brem / tot_frac
    pair_frac = sig_pair / tot_frac
    pn_frac = sig_pn / tot_frac
    
    # Generate a random number and determine interaction type based on cumulative fractions
    y = np.random.random()

    # Create a list of (cumulative fraction, interaction type) tuples
    interaction_types = [
        (decay_frac, 2),  # decay
        (decay_frac + brem_frac, 3),  # Bremmstrahlung
        (decay_frac + brem_frac + pair_frac, 4),  # pair-production
        (decay_frac + brem_frac + pair_frac + pn_frac, 5),  # photonuclear
    ]

    # Determine interaction type based on cumulative fractions
    int_type = 6  # Default to placeholder for CC/NC
    for frac, type_id in interaction_types:
        if y <= frac:
            int_type = type_id
            break

    return int_type

def find_y(energy, ixc_arr, ip, E_nu, E_lep, yvals):
    """
    Stochastic determination of neutrino/lepton inelasticity.

    Args:
        energy (float): Neutrino or charged lepton energy, in GeV.
        ixc_arr (np.ndarray): Neutrino or charged lepton integrated cross-section CDF values.
        ip (int): Type of neutrino-nucleon or lepton-nucleon interaction.
        E_nu (np.ndarray): Array of neutrino energies (for neutrinos).
        E_lep (np.ndarray): Array of charged lepton energies (for charged leptons).
        yvals (np.ndarray): Array of min. y values from which the cross-section CDF is calculated.

    Returns:
        float: Inelasticity, y = (E_init-E_final)/E_initial.
    """
    if ip == 0 or ip == 1:  # for neutrinos
        energy_index = searchsorted(E_nu, energy)
        ip_id = 0 if ip == 0 else 1 #ip_id=0 for CC and 1 for NC 
    else:  # for charged leptons
        energy_index = searchsorted(E_lep, energy)
        ip_id = ip - 2  # Convert ip (3, 4, 5) to index (1, 2, 3) for brem, PP, PN, respectively
        
    search_arr = ixc_arr[:, energy_index, ip_id]
    y = np.interp(dy, search_arr, yvals)
    # dy is the randomly sampled cross-section CDF value (between 0 & 1)
    # search_arr = cross-section CDF value array for energy_index
    # yvals = array of min. y values from which the cross-section CDF is calculated (see models.py for calculation details)
    # y is the interpolated (yvals) value corresponding to the cross-section CDF value = dy; this y is responsible for stochastic energy losses
    
    if y > 1.0:
        y = 1.0
              
    return y

    subroutine polarization(y, pin, theta_in, pout, theta_out)

      implicit none

      real(dp), intent(in) :: pin, theta_in, y
      real(dp), intent(out) :: pout, theta_out

      real(dp) :: pzout, theta, p0, rs, cth

      if (y<0.01) then
         theta_out = theta_in
         pout = pin
         return
      else
         call interpol(y, ypol, Pcthp, pzout)  !!interpolating Pcthp value for y
         call interpol(y, ypol, P, pout)   !!interpolating P value for y

      end if

      cth = pzout/pout
      theta_out = acos(cth)
      p0 = pin*pout
      pout = p0

      call random(rs)

      theta = theta_in + rs*theta_out
      theta_out = theta

    end subroutine polarization

def polarization(y, pin, theta_in):
    '''
    Calculate polarization after scattering.

    Args:
    y (float): Inelasticity, y = (E_init - E_final) / E_init.
    pin (float): Initial momentum magnitude.
    theta_in (float): Initial polar angle in radians.

    Returns:
    float: pout, Final momentum magnitude.
    float: theta_out, Final polar angle in radians.
    '''
    if y < 0.01:
        theta_out = theta_in
        pout = pin
        return pout, theta_out
    else:
        # Interpolate Pcthp and P values for given y
        pzout = np.interp(y, ypol, Pcthp)
        pout = np.interp(y, ypol, P)

    # Calculate cosine of scattering angle
    cth = pzout / pout
    theta_out = np.arccos(cth)

    # Update momentum
    p0 = pin * pout
    pout = p0

    # Random scattering angle adjustment
    rs = np.random.random()

    if (rs<0.5):
        r = 1
    else:
        r = -1
        
    theta = theta_in + r * theta_out
    theta_out = theta

    return pout, theta_out
