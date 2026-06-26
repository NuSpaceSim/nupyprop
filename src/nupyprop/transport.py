#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 15:02:27 2024

@author: Diksha Garg
Comment: This code contains all the functions that are required in the propagation process of the particles. 
         These functions will be called by the propagate.py code.
"""
import numpy as np
import nupyprop.constants as const

N_A = const.N_A # Avogadro's number
rho_rock = const.rho_rock # rock density
rho_iron = const.rho_iron # iron density

def cd2distd(xalong, cdalong, col_depth):
    '''
    Interpolate between column depth & distance in a medium.

    Parameters
    ----------
    xalong : float
        Array of distance in the medium, in km.
    cdalong : float
        Array of column depth at xalong, in g/cm^2.
    col_depth : float
        Column depth to interpolate at, in g/cm^2.

    Returns
    ----------
    final_distance : float
        Interpolated distance in the medium, in km.
    '''
    # Interpolate for values within bounds
    interp_result = np.interp(col_depth, cdalong, xalong)
    
    # For values below the lowest cdalong, linear extrapolation from the first segment
    below = (col_depth < np.min(cdalong))
    below_value = (col_depth[below] / cdalong[0]) * xalong[0]  # vectorized

    # For values above the highest cdalong, assign max xalong
    above = (col_depth > np.max(cdalong))
    above_value = np.full(np.sum(above), np.max(xalong))
    
    # Prepare output array
    final_distance = interp_result.copy()
    final_distance[below] = below_value
    final_distance[above] = above_value
    
    return final_distance

def get_rho_frac(rho):
    '''
    Calculates the correction/scaling fraction for material density between rock & iron.

    Parameters
    ----------
    rho : float
        Density of material, in g/cm^3.

    Returns
    -------
    frac : float
        Scaling factor for density change between rock & iron (for Bremmstrahlung & pair production).
    frac_pn : float
        Scaling factor for density change between rock & iron (for photonuclear).
   '''
    rho = np.asarray(rho, dtype=float)

    # Initialize output arrays with default values
    frac = np.ones_like(rho)
    frac_pn = np.ones_like(rho)

    # Mask for densities between rock and iron
    mask = (rho > rho_rock) & (rho < rho_iron)

    if np.any(mask):
        f_rock = (rho_iron - rho[mask]) / (rho_iron - rho_rock)
        frac[mask] = 1.97 - 0.97 * f_rock
        frac_pn[mask] = 0.91 + 0.09 * f_rock

    return frac, frac_pn

def int_xc_nu(energy, nu_xc, nu_bsm_xc, fac_nu, E_nu):
    '''
    Interpolate between neutrino energy & cross-section values.

    Parameters
    ----------
    energy : float array
        Energy to interpolate at, in GeV.
    nu_xc : numpy.ndarray
        2D array of neutrino cross-section values.
    nu_bsm_xc : np.ndarray 
        2D array containing neutrino BSM CC & NC cross-section values, in cm^2.
    fac_nu : float
        Rescaling factor for SM neutrino cross-sections.
    E_nu : numpy.ndarray
        Array of neutrino energy

    Returns
    ----------
    sig_cc : float array
        Interpolated CC cross-section value, in cm^2.
    sig_nc : float array
        Interpolated NC cross-section value, in cm^2.
    '''

    # Perform interpolation
    sig_cc = np.interp(energy, E_nu, nu_xc[:, 0])
    sig_nc = np.interp(energy, E_nu, nu_xc[:, 1])
    sig_bsm = np.interp(energy, E_nu, nu_bsm_xc[:, 0])

    # Apply rescaling factor
    sig_cc *= fac_nu
    sig_nc *= fac_nu
    sig_bsm *= fac_nu

    return sig_cc, sig_nc, sig_bsm

def int_xc_lep(energy, xc_arr, rho, E_lep):
    '''
    Interpolate between lepton energy & cross-section values.

    Parameters
    -----------
    energy : float
        Energy value to interpolate at, in GeV.
    xc_arr : numpy.ndarray
        2D array of N_A/A*charged lepton cross-section values, in cm^2/g.
    rho : float
        Density of the material, in g/cm^3.
    E_lep : np.ndarray
        Array of charged lepton energies, in GeV

    Returns
    -------
    sig_brem : float
        Interpolated cross-section value for Bremmstrahlung, in cm^2/g.
    sig_pair : float
        Interpolated cross-section value for pair production, in cm^2/g.
    sig_pn : float
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

def int_alpha(energy, alpha_sig, E_lep):
    '''
    Interpolate between charged lepton energy & ionization energy loss values.

    Parameters
    ----------
    energy : float
        Energy value to interpolate at, in GeV.
    alpha_sig : numpy.ndarray
        1D array of ionization energy loss, in (GeV*cm^2)/g.
    E_lep : np.ndarray
        Array of charged lepton energies, in GeV

    Returns
    ----------
    alpha : float
        Interpolated ionization energy loss value, in (GeV*cm^2)/g.
    '''
    # Perform interpolation
    alpha = np.interp(energy, E_lep, alpha_sig)
    return alpha

def int_beta(energy, beta_arr, rho, E_lep):
    '''
    Interpolate between charged lepton energy & beta (energy loss parameter) values.

    Parameters
    -----------
    energy : float
        Energy value to interpolate at, in GeV.
    beta_arr : numpy.ndarray
        2D array of beta values, in cm^2/g.
    rho : float
        Density of the material, in g/cm^3.
    E_lep : np.ndarray
        Array of charged lepton energies, in GeV

    Returns
    ---------
    tot : float
        Interpolated (& summed) value of beta, in cm^2/g.
    '''
    # Get scaling factors based on material density
    frac, frac_pn = get_rho_frac(rho)

    # Perform interpolation for each beta component
    brem = np.interp(energy, E_lep, beta_arr[:, 0])
    pair = np.interp(energy, E_lep, beta_arr[:, 1])
    pn = np.interp(energy, E_lep, beta_arr[:, 2])

    # Calculate total beta value with scaling
    tot = (frac * brem) + (frac * pair) + (frac_pn * pn)
    return tot

def idecay(energy, distance, m_le, c_tau):
    '''
    Calculate decay probability of lepton.

    Parameters
    ----------
    energy : float 
        Charged lepton energy, in GeV.
    distance : float 
        Distance of charged lepton travel, in cm.
    m_le : float
        Mass of charged lepton, in GeV.
    c_tau : float 
        Decay length of charged lepton, in cm.

    Returns
    -------
    Decay : int
        Decay status (0 means the charged lepton decayed, 1 means it did not decay).
    '''
    # Calculate the Lorentz factor
    gamma_val = energy / m_le

    # Calculate the decay probability
    prob_decay = 1.0 - np.exp(-distance / (gamma_val * c_tau))

    # Generate a random number between 0 and 1
    dy = np.random.random(size=prob_decay.shape)

    # Determine the decay status based on the probability
    decay = (dy >= prob_decay).astype(int)  # 0 = decayed, 1 = survived

    return decay

def em_cont_part(E_init, alpha_val, beta_val, x, m_le):
    '''
    Calculate the charged lepton electromagnetic energy loss (continuous part) a la MUSIC.

    Parameters
    ----------
    E_init : float 
        Initial charged lepton energy, in GeV.
    alpha_val : float
        Ionization energy loss value, in (GeV*cm^2)/g.
    beta_val : float 
        Energy loss parameter (brem + pair + pn), in cm^2/g.
    x : float
        Distance (column depth) of charged lepton travel, in g/cm^2.
    m_le : float
        Mass of charged lepton, in GeV.

    Returns
    -------
    E_fin : float
        Final charged lepton energy, in GeV.
    '''
    bx = beta_val * x

    # Use np.where to handle small-beta*x regime elementwise
    E_fin = np.where(
        bx < 1e-6,
        E_init * (1 - bx) - alpha_val * x,
        E_init * np.exp(-bx) - alpha_val / beta_val * (1 - np.exp(-bx))
    )
    
    # Ensure E_fin is not below lepton rest mass
    E_fin = np.maximum(E_fin, m_le)

    return E_fin

def int_depth_nu(energy, nu_xc, nu_bsm_xc, fac_nu, E_nu):
    '''
    Calculate neutrino interaction depth.
    int_depth = M/(N_A*sigma_tot)

    Parameters
    ----------
    energy : float array
        Neutrino energy, in GeV.
    fac_nu : float
        Rescaling factor for SM neutrino cross-sections.
    nu_xc : np.ndarray
        2D array containing neutrino CC & NC cross-section values, in cm^2.
    nu_bsm_xc : np.ndarray 
        2D array containing neutrino BSM CC & NC cross-section values, in cm^2.
    E_nu : numpy.ndarray
        Array of neutrino energy

    Returns
    ----------
    x_int : float array
        Neutrino interaction depth, in cm^2/g.
    '''
    sig_cc, sig_nc, sig_bsm = int_xc_nu(energy, nu_xc, nu_bsm_xc, fac_nu, E_nu)
    sig_weak = sig_cc  + sig_nc + sig_bsm # weak interactions 
    x_int = 1.0 / (N_A * sig_weak)

    return x_int

def int_depth_lep(energy, xc_arr, rho, m_le, c_tau, E_lep):
    '''
    Calculate charged lepton interaction depth.
    int_depth = M/((N_A/A)*sigma_tot + 1/(gamma*c*tau*rho)); here we need rho to convert decay distance to decay depth

    Parameters
    ----------
    energy : float
        Charged lepton energy, in GeV.
    xc_arr : np.ndarray
        2D array containing N_A/A*charged lepton-nucleon cross-section values, in cm^2/g.
    rho : float
         Density of material, in g/cm^3.
    m_le float
        Mass of charged lepton, in GeV.
    c_tau : float
        Decay length of charged lepton, in cm.
    E_lep : np.ndarray
        Array of charged lepton energies, in GeV

    Returns
    ----------
    x_int : float
        Charged lepton interaction length, in g/cm^2.
    '''
    # Initialize CC and NC cross-sections
    sig_cc, sig_nc = 0.0, 0.0

    # Calculate decay length and decay depth inverse
    decay_length = (energy / m_le) * c_tau
    decay_depth_inv = 1.0 / (decay_length * rho)

    # Get interpolated cross-section values for brem, pair production, and photonuclear interactions
    sig_brem, sig_pair, sig_pn = int_xc_lep(energy, xc_arr, rho, E_lep)

    # Calculate total EM and weak interactions
    sig_em = sig_brem + sig_pair + sig_pn
    sig_weak = sig_cc + sig_nc

    # Calculate the interaction depth
    x_int = 1.0 / (sig_em + sig_weak + decay_depth_inv)

    return x_int

def interaction_type_nu(energy, nu_xc, nu_bsm_xc, fac_nu, E_nu):
    '''
    Determine the type of neutrino-nucleon interaction.

    Parameters
    ----------
    energy : float array
        Neutrino energy, in GeV.
    nu_xc : np.ndarray
        2D array containing neutrino CC & NC cross-section values, in cm^2.
    nu_bsm_xc : np.ndarray 
        2D array containing neutrino BSM CC & NC cross-section values, in cm^2.
    fac_nu : float
        Rescaling factor for SM neutrino cross-sections.
    E_nu : numpy.ndarray
        Array of neutrino energy

    Returns
    ----------
    int_type : int array
        Type of neutrino interaction. 0=CC; 1=NC.
    '''
    # Interpolate CC & NC cross-section values
    sig_cc, sig_nc, sig_bsm = int_xc_nu(energy, nu_xc, nu_bsm_xc, fac_nu, E_nu)

    # Calculate the total and CC fraction of the cross-section
    tot_frac = sig_cc  + sig_bsm + sig_nc
    cc_frac = (sig_cc + sig_bsm)/ tot_frac
    bsm_frac = sig_bsm / tot_frac

    # Generate a random number
    x = np.random.random(size=energy.shape)

    # Determine the interaction type based on the random number
    int_type = np.empty_like(x, dtype=int)
    int_type[x < bsm_frac] = 2      # BSM [0, bsm_frac)
    int_type[(x >= bsm_frac) & (x < cc_frac)] = 0    # CC [bsm_frac, bsm_frac+cc_frac)
    int_type[x >= cc_frac] = 1      # NC [cc_frac, 1)
    #int_type = (x > cc_frac).astype(int)  # 0 for CC, 1 for NC (vectorized)

    return int_type

def interaction_type_lep(energy, xc_arr, rho, m_le, c_tau, E_lep):
    '''
    Determine the type of charged lepton-nucleon interaction.

    Parameters
    ----------
    energy : float
        Charged lepton energy, in GeV.
    xc_arr : np.ndarray
        2D array containing N_A/A*charged lepton-nucleon cross-section values, in cm^2/g.
    rho : float
        Density of material, in g/cm^3.
    m_le : float
        Mass of charged lepton, in GeV.
    c_tau : float
        Decay length of charged lepton, in cm.
    E_lep : np.ndarray
        Array of charged lepton energies, in GeV

    Returns
    -------
    int_type : int
        Type of lepton interaction. 2=decay; 3=Bremmstrahlung; 4=pair-production; 5=photonuclear; 6=CC/NC (placeholder).
    '''
    # Placeholders for CC and NC interactions (values can be read from a lookup table in the future)
    sig_cc, sig_nc = 0.0, 0.0

    # Interpolate cross-section values for Bremmstrahlung, pair-production, and photonuclear interactions
    sig_brem, sig_pair, sig_pn = int_xc_lep(energy, xc_arr, rho, E_lep)

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
    y = np.random.random(size=energy.shape)

    # cumulative fractions
    cum_decay = decay_frac
    cum_brem = decay_frac + brem_frac
    cum_pair = cum_brem + pair_frac
    cum_pn = cum_pair + pn_frac

    # default type (6 = CC/NC placeholder)
    int_type = np.full(energy.shape, 6, dtype=int)

    # assign based on thresholds
    int_type[y <= cum_decay] = 2
    int_type[(y > cum_decay) & (y <= cum_brem)] = 3
    int_type[(y > cum_brem) & (y <= cum_pair)] = 4
    int_type[(y > cum_pair) & (y <= cum_pn)] = 5

    return int_type

def find_y(energy, ixc_arr, ip, E_nu, E_lep, yvals, ixc_bsm_arr=None):
    '''
    Stochastic determination of neutrino/lepton inelasticity.

    Parameters
    ----------
    energy : float array
        Neutrino or charged lepton energy, in GeV.
    ixc_arr : np.ndarray
        Neutrino or charged lepton integrated cross-section CDF values.
    ip : int array
        Type of neutrino-nucleon or lepton-nucleon interaction.
        0=CC, 1=NC, 2=BSM, 3=brem, 4=pair, 5=pn
    E_nu : np.ndarray
        Array of neutrino energies (for neutrinos).
    E_lep : np.ndarray
        Array of charged lepton energies (for charged leptons).
    yvals : np.ndarray
        Array of min. y values from which the cross-section CDF is calculated.
    ixc_bsm_arr : np.ndarray or None
        Neutrino BSM integrated cross-section CDF values. None if no BSM.

    Returns
    ----------
    y : float array
        Inelasticity, y = (E_init-E_final)/E_initial.
    '''
    is_sm_nu  = (ip == 0) | (ip == 1)
    is_bsm_nu = (ip == 2)
    is_lep    = ~(is_sm_nu | is_bsm_nu)

    ny = yvals.size
    y = np.empty_like(energy, dtype=float)

    def invert_cdf(cdf, u):
        ''' Vectorized CDF inversion. cdf shape: (Nevents, Ny), u shape: (Nevents,) '''
        k = np.sum(cdf < u[:, None], axis=1).astype(int)
        k = np.clip(k, 1, ny - 1)
        row = np.arange(cdf.shape[0])
        c0 = cdf[row, k - 1]
        c1 = cdf[row, k]
        y0 = yvals[k - 1]
        y1 = yvals[k]
        den = c1 - c0
        t = np.where(den > 0.0, (u - c0) / den, 0.0)
        return y0 + t * (y1 - y0)

    # --- SM neutrino part (CC=0, NC=1) ---
    if np.any(is_sm_nu):
        e = np.clip(np.searchsorted(E_nu, energy[is_sm_nu]), 0, E_nu.size - 1)
        i = ip[is_sm_nu]  # 0 or 1
        cdf = ixc_arr[:, e, i].T  # shape (Nevents, Ny)
        u = np.random.random(size=e.shape)
        y[is_sm_nu] = invert_cdf(cdf, u)

    # --- BSM neutrino part (ip=2) ---
    if np.any(is_bsm_nu) and ixc_bsm_arr is not None:
        e = np.clip(np.searchsorted(E_nu, energy[is_bsm_nu]), 0, E_nu.size - 1)
        i = np.zeros_like(e, dtype=int)  # single BSM channel
        cdf = ixc_bsm_arr[:, e, i].T  # shape (Nevents, Ny)
        u = np.random.random(size=e.shape)
        y[is_bsm_nu] = invert_cdf(cdf, u)

    # --- Charged lepton part (brem=3, pair=4, pn=5) ---
    if np.any(is_lep):
        e = np.clip(np.searchsorted(E_lep, energy[is_lep]), 0, E_lep.size - 1)
        i = (ip[is_lep] - 3).astype(int)  # 0=brem, 1=pair, 2=pn
        cdf = ixc_arr[:, e, i].T  # shape (Nevents, Ny)
        u = np.random.random(size=e.shape)
        y[is_lep] = invert_cdf(cdf, u)

    y = np.minimum(y, 1.0)
    return y

def polarization(y, pin, theta_in, ypol, Pcthp, P):
    '''
    Calculate polarization after electromagnetic scattering of tau-leptons.
    Based on https://arxiv.org/pdf/2205.05629, D. Garg, et al.

    Parameters
    -----------
    y : float
        Inelasticity, y = (E_init - E_final) / E_init.
    pin : float
        Initial momentum magnitude.
    theta_in : float
        Initial polar angle in radians.
    ypol : float array
        Predefined inelasticity array
    Pcthp :
        np.cos(theta_P), where theta_P is the polar angle of the spin vector in tau rest frame
    P : float array
        Magnitude of polarization vector, defining the degree of polarization

    Returns
    ---------
    pout : float
        Final momentum magnitude.
    theta_out : float
        Final polar angle in radians.
    '''
    pzout = np.interp(y, ypol, Pcthp)
    pout_factor = np.interp(y, ypol, P)

    # Initialize output arrays
    pout = np.empty_like(pin)
    theta_out = np.empty_like(theta_in)

    # Case 1: y < 0.01  → no polarization change
    small_y_mask = (y < 0.01)
    pout[small_y_mask] = pin[small_y_mask]
    theta_out[small_y_mask] = theta_in[small_y_mask]

    # Case 2: y >= 0.01 → polarization update
    large_y_mask = ~small_y_mask
    if np.any(large_y_mask):
        pzout_large = pzout[large_y_mask]
        pout_large = pout_factor[large_y_mask]
        pin_large = pin[large_y_mask]
        theta_in_large = theta_in[large_y_mask]

        # Compute scattering angle
        cth = pzout_large / pout_large
        cth = np.clip(cth, -1.0, 1.0)
        theta_scatter = np.arccos(cth)

        # Final momentum magnitude
        pout_large = pin_large * pout_large

        # Random sign per event
        r = np.random.choice([1, -1], size=np.count_nonzero(large_y_mask))
        theta_out_large = theta_in_large + r * theta_scatter

        pout[large_y_mask] = pout_large
        theta_out[large_y_mask] = theta_out_large

    return pout, theta_out
