#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 3 14:17:24 2020

@author: sam
"""
# import matplotlib as mpl
# mpl.use('Agg') # for clusters

import nupyprop.data as Data
import nupyprop.cross_section as Cross_section

import numpy as np
import pandas as pd
import scipy.integrate as integrate
import scipy.constants as scc
import multiprocessing as mp
from multiprocessing import Pool
import matplotlib.pyplot as plt
import pickle
import random
# import quadpy
import time

m_e = scc.physical_constants["electron mass energy equivalent in MeV"][0]*1e-3 # GeV
m_mu = scc.physical_constants["muon mass energy equivalent in MeV"][0]*1e-3 # GeV
m_tau = scc.physical_constants["tau mass energy equivalent in MeV"][0]*1e-3 # GeV
m_pi = 139.57018e-3 # pion mass in GeV
alpha_fs = scc.fine_structure
le = 3.8616e-11 # electron Compton wavelength in cm

E_lep = Data.E_lep
m_p = Cross_section.m_p
N_A = Cross_section.N_A

def integrator(args):
    '''

    Parameters
    ----------
    args : list of lists
        Function arguments.

    Returns
    -------
    function
        Integrator function.

    '''
    fun = args[0]
    low_lim = args[1]
    up_lim = args[2]
    arg = args[3]
    return integrate.quad(fun,low_lim,up_lim,args=(arg,))

def nquad_integrator(args):
    '''

    Parameters
    ----------
    args : list of lists
        Function arguments.

    Returns
    -------
    function
        Integrator function.

    '''
    fun = args[0]
    lim_1 = args[1] # dep limit
    lim_2 = args[2] # indep limit
    arg = args[3]
    return integrate.nquad(fun, [lim_1,lim_2], args=(arg,))

def dbl_integrator(args):
    '''

    Parameters
    ----------
    args : list of lists
        Function arguments.

    Returns
    -------
    function
        Integrator function.

    '''
    fun = args[0]
    indep_lim_low = args[1] # has to be a value and not a function, but in the form of a function! See pair and PN limits
    indep_lim_high = args[2] # has to be a value and not a function, but in the form of a function!
    dep_lim_low = args[3]
    dep_lim_high = args[4]
    var = args[5] # the variable we want to loop over (for eg. energy)
    return integrate.dblquad(fun, indep_lim_low(var), indep_lim_high(var), dep_lim_low, dep_lim_high, args=(var,), epsabs=1e-12)[0]

def rep(val): # Remove non-physical values
    '''

    Parameters
    ----------
    val : float
        Usually negative or np.nan values, to set to 0.

    Returns
    -------
    float
        Returns 0 in case of non-physical entries or returns val.

    '''
    if (val<0) or (np.isnan(val)==True):
        return 0
    else:
        return val
    return 'Problem in rep function'

def short_int_nquad(fn, indep_fn_lim_high, dep_fn_lim, arr, arg):
    '''

    Parameters
    ----------
    fn : NONE
        Lamba function (integrating function).
    indep_fn_lim_high : NONE
        Lambda function (independent function, upper limit).
    dep_fn_lim : NONE
        Lambda function (dependent function limit).
    arr : ndarray
        1D array containing y-values you want to integrate over.
    arg : Float
        Integrating variable (usually energy).

    Returns
    -------
    ndarray
        1D array containing integrated results.

    '''
    ans = []
    result = integrate.nquad(fn, [dep_fn_lim,[arr[0],indep_fn_lim_high(arg)]],args=(arg,))[0]
    ans.append(result)
    for i in range(len(arr)-1):
        result+= integrate.nquad(fn, [dep_fn_lim,[arr[i+1],arr[i]]],args=(arg,))[0]
        ans.append(result)
    return np.asarray(ans)

def short_int_dblquad(fn, indep_fn_lim_high, dep_fn_lim_low, dep_fn_lim_high, arr, arg): #indep_fn_lim_high has to return a value and not a function!
    '''

    Parameters
    ----------
    fn : NONE
        Lamba function (integrating function).
    indep_fn_lim_high : float
        Independent function, upper limit value.
    dep_fn_lim_low : NONE
        Lambda function (dependent function, lower limit).
    dep_fn_lim_high : NONE
        Lambda function (dependent function, upper limit).
    arr : ndarray
        1D array containing y-values you want to integrate over.
    arg : Float
        Integrating variable (usually energy).

    Returns
    -------
    ndarray
        1D array containing integrated results.

    '''
    ans = []
    result = integrate.dblquad(fn, arr[0],indep_fn_lim_high(arg), dep_fn_lim_low, dep_fn_lim_high, args=(arg,), epsabs=1e-12)[0]
    ans.append(result)
    for i in range(len(arr)-1):
        result+= integrate.dblquad(fn, arr[i+1],arr[i], dep_fn_lim_low, dep_fn_lim_high, args=(arg,), epsabs=1e-12)[0]
        ans.append(result)
    return ans

# =============================================================================
# Ionization Energy Loss
# =============================================================================
def alpha_i(E):
    '''

    Parameters
    ----------
    E : float
        Energy of the lepton, in GeV.

    Returns
    -------
    dEdX : float
        Rate of change of energy with column depth, in (GeV*cm^2)/g.

    '''
    param = {'water':{'I':7.5e-8, 'C':-3.502, 'X_0':0.240, 'X_1':2.8, 'a':0.091, 'm':3.477}} # Note: I has to be in GeV
    param.update({'rock':{'I':1.364e-7, 'C':-3.774, 'X_0':0.049, 'X_1':3.055, 'a':0.083, 'm':3.412}})
    param.update({'iron':{'I':2.86e-7, 'C':-4.291, 'X_0':-0.0012, 'X_1':3.153, 'a':0.147, 'm':2.963}})

    if E<m_le:
        p = 0
    else:
        p = np.sqrt(E**2-m_le**2)

    E_m = 2*m_e*((p**2)/(m_e**2 + m_le**2 + 2*m_e*E))

    gamma = E/m_le
    beta = p/E

    X = np.log10(beta*gamma)
    X_0 = param[material]['X_0']
    X_1 = param[material]['X_1']
    a = param[material]['a']
    m = param[material]['m']
    C = param[material]['C']
    I = param[material]['I']

    if X_0<X<X_1:
        delta = 4.6052*X + a*(X_1-X)**m + C
    else:
        delta = 4.6052*X + C

    dEdX = alpha_fs**2 * 2*np.pi*N_A*le**2*(z*m_e)/(A*beta**2)*(np.log((2*m_e*beta**2*gamma**2*E_m)/I**2) - 2*beta**2 + E_m**2/(4*E**2) - delta)

    return dEdX


# =============================================================================
#      Bremmstrahlung Energy Loss
# =============================================================================
def brem(y, E): # eq. A5
    '''

    Parameters
    ----------
    y : float
        y-value.
    E : float
        Energy of the lepton, in GeV.

    Returns
    -------
    float
        Bremmstrahlung energy loss value.

    '''
    y_max = 1 - (3*m_le/(4*E)) * np.sqrt(np.e) * z**(1/3)
    if y>y_max:
        y_dsig_dy = 0
    else:
        delta = (m_le**2*y)/(2*E*(1-y))
        if z > 10:
            phi = np.log((2/3 * 189/m_e * m_le * z**(-2/3))/(1 + 189 * np.sqrt(np.e)/m_e * delta * z**(-1/3)))
        else:
            phi = np.log((189 * m_le/m_e * z**(-1/3))/(1 + 189* np.sqrt(np.e)/m_e * delta * z**(-1/3)))

        brem = alpha_fs**3 * (4*z*(z+1)*le**2) * (m_e**2/m_le**2) * (1/y) * (4/3 - (4*y)/3 + y**2)*phi

        y_dsig_dy = brem * y

    return (N_A/A) * y_dsig_dy

def cs_brem(y, E):
    '''

    Parameters
    ----------
    y : float
        Integration variable.
    E : float
        Energy, in GeV.

    Returns
    -------
    float
        Pair production cross-section value.

    '''
    return brem(y,E)/y

def brem_bb_high(E):
    '''

    Parameters
    ----------
    E : float
        Energy, in GeV.

    Returns
    -------
    float
        Bremsstrahlung integration upper limit.

    '''
    return 1 - (3*m_le/(4*E)) * np.sqrt(np.e) * z**(1/3)
# =============================================================================
#     Pair Production Energy Loss
# =============================================================================

def pair(rho, y, E): # eq. A9
    '''

    Parameters
    ----------
    rho : float
        Integration variable.
    y : float
        Integration variable.
    E : float
        Energy, in GeV.

    Returns
    -------
    float
        Pair production energy loss value.

    '''
    R = 189
    beta = y**2/(2*(1-y))
    xi = ((m_le * y)/(2*m_e))**2 * (1 - rho**2)/(1 - y)

    y_l = (4 + rho**2 + 3*beta*(1 + rho**2))/((1 + rho**2) * (3/2 + 2*beta) * np.log(3 + xi) + 1 - (3/2)*rho**2)

    y_e = (5 - rho**2 + 4*beta*(1 + rho**2))/(2*(1 + 3*beta) * np.log(3 + 1/xi) - rho**2 - 2*beta*(2 - rho**2))

    L_l = np.log((R * z**(-2/3) * 2/3 * m_le/m_e)/(1 + (2*m_e*np.sqrt(np.e)*R * z**(-1/3)) * (1 + xi) * (1 + y_l)/(E * y * (1 - rho**2))))

    L_e = np.log((R * z**(-1/3) * np.sqrt((1 + xi) * (1 + y_e)))/(1 + (2*m_e*np.sqrt(np.e)*R * z**(-1/3) * (1 + xi) * (1 + y_e))/(E * y * (1 - rho**2)))) - 1/2 * np.log(1 + (3/2 * m_e/m_le * z**(1/3))**2 * (1 + xi) * (1 + y_e))

    phi_l = (((1 + rho**2) * (1 + (3/2)*beta) - 1/xi * (1 + 2*beta) * (1 - rho**2)) * np.log(1 + xi) + (xi * (1 - rho**2 - beta))/(1 + xi) + (1 + 2*beta) * (1 - rho**2)) * L_l

    phi_e = (((2 + rho**2) * (1 + beta) + xi*(3 + rho**2)) * np.log(1 + 1/xi) + (1 - rho**2 - beta)/(1 + xi) - (3 + rho**2)) * L_e

    y_d2sigma_dydrho = y * alpha_fs**4 * 2/(3 * np.pi) *z*(z+1) *(le)**2 * (1-y)/y * (phi_e + (m_e**2/m_le**2) * phi_l)

    return (N_A/A) * y_d2sigma_dydrho

def cs_pair(rho, y, E):
    '''

    Parameters
    ----------
    rho : float
        Integration variable.
    y : float
        Integration variable.
    E : float
        Energy, in GeV.

    Returns
    -------
    float
        Pair production cross-section value.

    '''
    return pair(rho, y, E)/y

def pair_rho_cut(y, E):
    '''

    Parameters
    ----------
    y : float
        Integration variable.
    E : float
        Energy, in GeV.

    Returns
    -------
    float
        Pair production rho limit, for beta_cut.

    '''
    return [-(1 - (6*m_le**2)/(E**2 * (1-y))) * np.sqrt(1 - (4 * m_e)/(E * y)), (1 - (6*m_le**2)/(E**2 * (1-y))) * np.sqrt(1 - (4 * m_e)/(E * y))] # [low, high]

def pair_rho_tot(y, E):
    '''

    Parameters
    ----------
    y : float
        Integration variable.
    E : float
        Energy, in GeV.

    Returns
    -------
    float
        Pair production rho limit, for beta_total.

    '''
    return [-(1 - (6*m_le**2)/(E**2 * (1-y))) * np.sqrt(1 - (4 * m_e)/(E * y)), (1 - (6*m_le**2)/(E**2 * (1-y))) * np.sqrt(1 - (4 * m_e)/(E * y))] # [low, high]


def pair_y_cut(E):
    '''

    Parameters
    ----------
    E : float
        Energy, in GeV.

    Returns
    -------
    float
        Pair production y limit, for beta_cut.

    '''
    return [(4*m_e)/E, 1e-3] # [low, high]

def pair_y_tot(E):
    '''

    Parameters
    ----------
    E : float
        Energy, in GeV.

    Returns
    -------
    float
        Pair production y limit, for beta_total.

    '''
    return [(4*m_e)/E, 1 - (3*m_le/(4*E)) * np.sqrt(np.e) * z**(1/3)] # [low, high]


# =============================================================================
#     Photonuclear Energy Loss (Bezrukov-Bugaev)
# =============================================================================
def pn_bb(y, E): # eq. A12
    '''

    Parameters
    ----------
    y : float
        Integration variable.
    E : float
        Energy, in GeV.

    Returns
    -------
    float
        Photonuclear energy loss (Bezrukov-Bugaev) energy loss value.

    '''
    m1_sq = 0.54 # Gev^2
    m1 = np.sqrt(m1_sq)
    m2_sq = 1.8 # Gev^2
    m2 = np.sqrt(m2_sq)
    k = 1 - 2/y + 2/y**2
    t = (m_le**2 * y**2)/(1-y)
    sigma = 114.3 + 1.647 * (np.log(0.0213 * E))**2
    x = 0.00282 * A**(1/3) * sigma
    G = 3/(x**3) * ((x**2)/2 - 1 + np.e**(-x) * (1 + x))
    t1 = 3/4 * G * (k * np.log(1 + m1/t) - (k * m1_sq)/(m1_sq + t) - (2 * m_le**2)/t)
    t2 = 1/4 * (k * np.log(1 + m2/t) - (2 * m_le**2)/t)
    t3 = m_le**2/(2*t) * (3/4 * G * m1_sq/(m1_sq + t) + 1/4 * m2_sq/t * np.log(1 + t/m2_sq))
    y_dsig_dy = y * alpha_fs/(2*np.pi) * A * y * (sigma* (t1 + t2 + t3))*1e-30
    return (N_A/A) * y_dsig_dy

def cs_pn_bb(y, E):
    '''

    Parameters
    ----------
    y : float
        Integration variable.
    E : float
        Energy, in GeV.

    Returns
    -------
    float
        Photonuclear energy loss (Bezrukov-Bugaev) cross-section value.

    '''
    return pn_bb(y, E)/y

# =============================================================================
#     Photonuclear Energy Loss (ALLM/BDHM/CKMT/Custom)
# =============================================================================
def pn(lnq2, y, E):
    '''

    Parameters
    ----------
    lnq2 : float
        Natural log of Q^2 value, in GeV^2.
    y : float
        Integration variable.
    E : float
        Energy, in GeV.

    Returns
    -------
    float
        Photonuclear energy loss value.

    '''
    q2 = np.exp(lnq2)

    x = q2/(2*m_p*E*y) # x*y*S = Q^2; where S = 2*m_p*E

    def f2_allm(q2,y,E): # ALLM parameterization for f2.
        '''

        Parameters
        ----------
        q2 : float
            Q^2 value, in GeV^2.
        y : float
            Integration variable.
        E : float
            Energy, in GeV.

        Returns
        -------
        f2 : float
            F2 (EM) structure function for ALLM model.

        '''
        cp1, cp2, cp3 = 0.28067, 0.22291, 2.1979
        ap1, ap2, ap3 = -0.0808, -0.44812, 1.1709
        bp1, bp2, bp3 = 0.60243, 1.3754, 1.8439
        cr1, cr2, cr3 = 0.80107, 0.97307, 3.4942
        ar1, ar2, ar3 = 0.58400, 0.37888, 2.6063
        br1, br2, br3 = 0.10711, 1.9386, 0.49338
        m02, mr2, mp2 = 0.31985, 0.15052, 49.457 #GeV^2
        lambda2 = 0.06527 #GeV^2
        q02 = 0.46017 #GeV^2

        argt = np.log((q2+q02)/lambda2)/np.log(q02/lambda2)
        t = np.log(argt)
        cr = cr1 + cr2*(t**cr3)
        ar = ar1 + ar2*(t**ar3)
        cp = cp1 + (cp1-cp2)*((1/(1 + t**cp3)) - 1)
        ap = ap1 + (ap1-ap2)*((1/(1 + t**ap3)) - 1)
        br = br1**2 + (br2**2)*(t**br3)
        bp = bp1**2 + (bp2**2)*(t**bp3)

        p = 1 - 1.85*x + 2.45*x**2 - 2.35*x**3 + x**4
        W2 = q2*(1/x - 1) + m_p**2

        xp = (q2 + mp2)/(q2 + mp2 + W2 - m_p**2)
        xr = (q2 + mr2)/(q2 + mr2 + W2 - m_p**2)
        if xr<0:
            xr = 0.0

        f_2p = cp*(xp**ap)*(1 - x)**bp
        f_2r = cr*(xr**ar)*(1 - x)**br

        f2p = q2/(q2 + m02) * (f_2p + f_2r) # proton

        if (x>0 and x<0.0014):
            # f2 = A/2 * A**(-0.1)*(1 + p)*f2p
            f2 = A**(-0.1)*(z + (A-z)*p)*f2p # changed 17/12/2020
        elif (x>0.0014 and x<0.04):
            # f2 = A/2 * A**(0.069 * np.log10(x) + 0.097)*(1 + p)*f2p
            # f2 = A**(-0.1) * A**(0.069 * np.log10(x) + 0.097)*(z + (A-z)*p)*f2p  # changed 17/12/2020
            f2 = A**(0.069 * np.log10(x) + 0.097)*(z + (A-z)*p)*f2p  # changed 25/3/2021
        else:
            # f2 = A/2*(1 + p)*f2p
            f2 = (z + (A-z)*p)*f2p  # changed 17/12/2020

        return f2

    def f2_bdhm(q2,y,E): # BDHM parameterization for f2.
        '''

        Parameters
        ----------
        q2 : float
            Q^2 value, in GeV^2.
        y : float
            Integration variable.
        E : float
            Energy, in GeV.

        Returns
        -------
        f2 : float
            F2 (EM) structure function for BDHM model.

        '''
        a_0 = 8.205e-4
        a_1 = -5.148e-2
        a_2 = -4.725e-3
        b_0 = 2.217e-3
        b_1= 1.244e-2
        b_2 = 5.958e-4
        c_0 = 0.255
        c_1 = 1.475e-1
        n = 11.49
        lambda2 = 2.430 # not lambda^2! Python has a reserved keyword for lambda
        mu2 = 2.82 # GeV^2
        M2 = 0.753 # GeV^2
        # two_m_nu = q2/x


        if x > 0.1: # x<0.1 constraint; see in paper
            f2p = 0
        else:

            A_A = a_0 + a_1*np.log(1 + q2/mu2) + a_2*np.log(1+q2/mu2)**2
            B = b_0 + b_1*np.log(1+q2/mu2) + b_2*np.log(1+q2/mu2)**2
            C = c_0 + c_1*np.log(1+q2/mu2)
            D = (q2*(q2 + lambda2*M2))/(q2 + M2)**2

            # f2 = D*(1 - q2/two_m_nu)**n * (C + A*np.log(two_m_nu/(q2+mu2)) + B*np.log(two_m_nu/(q2+mu2))**2)

            f2p = D*(1 - x)**n * (C + A_A*np.log(1/x * q2/(q2 + mu2)) + B*np.log(1/x * q2/(q2 + mu2))**2)

        p = 1 - 1.85*x + 2.45*x**2 - 2.35*x**3 + x**4 # correction since neutrons != protons

        if (x>0 and x<0.0014):
            # f2 = A/2 * A**(-0.1)*(1 + p)*f2p
            f2 = A**(-0.1)*(z + (A-z)*p)*f2p # changed 17/12/2020
        elif (x>0.0014 and x<0.04):
            # f2 = A/2 * A**(0.069 * np.log10(x) + 0.097)*(1 + p)*f2p
            # f2 = A**(-0.1) * A**(0.069 * np.log10(x) + 0.097)*(z + (A-z)*p)*f2p  # changed 17/12/2020
            f2 = A**(0.069 * np.log10(x) + 0.097)*(z + (A-z)*p)*f2p  # changed 25/3/2021
        else:
            # f2 = A/2*(1 + p)*f2p
            f2 = (z + (A-z)*p)*f2p  # changed 17/12/2020

        return f2

    def f2_ckmt(q2,y,E): # CKMT parameterization for f2.
        '''

        Parameters
        ----------
        q2 : float
            Q^2 value, in GeV^2.
        y : float
            Integration variable.
        E : float
            Energy, in GeV.

        Returns
        -------
        f2 : float
            F2 (EM) structure function for CKMT model.

        '''
        A_A = 0.1502 # A is reserved for atomic mass number
        delta_0 = 0.07684
        B = 1.2064
        alpha_R = 0.4150
        # f = 0.15
        a = 0.2631 # GeV^2
        d = 1.1170 # GeV^2
        b = 0.6452 # GeV^2
        c = 3.5489 # GeV^2

        n_q2 = 3/2 * (1 + q2/(q2+c))
        delta_q2 = delta_0 * (1 + (2*q2)/(q2+d))

        f2_sea = A_A*x**(-delta_q2) * (1-x)**(n_q2+4) * (q2/(q2+a))**(1+delta_q2)
        # f2_val = B*x**(1-alpha_R) * (1-x)**n_q2 * (q2/(q2+b))**alpha_R * (1+f*(1-x))
        f2_val = B*x**(1-alpha_R) * (1-x)**n_q2 * (q2/(q2+b))**alpha_R

        f2p = f2_sea + f2_val

        p = 1 - 1.85*x + 2.45*x**2 - 2.35*x**3 + x**4 # correction since neutrons != protons

        if (x>0 and x<0.0014):
            # f2 = A/2 * A**(-0.1)*(1 + p)*f2p
            f2 = A**(-0.1)*(z + (A-z)*p)*f2p # changed 17/12/2020
        elif (x>0.0014 and x<0.04):
            # f2 = A/2 * A**(0.069 * np.log10(x) + 0.097)*(1 + p)*f2p
            # f2 = A**(-0.1) * A**(0.069 * np.log10(x) + 0.097)*(z + (A-z)*p)*f2p  # changed 17/12/2020
            f2 = A**(0.069 * np.log10(x) + 0.097)*(z + (A-z)*p)*f2p  # changed 25/3/2021
        else:
            # f2 = A/2*(1 + p)*f2p
            f2 = (z + (A-z)*p)*f2p  # changed 17/12/2020

        return f2

        # def f2_custom(q2,y,E):

            # your parameterization of F2. Even though the F2 doesn't depend on y & E, you will have to pass them as dummy variables in f2_custom(q2,y,E) as pass through values.

            # return f2

    R = 0
# only one power of q2 in denominator, since this is integral over log(q2)
    dsigmadq2dx = (4*np.pi*alpha_fs**2)/q2 * f2_allm(q2,y,E)/x * (1 - y - (m_p*x*y)/(2*E) + (1 - (2*m_le**2)/q2)*y**2*((1 + 4*m_p**2*x**2/q2)/(2*(1+R)))) # NB: f2_allm is used here by default

    y_dsigmadq2dy = y * dsigmadq2dx * q2/(2*m_p*E*y**2) * 0.389e-27 # changing differential variable from x to y

    return (N_A/A) * y_dsigmadq2dy

def cs_pn(lnq2, y, E):
    '''

    Parameters
    ----------
    lnq2 : float
        Natural log of Q^2 value, in GeV^2.
    y : float
        Integration variable.
    E : float
        Energy, in GeV.

    Returns
    -------
    float
        Photonuclear cross-section value.

    '''
    return pn(lnq2, y, E)/y

def pn_q2_cut(y, E):
    '''

    Parameters
    ----------
    y : float
        Integration variable.
    E : float
        Energy, in GeV.

    Returns
    -------
    float
        Photonuclear Q^2 integration limit, for beta_cut.

    '''
    return [np.log((m_le**2 * y**2)/(1 - y)), np.log(2*m_p*E*y - ((m_p+m_pi)**2-m_p**2))] # [low, high]

def pn_q2_tot(y, E):
    '''

    Parameters
    ----------
    y : float
        Integration variable.
    E : float
        Energy, in GeV.

    Returns
    -------
    float
        Photonuclear Q^2 integration limit, for beta_total.

    '''
    return [np.log((m_le**2 * y**2)/(1 - y)), np.log(2*m_p*E*y - ((m_p+m_pi)**2-m_p**2))] # [low, high]


def pn_y_cut(E):
    '''

    Parameters
    ----------
    y : float
        Integration variable.
    E : float
        Energy, in GeV.

    Returns
    -------
    float
        Photonuclear y integration limit, for beta_cut.

    '''
    return [((m_p+m_pi)**2-m_p**2)/(2*m_p*E), 1e-3] # [low, high]

def pn_y_tot(E):
    '''

    Parameters
    ----------
    y : float
        Integration variable.
    E : float
        Energy, in GeV.

    Returns
    -------
    float
        Photonuclear y integration limit, for beta_total.

    '''
    return [((m_p + m_pi)**2 - m_p**2)/(2*m_p*E), (1 - m_le/E)] # [low, high]


#==================================================
#     Add Lookup Table Entries For Alpha
# =============================================================================

def calc_alpha(lepton, material):
    '''

    Parameters
    ----------
    lepton : str
        Type of lepton. Can be tau or muon.
    material : str
        Material for lepton propagation.

    Returns
    -------
    NONE
        Calculates & adds lepton ionization energy loss lookup entries in lookup_tables.h5.

    '''
    alpha_arr = np.asarray([alpha_i(i) if i>m_le else 0 for i in E_lep])
    Data.add_alpha(alpha_arr, lepton, material)
    return 'Problem in calc_alpha function'

# =============================================================================
#     Add Lookup Table Entries For Beta Cut/Total
# =============================================================================

def calc_beta(lepton, material, process, beta_type): # energy_loss_type can be 'brem', 'pair', 'pn_bb' or 'pn_*'
    '''

    Parameters
    ----------
    lepton : str
        Type of lepton. Can be tau or muon.
    material : str
        Material for lepton propagation.
    process : str
        Lepton energy loss (non-ionization) process/model. Can be brem, pair, pn_bb or pn_*
    beta_type : str
        Type of beta to be used for calculation. Can be cut or total.

    Returns
    -------
    None
        Calculates & adds lepton (non-ionization) energy loss lookup entries in lookup_tables.h5.

    '''

    p = Pool(mp.cpu_count()) # use all available cores

    if process == 'brem' and beta_type == 'cut':
        brem_arr = np.asarray([p.map(integrator,[[brem, 0, 1e-3, i]])[0][0] for i in E_lep])
        brem_arr = np.asarray([rep(i) for i in brem_arr])
        Data.add_beta(brem_arr, lepton, material, 'brem', beta_type)

    elif process == 'brem' and beta_type == 'total':
        brem_arr = np.asarray([p.map(integrator,[[brem, 1e-7, 1, i]])[0][0] for i in E_lep])
        brem_arr = np.asarray([rep(i) for i in brem_arr])
        Data.add_beta(brem_arr, lepton, material, 'brem', beta_type)

    elif process == 'pair' and beta_type == 'cut':
        pair_arr = np.asarray([p.map(nquad_integrator,[[pair, pair_rho_cut, pair_y_cut, i]])[0][0] for i in E_lep])
        pair_arr = np.asarray([rep(i) for i in pair_arr])
        Data.add_beta(pair_arr, lepton, material, 'pair', beta_type)

    elif process == 'pair' and beta_type == 'total':
        pair_arr = np.asarray([p.map(nquad_integrator,[[pair, pair_rho_tot, pair_y_tot, i]])[0][0] for i in E_lep])
        pair_arr = np.asarray([rep(i) for i in pair_arr])
        Data.add_beta(pair_arr, lepton, material, 'pair', beta_type)

    elif process == 'pn_bb' and beta_type == 'cut':
        pn_bb_arr = np.asarray([p.map(integrator,[[pn_bb, 1e-5, 1e-3, i]])[0][0] for i in E_lep])
        pn_bb_arr = np.asarray([rep(i) for i in pn_bb_arr])
        Data.add_beta(pn_bb_arr, lepton, material, 'pn_bb', beta_type)

    elif process == 'pn_bb' and beta_type == 'total':
        pn_bb_arr = np.asarray([p.map(integrator,[[pn_bb, 1e-5, 1, i]])[0][0] for i in E_lep])
        pn_bb_arr = np.asarray([rep(i) for i in pn_bb_arr])
        Data.add_beta(pn_bb_arr, lepton, material, 'pn_bb', beta_type)

    elif beta_type == 'cut': # pn_*
        pn_arr = np.asarray([p.map(nquad_integrator,[[pn, pn_q2_cut, pn_y_cut, i]])[0][0] for i in E_lep])
        pn_arr = np.asarray([rep(i) for i in pn_arr])
        Data.add_beta(pn_arr, lepton, material, process, beta_type)

    elif beta_type == 'total': # pn_*
        pn_arr = np.asarray([p.map(nquad_integrator,[[pn, pn_q2_tot, pn_y_tot, i]])[0][0] for i in E_lep])
        pn_arr = np.asarray([rep(i) for i in pn_arr])
        Data.add_beta(pn_arr, lepton, material, process, beta_type)

    p.close()

    return None

# =============================================================================
# Calculate Integrated & Absolute Lepton Cross-Sections
# =============================================================================

def calc_ixc(lepton, material, process):
    '''

    Parameters
    ----------
    lepton : str
        Type of lepton. Can be tau or muon.
    material : str
        Material for lepton propagation.
    process : str
        Lepton energy loss (non-ionization) process/model. Can be brem, pair, pn_bb or pn_*

    Returns
    -------
    NONE
        Calculates & adds lepton (non-ionization) energy loss integrated (CDFs) & absolute cross-section lookup entries in lookup_tables.h5.

    '''

    yvals = 10**(-np.linspace(0.1,3,num=30)) # The integrated cross-section values should go from y = 10^(-0.1) to y = 10^(-3). This is a convention we chose to adopt.

    ixc_dict = {}
    xc_arr = []

    for E in E_lep:

        ixc_dict.update({E:{0:0}}) # pad with 0 at the beginning because it is the CDF after all

        if process == 'brem':
            for i in range(1,31):
                brem_val = integrator([cs_brem,yvals[i-1],brem_bb_high(E),E])[0]
                ixc_dict[E].update({i:rep(brem_val)})

                if i==30: # update absolute cross-section values
                    xc_arr.append(rep(brem_val))

        elif process == 'pair':
            def pair_y_high_ics(E):
                return 1 - (3*m_le/(4*E)) * np.sqrt(np.e) * z**(1/3)

            def pair_rho_ics(y, E):
                return [-(1 - (6*m_le**2)/(E**2 * (1-y))) * np.sqrt(1 - (4 * m_e)/(E * y)),(1 - (6*m_le**2)/(E**2 * (1-y))) * np.sqrt(1 - (4 * m_e)/(E * y))]

            pair_int =  short_int_nquad(cs_pair, pair_y_high_ics,pair_rho_ics, yvals, E)

            for i in range(1,31):
                ixc_dict[E].update({i:rep(pair_int[i-1])})

                if i==30: # update absolute cross-section values
                    xc_arr.append(rep(pair_int[-1]))

        elif process == 'pn_bb':
            for i in range(1,31):
                pn_bb_val = integrator([cs_pn_bb,yvals[i-1],1,E])[0]
                ixc_dict[E].update({i:rep(pn_bb_val)})

                if i==30: # update absolute cross-section values
                    xc_arr.append(rep(pn_bb_val))

        else: # pn_*
            def pn_q2(y, E):
                return [np.log((m_le**2 * y**2)/(1 - y)),np.log(2*m_p*E*y - ((m_p+m_pi)**2-m_p**2))]

            def pn_y_high_ics(E):
                return (1 - m_le/E)

            pn_int = short_int_nquad(cs_pn, pn_y_high_ics, pn_q2, yvals, E)
            # pn_int = quadpy.quad(cs_pn, pn_y_high_ics, pn_q2, yvals, E)[0]

            for i in range(1,31):
                ixc_dict[E].update({i:rep(pn_int[i-1])})

                if i==30: # update absolute cross-section values
                    xc_arr.append(rep(pn_int[-1]))

        for j in range(0,31): # to normalize ixc's
            if ixc_dict[E][30]!=0:
                ixc_dict[E].update({j:ixc_dict[E][j]/ixc_dict[E][30]})

    ixc_dframe = pd.DataFrame.from_dict(ixc_dict,orient='index').transpose()

    Data.add_ixc(lepton, ixc_dframe, process, material=material)
    Data.add_xc(lepton, xc_arr, process, material=material)
    return None

# =============================================================================
#     Plot Alpha
# =============================================================================
def plot_alpha():
    plt.figure(1)
    alpha = Data.get_alpha(lepton,material)['alpha']

    ener=np.genfromtxt('./fortran/muontables/alpha_he_muon_rock.dat',usecols=0)
    alpha_hls=np.genfromtxt('./fortran/tautables/alpha_he_tau_rock.dat',usecols=1)
    # alpha_hls=np.genfromtxt('./fortran/muontables/alpha_he_muon_rock.dat',usecols=1)

    plt.loglog(E_lep, alpha,ls=':',color='k',label='Alpha (NuSamProp)')
    plt.loglog(ener, alpha_hls,color='k',label='Alpha (HLS)')

    plt.legend(loc='best')
    plt.xlabel(r"$E_{\tau}$ [GeV]") # or E_{\tau}
    plt.ylabel(r"$\alpha$")
    plt.xlim(100, 1e9)
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth='0.2', color='black')
    # plt.savefig('Alpha.png', format='png', dpi = 300)

    return None
# =============================================================================
#   Plot Beta - Cut
# =============================================================================
def plot_c(model):
    plt.figure(2)

    brem = Data.get_beta(particle=lepton,material=material,pn_model=model,beta_type='cut')['brem']
    pair = Data.get_beta(particle=lepton,material=material,pn_model=model,beta_type='cut')['pair']
    pn = Data.get_beta(particle=lepton,material=material,pn_model=model,beta_type='cut')['pn']

    ener=np.genfromtxt('./fortran/muontables/betac_he_muon_rock.dat',usecols=0)

    brem_hls = np.genfromtxt('./fortran/tautables_2/betac_he_tau_rock.dat',usecols=1)
    # brem_hls = np.genfromtxt('./fortran/muontables/betac_he_muon_rock.dat',usecols=1)

    pair_hls = np.genfromtxt('./fortran/tautables_2/betac_he_tau_rock.dat',usecols=2)
    # pair_hls = np.genfromtxt('./fortran/muontables/betac_he_muon_rock.dat',usecols=2)

    pn_bb_hls = np.genfromtxt('./fortran/tautables_2/betac_he_tau_rock.dat',usecols=3)
    # pn_bb_hls = np.genfromtxt('./fortran/muontables/betac_he_muon_rock.dat',usecols=3)

    pn_hls = np.genfromtxt('./fortran/tautables_2/Allms/betac-tau-rock-highstat-2020.dat',usecols=3)
    # pn_hls = np.genfromtxt('./fortran/muontables/betac-muon-rock-allm.dat',usecols=3)


    plt.loglog(E_lep, [brem[i] for i in E_lep],ls=':',color='r', label = "Bremsstrahlung - NuSamProp")
    plt.loglog(E_lep, [pair[i] for i in E_lep],ls=':',color='b', label = "Pair Production - NuSamProp")
    # plt.loglog(E_lep, [pn_bb[i] for i in E_lep],ls=':',color='g', label = "Photonuclear (BB) - NuSamProp")
    plt.loglog(E_lep, [pn[i] for i in E_lep],ls=':',color='k', label = "Photonuclear - NuSamProp")

    plt.loglog(ener, brem_hls,color='r', label = "Bremsstrahlung - HLS")
    plt.loglog(ener, pair_hls,color='b', label = "Pair Production - HLS")
    plt.loglog(ener, pn_bb_hls,color='g', label = "Photonuclear (BB) - HLS")
    plt.loglog(ener, pn_hls,color='k', label = "Photonuclear (ALLM) - HLS")
    plt.title("Beta - Cut (beta)")
    plt.legend(loc='best')
    plt.xlabel(r"$E_{\tau}$ [GeV]") # or E_{\tau}
    plt.ylabel(r"$\beta_{std~rock}$ [cm$^2$/g]")
    plt.xlim(100, 1e9)
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth='0.2', color='black')
    # plt.savefig('Cut.png', format='png', dpi = 300)
    return None

# =============================================================================
#   Plot Beta - Total
# =============================================================================
def plot_t(model):
    plt.figure(3)

    brem_arr = Data.get_beta(particle=lepton,material=material,pn_model=model,beta_type='total')['brem']
    pair_arr = Data.get_beta(particle=lepton,material=material,pn_model=model,beta_type='total')['pair']
    pn_arr = Data.get_beta(particle=lepton,material=material,pn_model=model,beta_type='total')['pn']

    ener = np.genfromtxt('./fortran/tautables/fort.26',usecols=0)

    brem_hls = np.genfromtxt('./fortran/tautables/fort.26',usecols=2)

    pair_hls = np.genfromtxt('./fortran/tautables/fort.26',usecols=3)

    pn_bb_hls = np.genfromtxt('./fortran/tautables/fort.26',usecols=4)

    pn_hls = np.genfromtxt('./fortran/tautables/allms_betat.18',usecols=1)


    plt.loglog(E_lep, [brem_arr[i] for i in E_lep],ls=':',color='r', label = "Bremsstrahlung - NuSamProp")
    plt.loglog(E_lep, [pair_arr[i] for i in E_lep],ls=':',color='b', label = "Pair Production - NuSamProp")
    plt.loglog(E_lep, [pn_arr[i] for i in E_lep],ls=':',color='k', label = "Photonuclear - NuSamProp")

    plt.loglog(ener, brem_hls,color='r', label = "Bremsstrahlung - HLS")
    plt.loglog(ener, pair_hls,color='b', label = "Pair Production - HLS")
    plt.loglog(ener, pn_bb_hls,color='g', label = "Photonuclear (BB) - HLS")
    plt.loglog(ener, pn_hls,color='k', label = "Photonuclear (ALLM) - HLS")

    plt.title("Beta - Total (betat)")
    plt.legend(loc='best')
    plt.xlabel(r"$E_{\tau}$ [GeV]") # or E_{\tau}
    plt.ylabel(r"$\beta_{std~rock}$ [cm$^2$/g]")
    plt.xlim(100, 1e9)
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth='0.2', color='black')
    # plt.savefig('Total.png', format='png', dpi = 300)
    return None

# =============================================================================
# Plot CS
# =============================================================================
def plot_cs_old():
    plt.figure(4)

    brem_arr = Data.get_xc(lepton,'allm',material=material)['brem']
    pair_arr = Data.get_xc(lepton,'allm',material=material)['pair']
    pn_arr = Data.get_xc(lepton,'allm',material=material)['pn']

    # ener=np.genfromtxt('./fortran/muontables/allms_cs.17',usecols=0)
    ener=np.genfromtxt('./fortran/tautables/allms_cs.17',usecols=0)
    brem_hls = np.genfromtxt('./fortran/tautables/fort.17',usecols=1)
    # brem_hls = np.genfromtxt('./fortran/muontables/fort.17',usecols=1)
    pair_hls = np.genfromtxt('./fortran/tautables/fort.17',usecols=2)
    # pair_hls = np.genfromtxt('./fortran/muontables/fort.17',usecols=2)
    # pn_bb_hls = np.genfromtxt('./fortran/tautables/fort.17',usecols=3)
    # pn_bb_hls = np.genfromtxt('./fortran/muontables/fort.17',usecols=3)
    pn_hls = np.genfromtxt('./fortran/tautables/allms_cs.17',usecols=3)
    # pn_hls = np.genfromtxt('./fortran/muontables/allms_cs.17',usecols=3)

    models = ['brem', 'pair', 'nucl', 'allm']
    # xc_dict_hls = {}
    # for model in models:
    #     file = './fortran/muontables/%s_he_muon_rock.dat' % str(model)
    #     # file = './fortran/muontables/%s_he_muon_rock.dat' % str(model)
    #     data_dict = {}
    #     skp = 0
    #     for energy in E_lep:
    #         data_lst = sorted(np.concatenate((np.genfromtxt(file,usecols=(0,1,2,3,4,5),skip_header=skp,max_rows=5))),reverse=False) # NOTE: sigma is the last (30th) value (1e-3 -> 1) # units = cm^2
    #         data_dict.update({energy:data_lst[-1]})
    #         skp+=6
    #     xc_dict_hls.update({model:data_dict})

    # brem = xc_dict['brem']
    # pair = xc_dict['pair']
    # pn_bb = xc_dict['pn_bb']
    # pn = xc_dict['pn']

    # brem_hls = np.array(list(xc_dict_hls['brem'].values()))
    # pair_hls = np.array(list(xc_dict_hls['pair'].values()))
    # pn_bb_hls = np.array(list(xc_dict_hls['nucl'].values()))
    # pn_hls = np.array(list(xc_dict_hls['allm'].values()))

    plt.loglog(E_lep, [brem_arr[i] for i in E_lep],ls=':',color='r', label = "Bremsstrahlung - NuSamProp")
    plt.loglog(E_lep, [pair_arr[i] for i in E_lep],ls=':',color='b', label = "Pair Production - NuSamProp")
    # plt.loglog(E_lep, [pn_bb_arr[i] for i in E_lep],ls=':',color='g', label = "Photonuclear (BB) - NuSamProp")
    plt.loglog(E_lep, [pn_arr[i] for i in E_lep],ls=':',color='k', label = "Photonuclear - NuSamProp")

    plt.loglog(E_lep, brem_hls,color='r', label = "Bremsstrahlung - HLS")
    plt.loglog(E_lep, pair_hls,color='b', label = "Pair Production - HLS")
    # plt.loglog(E_lep, pn_bb_hls,color='g', label = "Photonuclear (BB) - HLS")
    plt.loglog(ener, pn_hls,color='k', label = "Photonuclear (ALLM) - HLS")

    plt.title("Cross-Section (CS)")
    plt.legend(loc='best')
    plt.xlabel(r"$E_{\tau}$ [GeV]") # or E_{\tau}
    plt.ylabel(r"$\sigma$")
    plt.xlim(100, 1e9)
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth='0.2', color='black')
    # plt.savefig('cs.png', format='png', dpi = 300)
    return None

# =============================================================================
# Plot cs NPP Models
# =============================================================================
def plot_cs(lepton, material):
    brem_arr = A/N_A * Data.get_xc(lepton, 'brem', material=material)
    pair_arr = A/N_A * Data.get_xc(lepton, 'pair', material=material)
    pn_bb = A/N_A * Data.get_xc(lepton, 'pn_bb', material=material)
    pn_allm = A/N_A * Data.get_xc(lepton, 'pn_allm', material=material)
    pn_bdhm = A/N_A * Data.get_xc(lepton, 'pn_bdhm', material=material)
    pn_ckmt = A/N_A * Data.get_xc(lepton, 'pn_ckmt', material=material)



    ener=np.genfromtxt('./fortran/tautables/allms_cs.17',usecols=0)
    # ener=np.genfromtxt('./fortran/tautables/fort.17',usecols=0)
    brem_hls = np.genfromtxt('./fortran/tautables/fort.17',usecols=1)
    # brem_hls = np.genfromtxt('./fortran/muontables/fort.17',usecols=1)
    pair_hls = np.genfromtxt('./fortran/tautables/fort.17',usecols=2)
    # pair_hls = np.genfromtxt('./fortran/muontables/fort.17',usecols=2)
    pn_bb_hls = np.genfromtxt('./fortran/tautables/fort.17',usecols=3)
    # pn_bb_hls = np.genfromtxt('./fortran/muontables/fort.17',usecols=3)
    pn_hls = np.genfromtxt('./fortran/tautables/allms_cs.17',usecols=3)
    # pn_hls = np.genfromtxt('./fortran/muontables/allms_cs.17',usecols=3)

    fig, ax = plt.subplots()

    # ax.loglog(E_lep, brem_arr,color='r', label = "Bremsstrahlung")
    # ax.loglog(E_lep, pair_arr,color='b', label = "Pair Production")
    ax.loglog(E_lep, pn_bb,color='g', label = "Photonuclear - BB")
    ax.loglog(E_lep, pn_allm,color='k', label = "Photonuclear - ALLM")
    ax.loglog(E_lep, pn_bdhm,color='y', label = "Photonuclear - BDHM")
    ax.loglog(E_lep, pn_ckmt,color='m', label = "Photonuclear - CKMT")

    # ax.loglog(ener, brem_hls, ls=':', color='r', label = "Bremsstrahlung - HLS")
    # ax.loglog(ener, pair_hls, ls=':', color='b', label = "Bremsstrahlung - HLS")
    # ax.loglog(ener, pn_bb_hls, ls=':', color='g', label = "Bremsstrahlung - HLS")
    # ax.loglog(ener, pn_hls, ls=':', color='k', label = "Bremsstrahlung - HLS")

    ax.set_title("Cross-Section (CS)")
    ax.legend(loc='best')
    ax.set_xlabel(r"$E_{\%s}$ [GeV]" % lepton) # or E_{\tau}
    ax.set_ylabel(r"$\sigma$ [cm$^2$]")
    ax.set_xlim(100, 1e9)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.grid(which='minor', linestyle=':', linewidth='0.2', color='black')
    ax.tick_params(axis='x', which='both', labelbottom = True, labeltop = True)
    ax.tick_params(axis='y', which='both', left = True, labelleft = True, labelright= True)
    # plt.savefig('cs.png', format='png', dpi = 300)

# =============================================================================
#   Plot Beta - Cut/Total
# =============================================================================
def plot_beta(lepton, material, beta_type):
    brem_arr = Data.get_beta(lepton,material,'brem',beta_type)
    pair_arr = Data.get_beta(lepton,material,'pair',beta_type)
    pn_bb_arr = Data.get_beta(lepton,material,'pn_bb',beta_type)
    pn_allm_arr = Data.get_beta(lepton,material,'pn_allm',beta_type)
    pn_bdhm_arr = Data.get_beta(lepton,material,'pn_bdhm',beta_type)
    pn_ckmt_arr = Data.get_beta(lepton,material,'pn_ckmt',beta_type)

    ener = np.genfromtxt('./fortran/tautables/fort.26',usecols=0)

    brem_hls = np.genfromtxt('./fortran/tautables/fort.26',usecols=2)

    pair_hls = np.genfromtxt('./fortran/tautables/fort.26',usecols=3)

    pn_bb_hls = np.genfromtxt('./fortran/tautables/fort.26',usecols=4)

    pn_hls = np.genfromtxt('./fortran/tautables/allms_betat.18',usecols=1)

    fig, ax = plt.subplots()

    # ax.loglog(E_lep, brem_arr, color='r', label = "Bremsstrahlung")
    # ax.loglog(E_lep, pair_arr, color='b', label = "Pair Production")
    ax.loglog(E_lep, pn_bb_arr, color='g', label = "Photonuclear - BB")
    ax.loglog(E_lep, pn_allm_arr, color='k', label = "Photonuclear - ALLM")
    ax.loglog(E_lep, pn_bdhm_arr, color='y', label = "Photonuclear - BDHM")
    ax.loglog(E_lep, pn_ckmt_arr, color='m', label = "Photonuclear - CKMT")

    # ax.loglog(ener, brem_hls,color='r', label = "Bremsstrahlung - HLS")
    # ax.loglog(ener, pair_hls,color='b', label = "Pair Production - HLS")
    # ax.loglog(ener, pn_bb_hls,color='g', label = "Photonuclear (BB) - HLS")
    # ax.loglog(ener, pn_hls,color='k', ls = ':', label = "Photonuclear (ALLM) - HLS")

    ax.set_title("Beta - %s" % beta_type)
    ax.legend(loc='best')
    if lepton == 'muon':lepton = 'mu'
    ax.set_xlabel(r"$E_{\%s}$ [GeV]" % lepton) # or E_{\tau}
    if material == 'rock':material = 'std~rock'
    if material == 'iso_water':material = 'iso~water'
    ax.set_ylabel(r"$\beta_{%s}$ [cm$^2$/g]" % material)
    ax.set_xlim(100, 1e9)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.grid(which='minor', linestyle=':', linewidth='0.2', color='black')
    ax.tick_params(axis='x', which='both', labelbottom = True, labeltop = True)
    ax.tick_params(axis='y', which='both', left = True, labelleft = True, labelright= True)
    # plt.savefig('Total.png', format='png', dpi = 300)
    return None


# =============================================================================
# Test
# =============================================================================
if __name__ == "__main__":
    start_time = time.time()

    leptons = ['tau']
    materials = ['rock']

    for lepton in leptons:
        for material in materials:

            if lepton=='tau':m_le=m_tau # mass of lepton
            else:m_le=m_mu  # mass of muon

            if material=='iso_water':z=7.0
            elif material=='water':z=6.6 # for NuTauSim Comparison
            elif material=='rock':z=11.0
            elif material=='iron':z=26.0
            else: z=float(input("Enter the atomic number of %s: " % material))

            if material=='iso_water':A=14.0
            elif material=='water':A=11.89 # for NuTauSim Comparison
            elif material=='rock':A=22.0
            elif material=='iron':A=55.84
            else: A=float(input("Enter the atomic mass of %s: " % material))

            # calc_beta(lepton, material, 'pn_ckmt', 'cut')
            # calc_beta(lepton, material, 'pn_ckmt', 'total')

            # calc_ixc(lepton, material, 'pn_allm')

    # plot_cs('tau', 'rock')
    # plot_beta('tau', 'rock', 'total')
    # plot_beta('tau', 'rock', 'cut')

    # plt_struc()

    end_time = time.time()
    print(f"It took {end_time-start_time:.2f} seconds to compute")