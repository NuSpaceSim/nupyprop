#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 3 14:17:24 2020

@author: sam
"""
# import matplotlib as mpl
# mpl.use('Agg') # for clusters

import data as Data
import cross_section as Cross_section

import numpy as np
import pandas as pd
import scipy.integrate as integrate
import scipy.constants as scc
import multiprocessing as mp
from multiprocessing import Pool
import matplotlib.pyplot as plt
# from numba import njit
import pickle

m_e = scc.physical_constants["electron mass energy equivalent in MeV"][0]*1e-3 # GeV
m_mu = scc.physical_constants["muon mass energy equivalent in MeV"][0]*1e-3 # GeV
m_tau = scc.physical_constants["tau mass energy equivalent in MeV"][0]*1e-3 # GeV
m_pi = 139.57018e-3 # pion mass in GeV
alpha = scc.fine_structure
le = 3.8616e-11 # electron Compton wavelength in cm

E_lep = Data.E_lep
m_p = Cross_section.m_p
N_A = Cross_section.N_A

def integrator(args):
    fun = args[0]
    low_lim = args[1]
    up_lim = args[2]
    arg = args[3]
    return integrate.quad(fun,low_lim,up_lim,args=(arg,))

def nquad_integrator(args):
    fun = args[0]
    lim_1 = args[1] # dep limit
    lim_2 = args[2] # indep limit
    arg = args[3]
    return integrate.nquad(fun, [lim_1,lim_2], args=(arg,))

def dbl_integrator(args):
    fun = args[0]
    indep_lim_low = args[1] # has to be a value and not a function, but in the form of a function! See pair and PN limits
    indep_lim_high = args[2] # has to be a value and not a function, but in the form of a function!
    dep_lim_low = args[3]
    dep_lim_high = args[4]
    var = args[5] # constant we want to loop over (for eg. energy)
    return integrate.dblquad(fun, indep_lim_low(var), indep_lim_high(var), dep_lim_low, dep_lim_high, args=(var,), epsabs=1e-12)[0]

def rep(val): # Remove non-physical values
    if (val<0) or (np.isnan(val)==True):
        return 0
    else:
        return val
    return 'Problem in rep function'

def short_int_nquad(fn, indep_fn_lim_high, dep_fn_lim, arr, arg):
    ans = []
    result = integrate.nquad(fn, [dep_fn_lim,[arr[0],indep_fn_lim_high(arg)]],args=(arg,))[0]
    ans.append(result)
    for i in range(len(arr)-1):
        result+= integrate.nquad(fn, [dep_fn_lim,[arr[i+1],arr[i]]],args=(arg,))[0]
        ans.append(result)
    return np.asarray(ans)

def short_int_dblquad(fn, indep_fn_lim_high, dep_fn_lim_low, dep_fn_lim_high, arr, arg): #indep_fn_lim_high has to return a value and not a function!
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

    dEdX = alpha**2 * 2*np.pi*N_A*le**2*(z*m_e)/(A*beta**2)*(np.log((2*m_e*beta**2*gamma**2*E_m)/I**2) - 2*beta**2 + E_m**2/(4*E**2) - delta)

    return dEdX


# =============================================================================
#      Bremmstrahlung Energy Loss
# =============================================================================
def brem(y, E): # eq. A5
    y_max = 1 - (3*m_le/(4*E)) * np.sqrt(np.e) * z**(1/3)
    if y>y_max:
        y_dsig_dy = 0
    else:
        delta = (m_le**2*y)/(2*E*(1-y))
        if z > 10:
            phi = np.log((2/3 * 189/m_e * m_le * z**(-2/3))/(1 + 189 * np.sqrt(np.e)/m_e * delta * z**(-1/3)))
        else:
            phi = np.log((189 * m_le/m_e * z**(-1/3))/(1 + 189* np.sqrt(np.e)/m_e * delta * z**(-1/3)))

        brem = alpha**3 * (4*z*(z+1)*le**2) * (m_e**2/m_le**2) * (1/y) * (4/3 - (4*y)/3 + y**2)*phi

        y_dsig_dy = brem * y

    return (N_A/A) * y_dsig_dy

def cs_brem(y, E):
    return brem(y,E)/y


def brem_bb_high(E):
    return 1 - (3*m_le/(4*E)) * np.sqrt(np.e) * z**(1/3)
# =============================================================================
#     Pair Production Energy Loss
# =============================================================================

def pair(rho, y, E): # eq. A9
    R = 189
    beta = y**2/(2*(1-y))
    xi = ((m_le * y)/(2*m_e))**2 * (1 - rho**2)/(1 - y)

    y_l = (4 + rho**2 + 3*beta*(1 + rho**2))/((1 + rho**2) * (3/2 + 2*beta) * np.log(3 + xi) + 1 - (3/2)*rho**2)

    y_e = (5 - rho**2 + 4*beta*(1 + rho**2))/(2*(1 + 3*beta) * np.log(3 + 1/xi) - rho**2 - 2*beta*(2 - rho**2))

    L_l = np.log((R * z**(-2/3) * 2/3 * m_le/m_e)/(1 + (2*m_e*np.sqrt(np.e)*R * z**(-1/3)) * (1 + xi) * (1 + y_l)/(E * y * (1 - rho**2))))

    L_e = np.log((R * z**(-1/3) * np.sqrt((1 + xi) * (1 + y_e)))/(1 + (2*m_e*np.sqrt(np.e)*R * z**(-1/3) * (1 + xi) * (1 + y_e))/(E * y * (1 - rho**2)))) - 1/2 * np.log(1 + (3/2 * m_e/m_le * z**(1/3))**2 * (1 + xi) * (1 + y_e))

    phi_l = (((1 + rho**2) * (1 + (3/2)*beta) - 1/xi * (1 + 2*beta) * (1 - rho**2)) * np.log(1 + xi) + (xi * (1 - rho**2 - beta))/(1 + xi) + (1 + 2*beta) * (1 - rho**2)) * L_l

    phi_e = (((2 + rho**2) * (1 + beta) + xi*(3 + rho**2)) * np.log(1 + 1/xi) + (1 - rho**2 - beta)/(1 + xi) - (3 + rho**2)) * L_e

    y_d2sigma_dydrho = y * alpha**4 * 2/(3 * np.pi) *z*(z+1) *(le)**2 * (1-y)/y * (phi_e + (m_e**2/m_le**2) * phi_l)

    return (N_A/A) * y_d2sigma_dydrho

def cs_pair(rho, y, E):
    return pair(rho, y, E)/y

def pair_rho_cont(y, E):
    return [-(1 - (6*m_le**2)/(E**2 * (1-y))) * np.sqrt(1 - (4 * m_e)/(E * y)), (1 - (6*m_le**2)/(E**2 * (1-y))) * np.sqrt(1 - (4 * m_e)/(E * y))] # [low, high]

def pair_rho_tot(y, E):
    return [-(1 - (6*m_le**2)/(E**2 * (1-y))) * np.sqrt(1 - (4 * m_e)/(E * y)), (1 - (6*m_le**2)/(E**2 * (1-y))) * np.sqrt(1 - (4 * m_e)/(E * y))] # [low, high]

def pair_rho_xc(y, E):
    return [-(1 - (6*m_le**2)/(E**2 * (1-y))) * np.sqrt(1 - (4 * m_e)/(E * y)), (1 - (6*m_le**2)/(E**2 * (1-y))) * np.sqrt(1 - (4 * m_e)/(E * y))] # [low, high]

def pair_y_cont(E):
    return [(4*m_e)/E, 1e-3] # [low, high]

def pair_y_tot(E):
    return [(4*m_e)/E, 1 - (3*m_le/(4*E)) * np.sqrt(np.e) * z**(1/3)] # [low, high]

def pair_y_xc(E):
    return [1e-3, 1 - (3*m_le/(4*E)) * np.sqrt(np.e) * z**(1/3)] # [low, high]

# =============================================================================
#     Photonuclear Energy Loss (Bezrukov-Bugaev)
# =============================================================================
def pn_bb(y, E): # eq. A12
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
    y_dsig_dy = y * alpha/(2*np.pi) * A * y * (sigma* (t1 + t2 + t3))*1e-30
    return (N_A/A) * y_dsig_dy

def cs_pn_bb(y, E):
    return pn_bb(y, E)/y

# =============================================================================
#     Photonuclear Energy Loss (ALLM/Custom)
# =============================================================================
def pn(lnq2, y, E):
    q2 = np.exp(lnq2)

    x = q2/(2*m_p*E*y)

    def f2_allm(q2,y,E): # ALLM parameterization for f2. Change this function for your own parameterization of f2.
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

        f2p = q2/(q2 + m02) * (f_2p + f_2r)

        if (x>0 and x<0.0014):
            # f2 = A/2 * A**(-0.1)*(1 + p)*f2p
            f2 = A**(-0.1)*(z + (A-z)*p)*f2p # changed 17/12/2020
        elif (x>0.014 and x<0.04):
            # f2 = A/2 * A**(0.069 * np.log10(x) + 0.097)*(1 + p)*f2p
            f2 = A**(-0.1) * A**(0.069 * np.log10(x) + 0.097)*(z + (A-z)*p)*f2p  # changed 17/12/2020
        else:
            # f2 = A/2*(1 + p)*f2p
            f2 = (z + (A-z)*p)*f2p  # changed 17/12/2020

        return f2

    R = 0
# only one power of q2 in denominator, since this is integral over log(q2)
    dsigmadq2dx = (4*np.pi*alpha**2)/q2 * f2_allm(q2,y,E)/x * (1 - y - (m_p*x*y)/(2*E) + (1 - (2*m_le**2)/q2)*y**2*((1 + 4*m_p**2*x**2/q2)/(2*(1+R))))

    y_dsigmadq2dy = y * dsigmadq2dx * q2/(2*m_p*E*y**2) * 0.389e-27 # changing differential variable from x to y

    return (N_A/A) * y_dsigmadq2dy

def cs_pn(lnq2, y, E):
    return pn(lnq2, y, E)/y

def pn_q2_cont(y, E):
    return [np.log((m_le**2 * y**2)/(1 - y)), np.log(2*m_p*E*y - ((m_p+m_pi)**2-m_p**2))] # [low, high]

def pn_q2_tot(y, E):
    return [np.log((m_le**2 * y**2)/(1 - y)), np.log(2*m_p*E*y - ((m_p+m_pi)**2-m_p**2))] # [low, high]

def pn_q2_xc(y, E):
    return [np.log((m_le**2 * y**2)/(1 - y)), np.log(2*m_p*E*y - ((m_p+m_pi)**2-m_p**2))] # [low, high]

def pn_y_cont(E):
    return [((m_p+m_pi)**2-m_p**2)/(2*m_p*E), 1e-3] # [low, high]

def pn_y_tot(E):
    return [((m_p + m_pi)**2 - m_p**2)/(2*m_p*E), (1 - m_le/E)] # [low, high]

def pn_y_xc(E):
    return [1e-3, (1-m_le/E)] # [low, high]

#==================================================
#     Add Lookup Table Entries For Alpha
# =============================================================================

def calc_alpha():
    alpha_arr = np.asarray([alpha_i(i) if i>m_le else 0 for i in E_lep])
    Data.add_alpha(alpha_arr, particle=lepton, material=material)
    return 'Problem in calc_alpha function'

# =============================================================================
#     Add Lookup Table Entries For Beta (Continuous)
# =============================================================================

def calc_beta_continuous():

    p = Pool(mp.cpu_count()) # use all available cores

    brem_arr = np.asarray([p.map(integrator,[[brem, 0, 1e-3, i]])[0][0] for i in E_lep]) # limits 0->1e-3
    brem_arr = np.asarray([rep(i) for i in brem_arr])

    pair_arr = np.asarray([p.map(nquad_integrator,[[pair, pair_rho_cont, pair_y_cont, i]])[0][0] for i in E_lep])
    pair_arr = np.asarray([rep(i) for i in pair_arr])

    pn_bb_arr = np.asarray([p.map(integrator,[[pn_bb, 1e-5, 1e-3, i]])[0][0] for i in E_lep]) # limits 1e-5->1e-3
    pn_bb_arr = np.asarray([rep(i) for i in pn_bb_arr])

    pn_arr = np.asarray([p.map(nquad_integrator,[[pn, pn_q2_cont, pn_y_cont, i]])[0][0] for i in E_lep])
    pn_arr = np.asarray([rep(i) for i in pn_arr])

    p.close()

    beta_dict = {'energy':E_lep,'brem':brem_arr,'pair':pair_arr,'pn_bb':pn_bb_arr,'pn':pn_arr}
    Data.add_beta(beta_dict, particle=lepton, material=material, beta_type='continuous')

    return 'Problem in calc_beta_continuous function'

# =============================================================================
# Beta Total
# =============================================================================

def calc_beta_total():

    p = Pool(mp.cpu_count()) # use all available cores

    brem_arr = np.asarray([p.map(integrator,[[brem, 1e-7, 1, i]])[0][0] for i in E_lep]) # limits 1e-7->1
    brem_arr = np.asarray([rep(i) for i in brem_arr])

    pair_arr = np.asarray([p.map(nquad_integrator,[[pair, pair_rho_tot, pair_y_tot, i]])[0][0] for i in E_lep])
    pair_arr = np.asarray([rep(i) for i in pair_arr])

    pn_bb_arr = np.asarray([p.map(integrator,[[pn_bb, 1e-5, 1, i]])[0][0] for i in E_lep]) # limits 1e-5->1
    pn_bb_arr = np.asarray([rep(i) for i in pn_bb_arr])

    pn_arr = np.asarray([p.map(nquad_integrator,[[pn, pn_q2_tot, pn_y_tot, i]])[0][0] for i in E_lep])
    pn_arr = np.asarray([rep(i) for i in pn_arr])

    p.close()

    beta_dict = {'energy':E_lep,'brem':brem_arr,'pair':pair_arr,'pn_bb':pn_bb_arr,'pn':pn_arr}
    Data.add_beta(beta_dict, particle=lepton, material=material, beta_type='total')

    return 'Problem in calc_beta_total function'

# # =============================================================================
# #     Calculate CS
# # =============================================================================

def calc_xc():
    p = Pool(mp.cpu_count()) # use all available cores

    brem_arr = np.asarray([p.map(integrator,[[cs_brem, 1e-3, brem_bb_high(i), i]])[0][0] for i in E_lep]) # limits 1e-7->brem_upper
    brem_cs = np.asarray([rep(i) for i in brem_arr])

    pair_arr = np.asarray([p.map(nquad_integrator,[[cs_pair, pair_rho_xc, pair_y_xc, i]])[0][0] for i in E_lep])
    pair_cs = np.asarray([rep(i) for i in pair_arr])

    pn_bb_arr = np.asarray([p.map(integrator,[[cs_pn_bb, 1e-3, 1, i]])[0][0] for i in E_lep]) # limits 1e-3->1
    pn_bb_cs = np.asarray([rep(i) for i in pn_bb_arr])

    pn_arr = np.asarray([p.map(nquad_integrator,[[cs_pn, pn_q2_xc, pn_y_xc, i]])[0][0] for i in E_lep])
    pn_cs = np.asarray([rep(i) for i in pn_arr])

    p.close()

    with open('brem_cs_%s_%s.npy' % (lepton,material), 'wb') as f:
        np.save(f, brem_cs)

    with open('pair_cs_%s_%s.npy' % (lepton,material), 'wb') as f:
        np.save(f, pair_cs)

    with open('pn_bb_cs_%s_%s.npy' % (lepton,material), 'wb') as f:
        np.save(f, pn_bb_cs)

    with open('pn_cs_%s_%s.npy' % (lepton,material), 'wb') as f:
        np.save(f, pn_cs)

    xc_dict = {'brem':brem_cs, 'pair':pair_cs, 'pn_bb':pn_bb_cs, 'pn':pn_cs}

    Data.add_xc(str(lepton), xc_dict, material = material)

    return 'Problem in calc_xc function'

# =============================================================================
#     Calculate Integrated Cross-Sections
# =============================================================================

def calc_ixc():
    yvals = 10**(-np.linspace(0.1,3,num=30)) # The integrated cross-section values should go from y = 10^(-0.1) to y = 10^(-3). This is a convention we chose to adopt.
    brem_dict = {}
    pair_dict = {}
    pn_bb_dict = {}
    pn_dict = {}

    brem_cs = []
    pair_cs = []
    pn_bb_cs = []
    pn_cs = []

    for E in E_lep:

        def pair_y_high_ics(E):
            return 1 - (3*m_le/(4*E)) * np.sqrt(np.e) * z**(1/3)

        def pair_rho_ics(y, E):
            return [-(1 - (6*m_le**2)/(E**2 * (1-y))) * np.sqrt(1 - (4 * m_e)/(E * y)),(1 - (6*m_le**2)/(E**2 * (1-y))) * np.sqrt(1 - (4 * m_e)/(E * y))]

        def pn_q2(y, E):
            return [np.log((m_le**2 * y**2)/(1 - y)),np.log(2*m_p*E*y - ((m_p+m_pi)**2-m_p**2))]

        def pn_y_high_ics(E):
            return (1 - m_le/E)

        brem_dict.update({E:{0:0}}) # pad with 0 at the beginning because it is the CDF after all
        pair_dict.update({E:{0:0}}) # pad with 0 at the beginning because it is the CDF after all
        pn_bb_dict.update({E:{0:0}}) # pad with 0 at the beginning because it is the CDF after all
        pn_dict.update({E:{0:0}}) # pad with 0 at the beginning because it is the CDF after all

        pair_int =  short_int_nquad(cs_pair, pair_y_high_ics,pair_rho_ics, yvals, E)

        pn_int = short_int_nquad(cs_pn, pn_y_high_ics, pn_q2, yvals, E)

        for i in range(1,31):
            brem_val = integrator([cs_brem,yvals[i-1],brem_bb_high(E),E])[0]
            brem_dict[E].update({i:rep(brem_val)})

            pair_dict[E].update({i:rep(pair_int[i-1])})

            pn_bb_val = integrator([cs_pn_bb,yvals[i-1],1,E])[0]
            pn_bb_dict[E].update({i:rep(pn_bb_val)})

            pn_dict[E].update({i:rep(pn_int[i-1])})

            if i==30: # update absolute cross-section values
                brem_cs.append(rep(brem_val))
                pair_cs.append(rep(pair_int[-1]))
                pn_bb_cs.append(rep(pn_bb_val))
                pn_cs.append(rep(pn_int[-1]))

        with open('brem_cs_%s_%s.npy' % (lepton,material), 'wb') as f:
            np.save(f, np.asarray(brem_cs))

        with open('pair_cs_%s_%s.npy' % (lepton,material), 'wb') as f:
            np.save(f, np.asarray(pair_cs))

        with open('pn_bb_cs_%s_%s.npy' % (lepton,material), 'wb') as f:
            np.save(f, np.asarray(pn_bb_cs))

        with open('pn_cs_%s_%s.npy' % (lepton,material), 'wb') as f:
            np.save(f, np.asarray(pn_cs))

        for j in range(0,31): # to normalize ixc's
            if brem_dict[E][30]!=0:
                brem_dict[E].update({j:brem_dict[E][j]/brem_dict[E][30]})
            if pair_dict[E][30]!=0:
                pair_dict[E].update({j:pair_dict[E][j]/pair_dict[E][30]})
            if pn_bb_dict[E][30]!=0:
                pn_bb_dict[E].update({j:pn_bb_dict[E][j]/pn_bb_dict[E][30]})
            if pn_dict[E][30]!=0:
                pn_dict[E].update({j:pn_dict[E][j]/pn_dict[E][30]})

    brem_dframe = pd.DataFrame.from_dict(brem_dict,orient='index').transpose()
    pair_dframe = pd.DataFrame.from_dict(pair_dict,orient='index').transpose()
    pn_bb_dframe = pd.DataFrame.from_dict(pn_bb_dict,orient='index').transpose()
    pn_dframe = pd.DataFrame.from_dict(pn_dict,orient='index').transpose()

    ixc_dict = {'brem':brem_dframe,'pair':pair_dframe,'pn_bb':pn_bb_dframe, 'pn':pn_dframe}

    try:
        for model in ixc_dict.keys():
            Data.add_ixc(str(lepton),ixc_dict,model=model, material = material)
    except ValueError or TypeError:
        pass

    f = open("ixc_%s_%s.pkl" % (lepton,material),"wb")
    pickle.dump(ixc_dict,f)
    f.close()

    xc_dict = {'brem':np.asarray(brem_cs), 'pair':np.asarray(pair_cs), 'pn_bb':np.asarray(pn_bb_cs), 'pn':np.asarray(pn_cs)}

    Data.add_xc(str(lepton), xc_dict, material=material)

    return 'Problem in calc_ixc function'

# @njit(nogil=True)
def em_cont_part(E_init, alpha_val, beta_val, x): # calculate continuous energy loss part for the stochastic process
    bx = beta_val * x
    if bx < 1e-6:
        E_fin = E_init * (1-bx) - alpha_val*x
    else:
        E_fin = E_init * np.exp(-bx) - alpha_val/beta_val*(1-np.exp(-bx))

    if E_fin<0:E_fin = m_tau
    return E_fin

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
#   Plot Beta - Continuous
# =============================================================================
def plot_c(model):
    plt.figure(2)

    brem = Data.get_beta(particle=lepton,material=material,pn_model=model,beta_type='continuous')['brem']
    pair = Data.get_beta(particle=lepton,material=material,pn_model=model,beta_type='continuous')['pair']
    pn = Data.get_beta(particle=lepton,material=material,pn_model=model,beta_type='continuous')['pn']

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
    plt.title("Beta - Continuous (beta)")
    plt.legend(loc='best')
    plt.xlabel(r"$E_{\tau}$ [GeV]") # or E_{\tau}
    plt.ylabel(r"$\beta_{std~rock}$ [cm$^2$/g]")
    plt.xlim(100, 1e9)
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth='0.2', color='black')
    # plt.savefig('Continuous.png', format='png', dpi = 300)
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
def plot_cs():
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
# Test
# =============================================================================
if __name__ == "__main__":
    lepton = 'tau'
    material = 'water'

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

    # print(em_cont_part(1e4, 1e-3, 1e-3, 1))
    # calc_alpha()
    betac = calc_beta_continuous()
    calc_beta_total()
    calc_ixc()
