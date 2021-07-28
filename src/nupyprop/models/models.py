#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 24 19:41:29 2021

@author: sam
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 3 14:17:24 2020

@author: sam
"""
import nupyprop.data as Data
# import data as Data

import numpy as np
import scipy.integrate as integrate
import scipy.constants as scc
import multiprocessing as mp
from multiprocessing import Pool
# import matplotlib.pyplot as plt
# mpl.use('Agg') # for clusters
from astropy.table import Table
from astropy.io import ascii
from collections import OrderedDict
import time

E_nu = Data.E_nu
E_lep = Data.E_lep

m_e = scc.physical_constants["electron mass energy equivalent in MeV"][0]*1e-3 # GeV
m_mu = scc.physical_constants["muon mass energy equivalent in MeV"][0]*1e-3 # GeV
m_tau = scc.physical_constants["tau mass energy equivalent in MeV"][0]*1e-3 # GeV
m_pi = 139.57018e-3 # pion mass in GeV
alpha_fs = scc.fine_structure
le = 3.8616e-11 # electron Compton wavelength in cm

m_p = scc.physical_constants["proton mass energy equivalent in MeV"][0]*1e-3 # GeV
G_F = scc.physical_constants["Fermi coupling constant"][0] # GeV^-2
N_A = scc.Avogadro

def rep(val): # Remove non-physical values
    '''

    Parameters
    ----------
    val : float
        Usually negative or np.nan values, to set to 0.

    Returns
    -------
    float
        Returns 0 in case of non-physical (NaN) entries or returns val.

    '''
    if (val<0) or (np.isnan(val)==True):
        return 0
    else:
        return val
    return 'Problem in rep function'

def short_int_nquad(fn, dep_fn_lim, indep_fn_lim_high, arr, args):
    '''

    Parameters
    ----------
    fn : NONE
        Lamba function (integrating function).
    dep_fn_lim : NONE
        Function (dependent function limit).
    indep_fn_lim_high : NONE
        Function (independent function, upper limit).
    arr : ndarray
        1D array containing y-values you want to integrate over.
    args : tuple
        Fixed parameters for integration.

    Returns
    -------
    ndarray
        1D array containing integrated results.

    '''
    ans = []
    result = integrate.nquad(fn, [dep_fn_lim,[arr[0],indep_fn_lim_high]],args=args)[0]
    ans.append(result)
    for i in range(len(arr)-1):
        result+= integrate.nquad(fn, [dep_fn_lim,[arr[i+1],arr[i]]],args=args)[0]
        ans.append(result)
    return np.asarray(ans)

# ============================================================================
# Connolly, Thorne & Waters (CTW) Neutrino & Anti-Neutrino Cross Section Model
# ============================================================================
def ctw_xc(): # CTW parameterization (eq. 7)
    '''

    Returns
    -------
    Creates astropy tables with the CTW neutrino/anti-neutrino-nucleon cross-section values.

    '''
    C_n_cc = np.array([-1.826,-17.31,-6.406,1.431,-17.91]) # neutrino; CC
    C_n_nc = np.array([-1.826,-17.31,-6.448,1.431,-18.61]) # neutrino; NC
    C_an_cc = np.array([-1.033,-15.95,-7.247,1.569,-17.72]) # anti-neutrino; CC
    C_an_nc = np.array([-1.033,-15.95,-7.296,1.569,-18.30]) # anti-neutrino; NC

    sigma = lambda E, c_0,c_1,c_2,c_3,c_4:c_1+c_2*np.log(np.log10(E)-c_0)+c_3*np.log(np.log10(E)-c_0)**2+c_4/(np.log(np.log10(E)-c_0))

    sigma_nu_cc = 10**np.asarray([sigma(i,C_n_cc[0],C_n_cc[1],C_n_cc[2],C_n_cc[3],C_n_cc[4]) for i in E_nu]) # neutrino; CC

    sigma_nu_nc = 10**np.asarray([sigma(i,C_n_nc[0],C_n_nc[1],C_n_nc[2],C_n_nc[3],C_n_nc[4]) for i in E_nu]) # neutrino; NC

    nu_xc_meta = OrderedDict({'Description':'Neutrino-nucleon cross-section values for CTW',
                               'energy':'Neutrino energy, in GeV',
                               'sigma_cc':'Charged current cross-section for CTW, in cm^2',
                               'sigma_nc':'Neutral current cross-section for CTW, in cm^2'})

    nu_xc_table = Table([E_nu, sigma_nu_cc, sigma_nu_nc], names=('energy','sigma_cc_ctw','sigma_nc_ctw'), meta=nu_xc_meta) # neutrino table

    fnm_nu = "xc_neutrino_ctw.ecsv"
    ascii.write(nu_xc_table, fnm_nu, format='ecsv', fast_writer=False, overwrite=True)


    sigma_anu_cc = 10**np.asarray([sigma(i,C_an_cc[0],C_an_cc[1],C_an_cc[2],C_an_cc[3],C_an_cc[4]) for i in E_nu]) # anti-neutrino; CC

    sigma_anu_nc = 10**np.asarray([sigma(i,C_an_nc[0],C_an_nc[1],C_an_nc[2],C_an_nc[3],C_an_nc[4]) for i in E_nu]) # anti-neutrino; NC

    anu_xc_meta = OrderedDict({'Description':'Anti_neutrino-nucleon cross-section values for CTW',
                               'energy':'Anti_neutrino energy, in GeV',
                               'sigma_cc':'Charged current cross-section for CTW, in cm^2',
                               'sigma_nc':'Neutral current cross-section for CTW, in cm^2'})

    anu_xc_table = Table([E_nu, sigma_anu_cc, sigma_anu_nc], names=('energy','sigma_cc_ctw','sigma_nc_ctw'), meta=anu_xc_meta) # anti-neutrino table

    fnm_anu = "xc_anti_neutrino_ctw.ecsv"
    ascii.write(anu_xc_table, fnm_anu, format='ecsv', fast_writer=False, overwrite=True)

    return None

# ===========================================================================
# Connolly, Thorne & Waters (CTW) Neutrino & Anti-Neutrino Cross Section CDFs
# ===========================================================================
def ctw_ixc(): # CTW parameterization (eqs. 12, 13)
    '''

    Returns
    -------
    Creates astropy tables with the CTW neutrino/anti-neutrino-nucleon cross-section CDF values.

    '''
    A_high_n_cc = np.array([-0.008,0.26,3.0,1.7])
    A_high_n_nc = np.array([-0.005,0.23,3.0,1.7]) # -0.005; changed 26/7/21
    A_high_an_cc = np.array([-0.0026,0.085,4.1,1.7])
    A_high_an_nc = np.array([-0.005,0.23,3.0,1.7])

    ixc_high = lambda y,E,A_0,A_1,A_2,A_3,ymin,ymax:(np.log((y - (A_0 - A_1 * np.exp(-(np.log10(E) - A_2)/A_3)))/(ymin - (A_0 - A_1 * np.exp(-(np.log10(E) - A_2)/A_3)))))/(np.log((ymax - (A_0 - A_1 * np.exp(-(np.log10(E) - A_2)/A_3)))/(ymin - (A_0 - A_1 * np.exp(-(np.log10(E) - A_2)/A_3))))) # eq. 13

    v2 = -np.linspace(0.1,3,num=30)
    yvals = 10**v2 # The integrated cross-section values should go from y = 0, 10^(-0.1),...,10^(-3). This is a convention we chose to adopt. (This is also why we don't use the CTW ixc_low)


    ymin = 1e-3 # high y region only
    ymax = 1 # high y region only

    C_n_cc, C_n_nc, C_an_cc, C_an_nc = [],[],[],[] # initialize CDF arrays

    for E in E_nu:

        nu_cc, nu_nc, anu_cc, anu_nc = np.zeros(1),np.zeros(1),np.zeros(1),np.zeros(1) # pad all arrays with 0 at the beginning

        for i in range(1,31):
            nu_cc = np.concatenate((nu_cc, np.asarray([1-ixc_high(yvals[i-1],E,A_high_n_cc[0],A_high_n_cc[1],A_high_n_cc[2],A_high_n_cc[3],ymin,ymax)])))

            nu_nc = np.concatenate((nu_nc, np.asarray([1-ixc_high(yvals[i-1],E,A_high_n_nc[0],A_high_n_nc[1],A_high_n_nc[2],A_high_n_nc[3],ymin,ymax)])))

            anu_cc = np.concatenate((anu_cc, np.asarray([1-ixc_high(yvals[i-1],E,A_high_an_cc[0],A_high_an_cc[1],A_high_an_cc[2],A_high_an_cc[3],ymin,ymax)])))

            anu_nc = np.concatenate((anu_nc, np.asarray([1-ixc_high(yvals[i-1],E,A_high_an_nc[0],A_high_an_nc[1],A_high_an_nc[2],A_high_an_nc[3],ymin,ymax)])))


        C_n_cc.append(nu_cc)
        C_n_nc.append(nu_nc)
        C_an_cc.append(anu_cc)
        C_an_nc.append(anu_nc)

    C_n_cc = np.asarray(C_n_cc).flatten()
    C_n_nc = np.asarray(C_n_nc).flatten()
    C_an_cc = np.asarray(C_an_cc).flatten()
    C_an_nc = np.asarray(C_an_nc).flatten()

    yvals_padded= np.insert(yvals,0,0)

    nu_ixc_meta = OrderedDict({'Description':'Neutrino-nucleon cross-section CDF values for CTW',
                            'energy':'Neutrino energy, in GeV',
                            'y':'Inelasticity; y = (E_initial-E_final)/E_initial',
                            'cc_cdf_ctw':'Charged current cross-section CDF values for CTW',
                            'nc_cdf_ctw':'Neutral current cross-section CDF values for CTW',
                            'Note':'The integrated cross-section CDF values should be integrated from y = 0, 10^(-0.1),...,10^(-3). This is a convention we chose to adopt'})

    nu_ixc_table = Table([np.repeat(E_nu,31), np.tile(yvals_padded,len(E_nu)), C_n_cc, C_n_nc], names=('energy','y','cc_cdf_ctw','nc_cdf_ctw'), meta=nu_ixc_meta)
    fnm_nu = "ixc_neutrino_ctw.ecsv"
    ascii.write(nu_ixc_table, fnm_nu, format='ecsv', fast_writer=False, overwrite=True)

    anu_ixc_meta = OrderedDict({'Description':'Anti_neutrino-nucleon cross-section CDF values for CTW',
                            'energy':'Anti_neutrino energy, in GeV',
                            'y':'Inelasticity; y = (E_initial-E_final)/E_initial',
                            'cc_cdf_ctw':'Charged current cross-section CDF values for CTW',
                            'nc_cdf_ctw':'Neutral current cross-section CDF values for CTW',
                            'Note':'The integrated cross-section CDF values should be integrated from y = 0, 10^(-0.1),...,10^(-3). This is a convention we chose to adopt'})

    anu_ixc_table = Table([np.repeat(E_nu,31), np.tile(yvals_padded,len(E_nu)), C_an_cc, C_an_nc], names=('energy','y','cc_cdf_ctw','nc_cdf_ctw'), meta=anu_ixc_meta)
    fnm_anu = "ixc_anti_neutrino_ctw.ecsv"
    ascii.write(anu_ixc_table, fnm_anu, format='ecsv', fast_writer=False, overwrite=True)

    return None

# =============================================================================
# BDHM F2 Structure Function
# =============================================================================
def f2(q2, y, E, z, A, model): # BDHM parameterization for f2; PhysRevD.89.094027
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
        F2 (EM) structure function for different models.

    '''
    x = q2/(2*m_p*E*y) # x*y*S = Q^2; where S = 2*m_p*E

    if model == 'bdhm': # PhysRevD.89.094027
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

        if x > 0.1: # x<0.1 constraint; see in paper
            f2p = 0
        else:
            A_A = a_0 + a_1*np.log(1 + q2/mu2) + a_2*np.log(1+q2/mu2)**2 # Since A is reserved for atomic mass number
            B = b_0 + b_1*np.log(1+q2/mu2) + b_2*np.log(1+q2/mu2)**2
            C = c_0 + c_1*np.log(1+q2/mu2)
            D = (q2*(q2 + lambda2*M2))/(q2 + M2)**2

            f2p = D*(1 - x)**n * (C + A_A*np.log(1/x * q2/(q2 + mu2)) + B*np.log(1/x * q2/(q2 + mu2))**2) # [eq.8]; f2_proton

    elif model=='ckmt': # DOI 10.1007/s100529900078
        A_A = 0.1502 # Since A is reserved for atomic mass number
        delta_0 = 0.07684
        B = 1.2064
        alpha_R = 0.4150
        a = 0.2631 # GeV^2
        d = 1.1170 # GeV^2
        b = 0.6452 # GeV^2
        c = 3.5489 # GeV^2

        n_q2 = 3/2 * (1 + q2/(q2+c))
        delta_q2 = delta_0 * (1 + (2*q2)/(q2+d))

        f2_sea = A_A*x**(-delta_q2) * (1-x)**(n_q2+4) * (q2/(q2+a))**(1+delta_q2) # [eq. 2]
        f2_val = B*x**(1-alpha_R) * (1-x)**n_q2 * (q2/(q2+b))**alpha_R # [eq. 4]

        f2p = f2_sea + f2_val # f2_proton

    else: # Your custom F2 model goes here
        # f2p = something
        pass # and uncomment this

    # NB: The following are nuclear shadowing corrections, which are material dependent.

    p = 1 - 1.85*x + 2.45*x**2 - 2.35*x**3 + x**4 # correction since neutrons != protons

    if (x>0 and x<0.0014):
        f2 = A**(-0.1)*(z + (A-z)*p)*f2p # changed 17/12/2020
    elif (x>0.0014 and x<0.04):
        f2 = A**(0.069 * np.log10(x) + 0.097)*(z + (A-z)*p)*f2p  # changed 25/3/2021
    else:
        f2 = (z + (A-z)*p)*f2p  # changed 17/12/2020

    return f2

# ========================================
# Differential Cross-Section for PN Models
# ========================================
def pn(lnq2, y, E, m_le, z, A, model): # a la PhysRevD.63.094020
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

    R = 0

# only one power of q2 in denominator, since this is integral over log(q2)
    dsigmadq2dx = (4*np.pi*alpha_fs**2)/q2 * f2(q2,y,E,z,A,model)/x * (1 - y - (m_p*x*y)/(2*E) + (1 - (2*m_le**2)/q2)*y**2*((1 + 4*m_p**2*x**2/q2)/(2*(1+R)))) # [eq. 3.4]

    y_dsigmadq2dy = y * dsigmadq2dx * q2/(2*m_p*E*y**2) * 0.389e-27 # changing differential variable from x to y

    return (N_A/A) * y_dsigmadq2dy # Note: This is N_A/A * y * dsigma/(dQ^2*dy)

def pn_y_cut(E, m_le, z, A, model):
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
    return [((m_p+m_pi)**2-m_p**2)/(2*m_p*E), 1e-3] # [lower, upper]

def pn_y_tot(E, m_le, z, A, model):
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
    return [((m_p + m_pi)**2 - m_p**2)/(2*m_p*E), (1 - m_le/E)] # [lower, upper]

def pn_q2(y, E, m_le, z, A, model):
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
        Photonuclear Q^2 integration limit

    '''
    return [np.log((m_le**2 * y**2)/(1 - y)), np.log(2*m_p*E*y - ((m_p+m_pi)**2-m_p**2))] # [lower, upper]

def cs_pn(lnq2, y, E, m_le, z, A, model):
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
    dsigmadq2dy = pn(lnq2, y, E, m_le, z, A, model)/y # Note: This is N_A/A * dsigma/(dQ^2*dy)
    return dsigmadq2dy

# ============================================================================
# Calculate beta, cross-section and cross-section CDF values for PN BDHM Model
# ============================================================================

def calc_beta(lepton, material, model):

    if lepton=='tau':m_le = m_tau
    elif lepton=='muon':m_le=m_mu

    if material=='water':
        z = 6.6
        A = 11.89
    elif material=='rock':
        z = 11.0
        A = 22.0

    beta_cut = np.asarray([integrate.nquad(pn, [pn_q2, pn_y_cut], args=(i, m_le, z, A, model))[0] for i in E_lep])
    beta_cut = np.asarray([rep(i) for i in beta_cut])

    beta_total = np.asarray([integrate.nquad(pn, [pn_q2, pn_y_tot], args=(i, m_le, z, A, model))[0] for i in E_lep])
    beta_total = np.asarray([rep(i) for i in beta_total])

    beta_meta = OrderedDict({'Description':'Model/parameterization dependent photonuclear energy loss lookup table for %s in %s' % (lepton.capitalize(),material),
                             'energy':'%s energy, in GeV' % lepton.capitalize(),
                             'beta_pn_%s_cut' % model:'Photonuclear %s energy loss model beta values integrated from y_min to y_max=1e-3, in cm^2/g' % str.upper(model),
                             'beta_pn_%s_total' % model:'Photonuclear %s energy loss model beta values integrated from y_min to y_max, in cm^2/g' % str.upper(model)})
    beta_table = Table([E_lep, beta_cut, beta_total], names=('energy','beta_pn_%s_cut' % model,'beta_pn_%s_total' % model), meta=beta_meta)

    fnm = "beta_%s_pn_%s_%s.ecsv" % (lepton,model,material)
    ascii.write(beta_table, fnm, format='ecsv', fast_writer=True, overwrite=True)

    return None

def calc_xc(lepton, material, model):

    yvals = 10**(-np.linspace(0.1,3,num=30)) # The integrated cross-section values should go from y = 0, 10^(-0.1),...,10^(-3). This is a convention we chose to adopt.

    yvals_padded = np.insert(yvals, 0, 0)

    if lepton=='tau':m_le = m_tau
    elif lepton=='muon':m_le=m_mu

    if material=='water':
        z = 6.6
        A = 11.89
    elif material=='rock':
        z = 11.0
        A = 22.0

    ixc_lst = [] # initialize empty array to put cross-section CDF values in
    xc_arr = []  # initialize empty array to put (absolute) cross-section values in


    for E in E_lep:

        args = (E, m_le, z, A, model)

        pn_int = short_int_nquad(cs_pn, pn_q2, (1-m_le/E), yvals, args=args) # pn_int will contain values integrated from y =  10^(-0.1) to y = 10^(-3); (1-m_le/E) is the y_max value to integrate upto

        pn_int = np.asarray([rep(i) for i in pn_int]) # remove non-physical values

        xc_arr.append(pn_int[-1]) # absolute cross-section value is the last element in pn_int since that's integrated from y_min to y_max

        pn_int = np.asarray([i/pn_int[-1] if pn_int[-1]!=0 else i for i in pn_int]) # normalize to CDF values
        # pn_int = pn_int/np.max(pn_int) # normalize to CDF values

        pn_int = np.insert(pn_int, 0, 0) # pad with 0 at the beginning of the array

        ixc_lst.append(pn_int)

    cdf = np.asarray(ixc_lst).flatten()

    xc_meta = OrderedDict({'Description':'%s-nucleon cross-section values for PN_%s in %s' % (lepton.capitalize(),str.upper(model),material),
                           'energy':'%s energy, in GeV' % lepton.capitalize(),
                           'sigma_%s' % model:'N_A/A * cross-section for PN_%s in %s, in cm^2/g' % (str.upper(model),material)})
    xc_table = Table([E_lep, xc_arr], names=('energy','sigma_pn_%s' % model), meta=xc_meta)

    fnm_xc = "xc_%s_pn_%s_%s.ecsv" % (lepton,model,material)
    ascii.write(xc_table, fnm_xc, format='ecsv', fast_writer=True, overwrite=True)

    ixc_meta = OrderedDict({'Description':'%s-nucleon cross-section CDF values for PN_%s in %s' % (lepton.capitalize(),str.upper(model),material),
                            'energy':'%s energy, in GeV' % lepton.capitalize(),
                            'y':'Inelasticity; y = (E_initial-E_final)/E_initial',
                            'cdf':'Cross-section CDF values for PN_%s in %s' % (str.upper(model),material),
                            'Note':'The integrated cross-section CDF values should be integrated from y = 0, 10^(-0.1),...,10^(-3). This is a convention we chose to adopt'})
    ixc_table = Table([np.repeat(E_lep,31), np.tile(yvals_padded,len(E_lep)), cdf], names=('energy','y','cdf_pn_%s' % model), meta=ixc_meta)

    fnm_ixc = "ixc_%s_pn_%s_%s.ecsv" % (lepton,model,material)
    ascii.write(ixc_table, fnm_ixc, format='ecsv', fast_writer=True, overwrite=True)

    return None

# =============================================================================
# Test
# =============================================================================
if __name__ == "__main__":
    pass
    # start_time = time.time()
    # ctw_xc()
    # ctw_ixc()

    # input_list = [['muon','water','bdhm'],['muon','rock','bdhm'],['muon','water','ckmt'],['muon','rock','ckmt'],['tau','water','bdhm'],['tau','rock','bdhm'],['tau','water','ckmt'],['tau','rock','ckmt']]
    # p = Pool(mp.cpu_count()) # use all available cores
    # p.starmap(calc_beta, input_list)
    # p.starmap(calc_xc, input_list)
    # p.close()

    # end_time = time.time()
    # print(f"It took a total of {end_time-start_time:.2f} seconds to compute")