#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 18:37:16 2020

@author: sam
"""
import data as Data
import numpy as np
import pandas as pd
import scipy.constants as scc
pd.set_option("display.max_rows", None, "display.max_columns", None)

E_nu = Data.E_nu
m_p = scc.physical_constants["proton mass energy equivalent in MeV"][0]*1e-3 # GeV
G_F = scc.physical_constants["Fermi coupling constant"][0] # GeV^-2
N_A = scc.Avogadro

def hls_xc(model): # Hallsie data files
    sigma_n_cc = np.asarray(np.genfromtxt('./nutables/%s/nucc-xc.dat' % str(model),usecols=1)*1e-38) # neutrino; CC

    sigma_n_nc = np.asarray(np.genfromtxt('./nutables/%s/nunc-xc.dat' % str(model),usecols=1)*1e-38) # neutrino; NC

    sigma_an_cc = np.asarray(np.genfromtxt('./nutables/%s/anucc-xc.dat' % str(model),usecols=1)*1e-38) # anti-neutrino; CC

    sigma_an_nc = np.asarray(np.genfromtxt('./nutables/%s/anunc-xc.dat' % str(model),usecols=1)*1e-38) # anti-neutrino; NC

    xc_dict = {'nu':{'cc':sigma_n_cc,'nc':sigma_n_nc}, 'anu':{'cc':sigma_an_cc,'nc':sigma_an_nc}}

    Data.add_xc('nu', xc_dict, model=model)
    return None


def hls_ixc(model): # Hallsie data files
    particle_current = ['anucc','anunc','nucc','nunc']
    ixc_dict={}
    for particle in particle_current:
        file = './nutables/%s/%s-ymin.dat' % (str(model),str(particle))
        data_dict = {}
        skp = 0
        for energy in E_nu:
            data_lst = sorted(np.concatenate((np.genfromtxt(file,usecols=(0,1,2,3,4,5),skip_header=skp,max_rows=5)*1e-38)),reverse=False) # NOTE: sigma is the last (30th) value (1e-3 -> 1) # units = cm^2
            data_lst = [data_lst[i]/data_lst[-1] for i in range(30)] # normalize all values
            data_lst.insert(0,0) # pad with 0 at the beginning because it is the CDF after all
            data_dict.update({energy:data_lst})
            skp+=6
        dframe = pd.DataFrame.from_dict(data_dict,orient='index',columns=[i for i in range(0,31)])
        ixc_dict.update({particle:dframe.transpose()})

    Data.add_ixc('nu', ixc_dict, model=model)
    return None


def ctw_xc(): # CTW parameterization (eq. 7)

    C_n_cc = np.array([-1.826,-17.31,-6.406,1.431,-17.91]) # neutrino; CC
    C_n_nc = np.array([-1.826,-17.31,-6.448,1.431,-18.61]) # neutrino; NC
    C_an_cc = np.array([-1.033,-15.95,-7.247,1.569,-17.72]) # anti-neutrino; CC
    C_an_nc = np.array([-1.033,-15.95,-7.296,1.569,-18.30]) # anti-neutrino; NC

    sigma = lambda E, c_0,c_1,c_2,c_3,c_4:c_1+c_2*np.log(np.log10(E)-c_0)+c_3*np.log(np.log10(E)-c_0)**2+c_4/(np.log(np.log10(E)-c_0))

    sigma_n_cc = 10**np.asarray([sigma(i,C_n_cc[0],C_n_cc[1],C_n_cc[2],C_n_cc[3],C_n_cc[4]) for i in E_nu]) # neutrino; CC

    sigma_n_nc = 10**np.asarray([sigma(i,C_n_nc[0],C_n_nc[1],C_n_nc[2],C_n_nc[3],C_n_nc[4]) for i in E_nu]) # neutrino; NC

    sigma_an_cc = 10**np.asarray([sigma(i,C_an_cc[0],C_an_cc[1],C_an_cc[2],C_an_cc[3],C_an_cc[4]) for i in E_nu]) # anti-neutrino; CC

    sigma_an_nc = 10**np.asarray([sigma(i,C_an_nc[0],C_an_nc[1],C_an_nc[2],C_an_nc[3],C_an_nc[4]) for i in E_nu]) # anti-neutrino; NC

    xc_dict = {'nu':{'cc':sigma_n_cc,'nc':sigma_n_nc}, 'anu':{'cc':sigma_an_cc,'nc':sigma_an_nc}}

    Data.add_xc('nu', xc_dict, model='ctw')

    return None

def ctw_ixc(): # CTW parameterization (eqs. 12, 13)

    A_high_n_cc = np.array([-0.008,0.26,3.0,1.7])
    A_high_n_nc = np.array([0.005,0.23,3.0,1.7])
    A_high_an_cc = np.array([-0.0026,0.085,4.1,1.7])
    A_high_an_nc = np.array([-0.005,0.23,3.0,1.7])


    ixc_high = lambda y,E,A_0,A_1,A_2,A_3,ymin,ymax:(np.log((y - (A_0 - A_1 * np.exp(-(np.log10(E) - A_2)/A_3)))/(ymin - (A_0 - A_1 * np.exp(-(np.log10(E) - A_2)/A_3)))))/(np.log((ymax - (A_0 - A_1 * np.exp(-(np.log10(E) - A_2)/A_3)))/(ymin - (A_0 - A_1 * np.exp(-(np.log10(E) - A_2)/A_3))))) # CTW eq. 13

    v2 = -np.linspace(0.1,3,num=30)
    yvals = 10**v2 # The integrated cross-section values should go from y = 10^(-0.1) to y = 10^(-3). This is a convention we chose to adopt. (This is also why we don't use the CTW ixc_low)

    nucc_dict = {}
    nunc_dict = {}
    anucc_dict = {}
    anunc_dict = {}

    ymin = 1e-3 # high y region only
    ymax = 1 # high y region only

    for E in E_nu:
        nucc_dict.update({E:{0:0}}) # pad with 0 at the beginning because it is the CDF after all
        nunc_dict.update({E:{0:0}}) # pad with 0 at the beginning because it is the CDF after all
        anucc_dict.update({E:{0:0}}) # pad with 0 at the beginning because it is the CDF after all
        anunc_dict.update({E:{0:0}}) # pad with 0 at the beginning because it is the CDF after all

        for i in range(1,31):
            nucc_dict[E].update({i:1-ixc_high(yvals[i-1],E,A_high_n_cc[0],A_high_n_cc[1],A_high_n_cc[2],A_high_n_cc[3],ymin,ymax)})

            nunc_dict[E].update({i:1-ixc_high(yvals[i-1],E,A_high_n_nc[0],A_high_n_nc[1],A_high_n_nc[2],A_high_n_nc[3],ymin,ymax)})

            anucc_dict[E].update({i:1-ixc_high(yvals[i-1],E,A_high_an_cc[0],A_high_an_cc[1],A_high_an_cc[2],A_high_an_cc[3],ymin,ymax)})

            anunc_dict[E].update({i:1-ixc_high(yvals[i-1],E,A_high_an_nc[0],A_high_an_nc[1],A_high_an_nc[2],A_high_an_nc[3],ymin,ymax)})

    nucc_dframe = pd.DataFrame.from_dict(nucc_dict,orient='index').transpose()
    nunc_dframe = pd.DataFrame.from_dict(nunc_dict,orient='index').transpose()
    anucc_dframe = pd.DataFrame.from_dict(anucc_dict,orient='index').transpose()
    anunc_dframe = pd.DataFrame.from_dict(anunc_dict,orient='index').transpose()

    ixc_dict = {'nucc':nucc_dframe,'nunc':nunc_dframe,'anucc':anucc_dframe,'anunc':anunc_dframe}

    Data.add_ixc('nu', ixc_dict, model='ctw')

    return None

# def Bulmahn_xc(): # Bulmahn parameterization
#     E_c = 3.5e4

    # dsigdy_low = [partial(lambda E,y: (0.389e-27* 2 * m_p * G_F**2/np.pi)*E*((0.19-0.0265*(2.214-np.log10(E_c/E))**2) + (0.036-0.0344*(1.994-np.log10(E_c/E))**2)*(1-y)**2)*(1/(y**(2.3e-2))),E) for E in E_nu[E_nu<E_c]]

#     dsigdy_high = [partial(lambda E,y: (0.389e-27* 2 * m_p * G_F**2/np.pi)*E*((0.060 * (E_c/E)**0.675) + (0.169*(E_c/E)**0.73)*(1-y)**2)*(1/(y**(0.66 * 10**(-1.453 * (np.log10(E_c)/np.log10(E))**6.24)))),E) for E in E_nu[E_nu>E_c]]

#     sig_low = [integrate.quad(i,1e-3,1)[0] for i in dsigdy_low]
#     sig_high = [integrate.quad(i,1e-3,1)[0] for i in dsigdy_high]
#     return None

# def nct15_lhapdf_param_xc(): # Our parameterization
#     C_n_cc = np.array([-2.01883996e+00, -2.03402144e+00, -1.43788314e+01,
#      2.81339770e+00, -2.79130565e+01]) # neutrino; CC
#     C_n_nc = np.array([-4.42333235e+00,  8.69479894e+01, -4.71666188e+01,
#      6.84068713e+00, -1.11437024e+02]) # neutrino; NC
#     C_an_cc = np.array([2.15833824e+00, -3.52953798e+01,  1.06107568e+00,
#      3.55233536e-01, -5.49129575e-04]) # anti-neutrino; CC
#     C_an_nc = np.array([2.37633550e+00, -3.55257626e+01,  1.02102439e+00,
#      3.64798307e-01, -4.40672475e-04]) # anti-neutrino; NC

#     sigma = lambda E, c_0,c_1,c_2,c_3,c_4:c_1+c_2*np.log(np.log10(E)-c_0)+c_3*np.log(np.log10(E)-c_0)**2+c_4/(np.log(np.log10(E)-c_0)) # Our parameterization à la CTW

#     sigma_n_cc = 10**np.asarray([sigma(i,C_n_cc[0],C_n_cc[1],C_n_cc[2],C_n_cc[3],C_n_cc[4]) for i in E_nu]) # neutrino; CC

#     sigma_n_nc = 10**np.asarray([sigma(i,C_n_nc[0],C_n_nc[1],C_n_nc[2],C_n_nc[3],C_n_nc[4]) for i in E_nu]) # neutrino; NC

#     sigma_an_cc = 10**np.asarray([sigma(i,C_an_cc[0],C_an_cc[1],C_an_cc[2],C_an_cc[3],C_an_cc[4]) for i in E_nu]) # anti-neutrino; CC

#     sigma_an_nc = 10**np.asarray([sigma(i,C_an_nc[0],C_an_nc[1],C_an_nc[2],C_an_nc[3],C_an_nc[4]) for i in E_nu]) # anti-neutrino; NC

#     cross_dict = {'nucc':sigma_n_cc, 'nunc':sigma_n_nc, 'anucc':sigma_an_cc, 'anunc':sigma_an_nc}

#     Data.add_xc('nct15_lhapdf_param',cross_dict)

#     return None

def ncteq15_lhapdf_xc(): # Directly from data files
    sigma_n_cc = np.genfromtxt('./nutables/nutables-nCTEQ15-LHAPDF/nct15-nu-cc-xc.dat',usecols=1) # neutrino; CC
    sigma_n_nc = np.genfromtxt('./nutables/nutables-nCTEQ15-LHAPDF/nct15-nu-nc-xc.dat',usecols=1) # neutrino; NC
    sigma_an_cc = np.genfromtxt('./nutables/nutables-nCTEQ15-LHAPDF/nct15-anu-cc-xc.dat',usecols=1) # anti-neutrino; CC
    sigma_an_nc = np.genfromtxt('./nutables/nutables-nCTEQ15-LHAPDF/nct15-anu-nc-xc.dat',usecols=1) # anti-neutrino; NC

    xc_dict = {'nu':{'cc':sigma_n_cc,'nc':sigma_n_nc}, 'anu':{'cc':sigma_an_cc,'nc':sigma_an_nc}}

    Data.add_xc('nu', xc_dict, model='ncteq15_lhapdf')

    return None

def ncteq15_lhapdf_ixc(): # Directly from data files
    particle_current = ['anucc','anunc','nucc','nunc']
    ixc_dict={}
    for particle in particle_current:
        file = './nutables/nutables-nCTEQ15-LHAPDF/%s-ymin.dat' % str(particle)
        data_dict = {}
        skp = 0
        for energy in E_nu:
            data_lst = sorted(np.concatenate((np.genfromtxt(file,usecols=(0,1,2,3,4,5),skip_header=skp,max_rows=5)*1e-38)),reverse=False) # NOTE: sigma is the last (30th) value (1e-3 -> 1) # units = cm^2
            data_lst = [data_lst[i]/data_lst[-1] for i in range(30)] # normalize all values
            data_lst.insert(0,0) # pad with 0 at the beginning because it is the CDF after all
            data_dict.update({energy:data_lst})
            skp+=6
        dframe = pd.DataFrame.from_dict(data_dict,orient='index',columns=[i for i in range(0,31)])
        ixc_dict.update({particle:dframe.transpose()})

    Data.add_ixc('nu', ixc_dict, model='ncteq15_lhapdf')

    return None


# =============================================================================
# Test
# =============================================================================
if __name__ == "__main__":
    # idepth = 4
    pass