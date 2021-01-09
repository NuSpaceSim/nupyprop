#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 00:30:47 2020

@author: sam
"""
import data as Data
import energy_loss as Energy_loss
import my_interpolation as Interpolation
import numpy as np
import time
import random
from numba import njit, prange

rho_rock = 2.65

alpha_rock = Data.get_alpha('tau', material='rock')

beta_rock = Data.get_beta('tau', 'rock', 'total', 'allm')

@njit(nogil=True)
def my_rand(): # because we need a random no. in the range of (0,1) and not [0,1)
    random_no = random.random()
    while random_no == 0:
        random_no = random.random()
    return random_no


@njit(nogil=True)
def dedx(energy):
    alpha_val = Interpolation.int_alpha(energy, alpha_rock)
    beta_val = Interpolation.int_beta(energy, beta_rock)
    # print("alpha = ", alpha_val)
    # print("beta = ", beta_val)
    dedx_val = alpha_val + beta_val*energy

    return dedx_val

@njit(nogil=True)
def idecay(energy, x):
    ctau = 87.11e-4
    gamma = energy/1.777
    pdec = 1 - np.exp(-x/(gamma*ctau))
    dy = my_rand()
    if dy < pdec:
        idecay_val = "decayed"
    else:
        idecay_val = "not_decayed"
    return idecay_val

@njit(nogil=True, parallel=True)
# def em_continuous_loss(E_init, alpha_val, beta_val, x): # calculate continuous energy loss part for the stochastic process
def em_continuous_loss(): # calculate continuous energy loss part for the stochastic process
    energy0 = 1e8
    gct = energy0/1.777*87.11e-4
    distance = 20*1e5 # km -> cm

    beta0 = 1e-6 # cm^2/g
    deltax = 1e-2/(beta0*rho_rock) # 1/(1/cm) = cm
    jmax = distance/deltax
    # print("jmax = ", jmax, deltax)
    deltaE = dedx(energy0*1.15)*deltax*rho_rock
    # print("deltaE = ", deltaE)

    nmax = 100000

    one = 1
    zero = 0
    print(energy0, zero, one, zero)

    for k in range(1,31):
        distance = float(k)*1e5*0.5
        deltax = 1e-2/(beta0 * rho_rock)
        jmax = distance/deltax
        # print("distance and jmax = ", distance, jmax)
        psurv = 0

        for n in prange(1, nmax+1):
            d = 0
            p_id = 'not_decayed'
            energy0 = 1e8
            energy = energy0

            for j in range(1, int(jmax)+1):
                d += deltax
                if p_id == 'not_decayed':
                    p_id = idecay(energy, deltax)
                else:
                    break # need to break out of j loop
                energy = energy - dedx(energy)*deltax*rho_rock

            if p_id == 'not_decayed': psurv +=1
        err1 = np.sqrt(psurv)/float(nmax)
        psurv = psurv/float(nmax)
        test = np.exp(-distance/gct)
        # print("distance in km = ", distance*1e-5)
        # print("psurv = ", psurv, test, err1)
        # print("%.3e \t %.3e \t %.3e \t %.3e" % (energy0, distance, psurv, err1)) # write in fort.16
        print(energy0, distance, psurv, err1) # write in fort.16

    return None



# =============================================================================
#
# =============================================================================
if __name__ == "__main__":
    start_time = time.time()
    print(em_continuous_loss())

    end_time = time.time()
    print(f"It took {end_time-start_time:.2f} seconds to compute")