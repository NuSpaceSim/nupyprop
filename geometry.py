#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 8 11:58:16 2020

@author: sam
"""
# import h5py
import data as Data

import numpy as np
from scipy import integrate
import multiprocessing as mp
from multiprocessing import Pool
from tables import *
import pandas as pd
from pandas import HDFStore,DataFrame
from scipy import interpolate
import sympy as sp
import warnings
from numba import njit,prange
from math import isclose


pd.set_option("display.max_rows", None, "display.max_columns", None)
warnings.filterwarnings('ignore')


Re = 6371.0 # radius of the earth in km
Rlay = np.array([1221.5, 3480.0, 5701.0, 5771.0, 5971.0, 6151.0, 6346.6, 6356.0, 6368.0, 6371.0]) # PREM layers based on R_earth
rho_water = 1.02 # density of water in g/cm^3
beta_arr = np.asarray([float('{:.1f}'.format(i)) for i in np.concatenate((np.linspace(0.1,5,50), np.linspace(6,90,85)))])


    # if idepth is None:idepth = Geometry.def_idepth # use the def_idepth of 4
    # else:idepth = Geometry.def_idepth = idepth # so that the change is carried out in other classes

    # Re = Geometry.Re
    # Rlay = Geometry.Rlay
    # rho_water = Geometry.rho_water
    # betad = Geometry.betad

# def sagitta(tnadir):
#     '''

#     Parameters
#     ----------
#     tnadir : float
#         tnadir in rad.

#     Returns
#     -------
#     sagitta : float
#         The sagitta in km.

#     '''
#     sagitta = Re*(1.0 - np.sin(tnadir))
#     return sagitta

def sagitta_deg(beta_deg):
    '''

    Parameters
    ----------
    beta_deg : float
        Beta in degrees

    Returns
    -------
    sagitta : float
        The sagitta in km.

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
        tnadir in rad.

    Returns
    -------
    pathlength : float
        Path length in km.

    '''
    pathlength = 2*Re*np.cos(tnadir)
    return pathlength

def trajlength(beta_deg):
    '''

    Parameters
    ----------
    beta_deg : float
        Beta in degrees.

    Returns
    -------
    trajlength : float
        Trajectory length in km.

    '''
    tnadir = (90.0-beta_deg)*(np.pi/180.0)
    trajlength = Re*np.cos(tnadir)*2
    return trajlength

@njit(nogil=True)
def PREMdensity(Rin):
    '''

    Parameters
    ----------
    Rin : float
        Distance from the Earth Center (at sagitta), in km.

    Returns
    -------
    rhoOUT : float
        Density in g/cm^3.

    '''
    # global Rlay # be very careful here!
    # Rlay2 = Rlay
    Rlay_2 = np.copy(Rlay)
    Rlay_2[8] = 6368.0+(3.0-float(idepth))

    x=Rin
    y=x/Re

    if (x<=Rlay_2[0]):
        # edens=13.0885-8.8381*y**2
        return 13.0885-8.8381*y**2
    elif (x<=Rlay_2[1]):
        # edens=12.5815-1.2638*y-3.6426*y**2-5.5281*y**3
        return 12.5815-1.2638*y-3.6426*y**2-5.5281*y**3
    elif (x<=Rlay_2[2]):
        # edens=7.9565-6.4761*y+5.5283*y**2-3.0807*y**3
        return 7.9565-6.4761*y+5.5283*y**2-3.0807*y**3
    elif (x<=Rlay_2[3]):
        # edens=5.3197-1.4836*y
        return 5.3197-1.4836*y
    elif (x<=Rlay_2[4]):
        # edens=11.2494-8.0298*y
        return 11.2494-8.0298*y
    elif (x<=Rlay_2[5]):
        # edens=7.1089-3.8045*y
        return 7.1089-3.8045*y
    elif (x<=Rlay_2[6]):
        # edens=2.6910+0.6924*y
        return 2.6910+0.6924*y
    elif (x<=Rlay_2[7]):
        # edens=2.900
        return 2.900
    elif (x<=Rlay_2[8]):
        # edens=2.600
        return 2.600
    elif (x<=Rlay_2[9]):
        # edens=1.020
        return 1.020
    elif (x<=Rlay_2[9]*1.001): # too close to call!
        return 1.020
    else:
        # edens=0.
        return 0.
    # rhoOUT=edens
    # return rhoOUT

# def PREMgramVSang(z):
#     '''

#     Parameters
#     ----------
#     z : float or integer
#         Angle in degrees of trajectory (relative to tangent to surface).

#     Returns
#     -------
#     gramlen : float
#         The "grammage", column density in g/cm2.

#     '''
#     Rlay_3 = np.copy(Rlay)
#     Rlay_3[8] = 6368.0+(3.0-float(idepth))
#     y = z*(np.pi/180)
#     Chord = 2.0*Re*np.sin(y)
#     Depth = Re-(0.5*np.sqrt(4.0*Re**2-Chord**2))
#     Rdep = Re-Depth

#     Rlen = np.zeros(10)
#     RlenS = np.zeros(10)
#     clen = np.zeros(10)
#     glen = np.zeros(10)

#     ilay = np.asarray([1 if Rdep<Rlay_3[i] else 0 for i in range(len(Rlay_3))])
#     gramlen = 0
#     j = 0
#     for i in range(10):
#         if Rdep < Rlay_3[i] and ilay[i]==1 and j==0: # run once only; check by j==0
#             j=i+1
#             Rlen[i] = Rlay_3[i] - Rdep
#             clen[i] = 2*np.sqrt(Rlay_3[i]**2 - (Rlay_3[i] - Rlen[i])**2)
#             # print(clen[i])
#             Rin = Rlay_3[i] - (Rlen[i]/2.0)
#             rho = PREMdensity(Rin)
#             glen[i] = clen[i]*1e5*rho
#             # print(glen)
#         # elif Rdep < Rlay_3[i]:
#         elif ilay[i]>0: # changed 21/12/2020
#             Rlen[i] = Rlay_3[i] - Rlay_3[i-1]
#             print('i = ', i)
#             # print(Rlen) # works until here!
#             print(Rlen[i])
#             RlenS[i] = Rlen[i]
#             print(RlenS)

#             for k in range(j,i):
#                 RlenS[i] = RlenS[i] + Rlen[k]
#                 # print(RlenS[i])
#             clen[i] = 2.0*np.sqrt(Rlay_3[i]**2 - (Rlay_3[i] - RlenS[i])**2)
#             # print(clen[i])
#             for k in range(j,i):
#                 clen[i] = clen[i]-clen[k]
#             Rin = Rlay_3[i] - (Rlen[i]/2.0)
#             rho = PREMdensity(Rin)
#             glen[i] = clen[i]*1e5*rho
#         gramlen = gramlen + glen[i]
#     return gramlen


def PREMgramVSang(z):
    '''

    Parameters
    ----------
    z : float or integer
        Angle in degrees of trajectory (relative to tangent to surface).

    Returns
    -------
    gramlen : float
        The "grammage", column density in g/cm2.

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
            rho = PREMdensity(Rin)
            glen[i] = clen[i]*1e5*rho

        elif ilay[i] > 0:
            Rlen[i] = Rlay_3[i] - Rlay_3[i-1]
            RlenS[i] = Rlay_3[i] - Rlay_3[i-1]

            for j in range(ifirst,i):
                RlenS[i] = RlenS[i] + Rlen[j]
                # print(RlenS[i])

            clen[i] = 2.0*np.sqrt(Rlay_3[i]**2 - (Rlay_3[i] - RlenS[i])**2)
            # print(clen[i])

            for j in range(ifirst,i):
                clen[i] = clen[i] - clen[j]
                # print(clen[i])

            Rin = Rlay_3[i] - (Rlen[i]/2.0)
            rhoOUT = PREMdensity(Rin)
            # print(rhoOUT)
            glen[i] = clen[i]*1e5*rhoOUT
            # print(glen[i])

        gramlen += glen[i]
        # print(gramlen)

    return gramlen

def columndepth(beta_deg):
# def columndepth(tnadir):
    '''

    Parameters
    ----------
    tnadir : float
        tnadir in rad.

    Returns
    -------
    columndepth : float
        Column density in g/cm^2.

    '''
    # tnadird = tnadir * (180.0/np.pi)
    # z = 90.0 - tnadird
    z = beta_deg
    # print('z = ', z)
    if z<0.5:
        z1 = 1.0
        c = PREMgramVSang(z1)
        # columndepth = c*(z/z1)
        columndepth = c*(z/z1-1.0/6.0*(z/z1)**3)
    else:
        columndepth = PREMgramVSang(z)
    return columndepth

@njit(nogil=True)
def densityatx(x, beta_deg):
    '''

    Parameters
    ----------
    x : float
        Position to find the density at, in km.
    beta_deg : float
        beta in degrees.

    Returns
    -------
    rho_at_x : float
        Density at the position x, in g/cm^3.

    '''
    # ell = trajlength(beta_deg)
    tnadir = (90.0-beta_deg)*(np.pi/180.0)
    ell = Re*np.cos(tnadir)*2
    r2 = x**2 - (ell*x) + Re**2
    if beta_deg<5.0:
        r = Re*(1.0 + 0.5*(x**2-ell*x)/Re**2)
    else:
        r = np.sqrt(r2)
    rho_at_x = PREMdensity(r)

    return r, rho_at_x

def cdtot(x_v, beta_deg):
    '''

    Parameters
    ----------
    x_v : lambda function
        Integration variable.
    beta_deg : float
        beta in degrees.

    Returns
    -------
    cdtot : float
        Integrated column density

    '''
    r, rho = densityatx(x_v, beta_deg)
    cdtot = rho*1e5
    return cdtot

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

    p = Pool(mp.cpu_count()) # use all available cores

    betad_fix = np.asarray([np.ones(100)*i for i in beta_arr]).flatten() # used for storing in hdf5 and matching len of cdalong array which uses list comprehension for faster calcs.
    chord = np.asarray([trajlength(float(i)) for i in beta_arr])
    dx = np.asarray([i/100.0 for i in chord])
    xalong = np.asarray([i*float(j) for i in dx for j in range(1,101)])
    cdalong = np.asarray([p.map(integrator,[[cdtot, i, 0, j]])[0][0] for i,j in zip(betad_fix,xalong)])
    p.close()
    return betad_fix, xalong, cdalong

def gen_water_trajs():
    dw = idepth
    Rrock = Re-dw
    chord = np.asarray([trajlength(i) for i in beta_arr])
    water = np.asarray([0.5*(i - np.sqrt(i**2 - 4.0 * (dw**2 + 2.0*Rrock*dw))) if sagitta_deg(j)>dw else i for i,j in zip(chord, beta_arr)])
    return beta_arr, chord, water

def find_interface():
    # global idepth
    y = sp.Symbol('y')
    eqn = sp.Eq(Re*(1-sp.sin((90-y)*(np.pi/180))), idepth)
    return sp.solve(eqn)

def create_traj_table():
    try:
        pd.read_hdf('lookup_tables.h5','Earth/traj_%s/Column_Trajectories' % str(int(idepth)))[0:2]
        return print("idepth = %s already exists in the lookup table. Will initialize that data." % (int(idepth)))
    except (KeyError, FileNotFoundError) as e:
        beta_col, xalong, cdalong = gen_col_trajs()
        dataset_col = pd.DataFrame({'beta':beta_col, 'xalong':xalong, 'cdalong':cdalong}).sort_values(by=['beta','xalong','cdalong'])

        beta_water, chord, water = gen_water_trajs()
        dataset_water = pd.DataFrame({'beta':beta_water, 'chord':chord, 'water':water})

        # hdf = HDFStore('lookup_tables.h5','a') # create an add_trajs function in data file
        # hdf.put('Earth/traj_%s/Column_Trajectories' % str(int(idepth)),dataset_col)
        # hdf.put('Earth/traj_%s/Water_Trajectories' % str(int(idepth)),dataset_water)
        # hdf.close()

        Data.add_trajs('col', int(idepth), dataset_col)
        Data.add_trajs('water', int(idepth), dataset_water) # yikes! fixed!
        return None


# def interpcd2distd(angle, depth):
#     dataset = pd.read_hdf('lookup_tables.h5','Earth/traj_%s/Column_Trajectories' % str(int(idepth))) # create an get_trajs function in data file
#     dataset = Data.get_trajs('col', int(idepth))

#     dataset_sliced = dataset[dataset['beta']==angle] # filter according to i/p angle
#     xalong = np.asarray(dataset_sliced['xalong'])
#     cdalong= np.asarray(dataset_sliced['cdalong'])

#     f3 = interpolate.interp1d(cdalong, xalong, kind = 'quadratic')
#     if depth < min(cdalong):
#         x = (depth/cdalong[0])*xalong[0]
#     elif depth > max(cdalong):
#         x = max(xalong)
#     else:
#         x = f3(depth)
#     # return xalong,cdalong,x
#     return x

# def get_traj(**kwargs):
#     dataset = pd.read_hdf('lookup_tables.h5','Earth/traj_%s/Column_Trajectories' % str(int(idepth)))

#     if "angle" not in kwargs:
#         dataset_sliced = dataset
#         try:
#             if kwargs['output'] == 'dat':
#                 np.savetxt('column_%s.dat' % str(int(idepth)), dataset.values, fmt='%.4e', delimiter="\t", header="beta\txalong\tcdalong")
#                 return print("column_%s.dat created" % str(int(idepth)))
#             else:
#                 dataset.to_html('column_%s.html' % str(int(idepth)))
#                 return print("column_%s.html created" % str(int(idepth)))

#         except KeyError:
#             dataset.to_html('column_%s.html' % str(int(idepth)))
#             return print("column_%s.html created" % str(int(idepth)))
#     else:
#         dataset_sliced = dataset[dataset['beta']==float(kwargs['angle'])]
#         try:
#             if kwargs['output'] == 'dat':
#                 np.savetxt('column_%s_%sdeg.dat' % (str(int(idepth)), str(kwargs['angle'])), dataset_sliced.values, fmt='%.4e', delimiter="\t", header="beta\txalong\tcdalong")
#                 return print("column_%s_%sdeg.dat created" % (str(int(idepth)), str(kwargs['angle'])))
#             else:
#                 dataset_sliced.to_html('column_%s_%sdeg.html' % (str(int(idepth)), str(kwargs['angle'])))
#                 return print("column_%s_%sdeg.html created" % (str(int(idepth)), str(kwargs['angle'])))
#         except KeyError:
#             dataset_sliced.to_html('column_%s_%sdeg.html' % (str(int(idepth)), str(kwargs['angle'])))
#             return print("column_%s_%sdeg.html created" % (str(int(idepth)), str(kwargs['angle'])))

#     return print("Error! Please enter an acceptable value of beta or leave blank for all angles")

# def get_water(**kwargs):
#     dataset = pd.read_hdf('lookup_tables.h5','Earth/traj_%s/Water_Trajectories' % str(int(idepth)))

#     if "angle" not in kwargs:
#         dataset_sliced = dataset
#         chord = np.asarray(dataset['chord'])
#         water= np.asarray(dataset['water'])
#         try:
#             if kwargs['output'] == 'dat':
#                 np.savetxt('water_%s.dat' % str(int(idepth)), dataset.values, fmt='%.4e', delimiter="\t", header="beta\tchord\twater")
#                 return print("water_%s.dat created" % str(int(idepth)))
#             else:
#                 dataset.to_html('column_%s.html' % str(int(idepth)))
#                 return print("water_%s.html created" % str(int(idepth)))

#         except KeyError:
#             dataset.to_html('water_%s.html' % str(int(idepth)))
#             return print("water_%s.html created" % str(int(idepth)))
#     else:
#         dataset_sliced = dataset[dataset['beta']==float(kwargs['angle'])]
#         chord = float(dataset_sliced['chord'])
#         water = np.float(dataset_sliced['water'])
#         s = '''\
#             beta = {beta} deg,
#             chord = {chord},
#             water = {water}\
#             '''.format(beta=str(kwargs['angle']), chord=str(chord), water=str(water))
#         return print(s)
#         # return print("beta = %.4e deg,\n chord = %.4e,\n water = %.4e", float(kwargs['angle']), chord, water)

#     return print("Error! Please enter acceptable an value of beta or leave blank for all angles")

# =============================================================================
# Test Zone!!
# =============================================================================
# print(find_interface())
# idepths = np.asarray([0,1,2,3,4,5,6,7,8,9,10])
# for i in idepths:
#     a = Geometry(idepth=i)
#     a.create_traj_table()

# a = Data()
# b = Geometry(idepth = 4)
# b.create_traj_table()
# b.get_traj(angle='0',output='dat')
# b.get_water(output='dat')
# b.my_fun(arg1 = 'first')
# print("For idepth = %d, " % a.idepth +  "the water/rock transition occurs at",'{:.2f}'.format(b.find_interface()[0]) + str(' deg'))
# x=b.interpcd2distd(1,2.2684E+07)

# =============================================================================
#
# =============================================================================
if __name__ == "__main__":
    # import data as Data

    idepth = 4
    print(find_interface()[0])

    # PREMdensity(1203)
    # r, rho = densityatx(6752.231264859498,32.0)
    # print(r,rho)
    # r,rho = densityatx(719814756.1985742, 10)
    # print(r,rho)
    # print(columndepth(10))
    # thn = 1.3962622222222221
    # print(columndepth(thn))
    # print(PREMgramVSang(10))