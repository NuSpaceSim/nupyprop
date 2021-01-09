#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 21:02:17 2020

@author: sam
"""

import argparse
import data as Data
import geometry as Geometry
import energy_loss as Energy_loss
import my_interpolation as Interpolation
import transport as Transport
import main as Main
import random
import time
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument('-e', '--energy', dest='energy_val', nargs='?', type=float, const=np.array([1e7, 1e8, 1e9, 1e10, 1e11]), default=np.array([1e7, 1e8, 1e9, 1e10, 1e11]), help='value of incoming neutrino energy; defaults are 10^7-10^11 GeV')

# parser.add_argument('-a', '--angle', dest='angle_val', nargs='?', const=np.arange(1,41), default=np.arange(1,41), help='value of angle; defaults are 1-40 degrees')
parser.add_argument('-a', '--angle', dest='angle_val', nargs='?', const=np.array([1,3,5,7,10,12,15,17,20,25,30,35]), default=np.array([1,3,5,7,10,12,15,17,20,25,30,35]), help='value of angle; defaults are 1,3,5,7,10,12,15,17,20,25,30,35 degrees')

parser.add_argument('-i', '--idepth', dest='idepth_val', nargs='?', type=float, const=4, default=4, help='value of idepth; default is 4')

parser.add_argument('-l', '--lepton', dest='lepton_id', nargs='?', type=str, const='tau', default='tau', help='particle for energy loss and propagation - can be tau or muon; default is tau')

parser.add_argument('-m', '--material', dest='material_id', nargs='?', type=str, const='rock', default='rock', help='material for energy loss - can be rock or water; default is rock')

parser.add_argument('-x', '--xc_model', dest='xc_model_id', nargs='?', type=str, const='ncteq15', default='ncteq15', help='neutrino cross-section model; default is ncteq15')

parser.add_argument('-p', '--pn_model', dest='pn_model_id', nargs='?', type=str, const='allm', default='allm', help='lepton photonuclear energy loss model; default is allm')

parser.add_argument('-f', '--fac_nu', dest='fac_nu_val', nargs='?', type=float, const=1.0, default=1.0, help='rescaling for SM cross-sections; default is 1.0')

parser.add_argument('-s', '--stat', dest='stat_val', nargs='?', type=float, const=1e7, default=1e7, help='statistics; default is 1e7')


args = parser.parse_args()
energies = np.asarray([args.energy_val])
angles = args.angle_val
if type(angles) is str:
    angles = angles.split(",") if angles else []
    angles = np.asarray([int(i) for i in angles])

idepth = int(args.idepth_val)
lepton = str(args.lepton_id)
material = str(args.material_id)
cross_section_model = str(args.xc_model_id)
pn_model = str(args.pn_model_id)
fac_nu = float(args.fac_nu_val)
stat = int(args.stat_val)

Data.idepth = Geometry.idepth = Transport.idepth = Main.idepth = idepth
Energy_loss.lepton = Transport.lepton = Interpolation.lepton = Main.lepton = lepton
Energy_loss.material = Transport.material = Interpolation.material = material

if lepton == 'tau':Energy_loss.m_le = Transport.m_le = Energy_loss.m_tau # mass of lepton
else:Energy_loss.m_le = Transport.m_le = Energy_loss.m_mu  # mass of muon

if material=='iso_water':Energy_loss.z=7.0
elif material=='water':Energy_loss.z=6.6 # for NuTauSim Comparison
elif material=='rock':Energy_loss.z=11.0
elif material=='iron':Energy_loss.z=26.0
else: Energy_loss.z=float(input("Enter the atomic number of %s: " % material))

if material=='iso_water':
    Energy_loss.A=14.0
    # Transport.rho=1.02 # g/cm^3
elif material=='water':
    Energy_loss.A=11.89 # for NuTauSim Comparison
    # Transport.rho=1.02 # g/cm^3
elif material=='rock':
    Energy_loss.A=22.0
    # Transport.rho=2.65 # g/cm^3
elif material=='iron':
    Energy_loss.A=55.84
    # Transport.rho=7.6 # g/cm^3
else:
    Energy_loss.A=float(input("Enter the atomic mass of %s: " % material))

Transport.cross_section_model = Main.cross_section_model = cross_section_model
Transport.pn_model = Main.pn_model = pn_model
Transport.fac_nu = Main.fac_nu = fac_nu

Main.stat = stat

if lepton == 'tau':Energy_loss.c_tau = Transport.c_tau = Main.c_tau = 8.703e-3 # c*lifetime, in cm, for taus (taken from PDB 2020)
else:Energy_loss.c_tau = Transport.c_tau = Main.c_tau = 6.586384e4 # c*lifetime, in cm, for muons (taken from PDB 2020)

Main.E_prop = energies
Main.angles = angles

if __name__ == "__main__":
    start_time = time.time()
    prob_dict = Main.main()
    end_time = time.time()
    print(f"It took {end_time-start_time:.2f} seconds to compute")
