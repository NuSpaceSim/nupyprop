#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 21:02:17 2020

@author: sam
"""

import argparse
# import data as Data
# import geometry as Geometry
# import energy_loss as Energy_loss
# import my_interpolation as Interpolation
# import transport as Transport
import main as Main
import random
import time
import numpy as np
import os
import glob

parser = argparse.ArgumentParser()

parser.add_argument('-e', '--energy', dest='energy_val', nargs='?', const=np.array([1e7, 1e8, 1e9, 1e10, 1e11]), default=np.array([1e7, 1e8, 1e9, 1e10, 1e11]), help='value of incoming neutrino energy; defaults are 10^7-10^11 GeV')

parser.add_argument('-a', '--angle', dest='angle_val', nargs='?', const=np.arange(1,41), default=np.arange(1,36), help='value of angle; defaults are 1-35 degrees')

# parser.add_argument('-a', '--angle', dest='angle_val', nargs='?', const=np.arange(1,36), default=np.arange(1,36), help='value of angle; defaults are 1,3,5,7,10,12,15,17,20,25,30,35 degrees')

# parser.add_argument('-a', '--angle', dest='angle_val', nargs='?', const=np.array([1,2,3,5,7,10,12,15,17,20,25,30,35]), default=np.array([1,2,3,5,7,10,12,15,17,20,25,30,35]), help='value of angle; defaults are 1,3,5,7,10,12,15,17,20,25,30,35 degrees')

parser.add_argument('-i', '--idepth', dest='idepth_val', nargs='?', type=float, const=4, default=4, help='value of water layer depth; default is 4')

parser.add_argument('-l', '--lepton', dest='lepton_id', nargs='?', type=str, const='tau', default='tau', help='particle for energy loss and propagation - can be tau or muon; default is tau')

parser.add_argument('-t', '--energy_loss', dest='loss_type', nargs='?', type=str, const='stochastic', default='stochastic', help='energy loss type for lepton - can be stochastic or continuous; default is stochastic')

parser.add_argument('-m', '--material', dest='material_id', nargs='?', type=str, const='rock', default='rock', help='material for energy loss - can be rock or water; default is rock')

parser.add_argument('-x', '--xc_model', dest='xc_model_id', nargs='?', type=str, const='cteq18_nlo', default='cteq18_nlo', help='neutrino cross-section model; default is CT18-NLO')

parser.add_argument('-p', '--pn_model', dest='pn_model_id', nargs='?', type=str, const='allm', default='allm', help='lepton photonuclear energy loss model; default is allm')

parser.add_argument('-f', '--fac_nu', dest='fac_nu_val', nargs='?', type=float, const=1.0, default=1.0, help='rescaling for SM cross-sections; default is 1.0')

parser.add_argument('-s', '--stat', dest='stat_val', nargs='?', type=float, const=1e7, default=1e7, help='statistics; default is 1e7')


args = parser.parse_args()
# energies = np.asarray([args.energy_val])

energies = args.energy_val
if type(energies) is str:
    energies = energies.split(",") if energies else []
    energies = np.asarray([float(i) for i in energies])

angles = args.angle_val
if type(angles) is str:
    angles = angles.split(",") if angles else []
    angles = np.asarray([int(i) for i in angles])

idepth = int(args.idepth_val)
lepton = str(args.lepton_id)
type_loss = str(args.loss_type)
material = str(args.material_id)
cross_section_model = str(args.xc_model_id)
pn_model = str(args.pn_model_id)
fac_nu = float(args.fac_nu_val)
stat = int(args.stat_val)

start_time = time.time()

prob_dict = Main.main(energies, angles, cross_section_model, pn_model, idepth, lepton, fac_nu, stat, type_loss)

e_out_files = glob.glob("e_out_*") # cleanup of Fortran e_out files

for file in e_out_files:
    os.remove(file)

end_time = time.time()
print(f"It took {end_time-start_time:.2f} seconds to compute")

# if __name__ == "__main__":
#     start_time = time.time()
#     # import main as Main
#     # Geometry.idepth = 7
#     prob_dict, e_out = Main.main(E_prop, angles, cross_section_model, pn_model, idepth, lepton, fac_nu, stat, prop_type)
#     print(prob_dict)
#     # print(Geometry.PREMdensity(1203))
#     end_time = time.time()
#     print(f"It took {end_time-start_time:.2f} seconds to compute")
