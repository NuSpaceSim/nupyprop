############################################################
#date: April 17, 2024
#author: Diksha Garg
#comment: This file contains all the required constants for nupyprop
############################################################
import numpy as np
import scipy.constants as scc

#useful constants
N_A = 6.0221409e+23 # Avogadro's number in units of 1/mole

R_earth = 6378.137 #radius of Earth in km

rho_water, rho_rock, rho_iron = 1.02, 2.65, 7.87 #density of water, rock and iron in g/cm^3

m_mu = scc.physical_constants["muon mass energy equivalent in MeV"][0]*1e-3 # muon mass in GeV
m_tau = scc.physical_constants["tau mass energy equivalent in MeV"][0]*1e-3 # tau mass in GeV
m_pi = 139.57018e-3 # pion mass in GeV
m_p = scc.physical_constants["proton mass energy equivalent in MeV"][0]*1e-3 # proton mass in GeV

alpha_fs = scc.fine_structure #fine structure constant

#Earth emergence angles
beta_arr = np.asarray([float('{:.1f}'.format(i)) for i in np.arange(0.1,90.1,step=0.1)]) # finer steps of 0.1 degerees

#Cross-section models that exsist in src/nupyprop/datafiles/lookup_table.h5
nu_models = ['allm', 'bdhm', 'ct18nlo', 'nct15'] #for weak interaction
pn_models = ['allm', 'bb'] #for Photo-nuclear electromagnetic interaction

#step size for continuous energy loss
step_size = 4500 #in cm

#array of energies used in interpolation for cross-section and CDFs
E_nu = np.logspace(3,12,91,base=10).astype(np.float64) #for neutrinos in GeV
E_lep = np.logspace(0,12,121,base=10).astype(np.float64) #for charged leptons in GeV

#inelasticity values used in interpolation for energy CDFs for neutrinos and charged leptons
yvals = np.logspace(-3,0,31).astype(np.float64) #goes from 1e-3 to 1

#Earth layers considered in PREM - Earth density model
Rlay = np.array([1221.5, 3480.0, 5701.0, 5771.0, 5971.0, 6151.0, 6346.6, 6356.0, 6368.0, 6371.0])


