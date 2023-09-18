#/usr/bin/env python
from nupyprop import data
import numpy as np
import os
import shutil
from collections.abc import Iterable

energy = 10
#print((energy))

stats = 1e8

angles = np.arange(1,43,1)
#print(angles)
#angles=[1,2]

job_num = np.arange(0,100,1)
#job_num = [0,1]

path = str(energy) + "/"
files_path = "/home/dikgarg/Research/nupyprop/"

make_array = lambda x : x if isinstance(x, Iterable) else np.array([x]) # to avoid errors with single or no columns

'''for i in range(len(angles)):
    eout_arr=np.array([])
    Pout_arr=np.array([])
    pexit_noregen_arr=np.array([])
    pexit_regen_arr=np.array([])
    for job in job_num:
        eout = make_array(np.genfromtxt(("eout_{:.2f}_{:.1f}_{}.dat".format(energy, angles[i], job))))
        #print(eout)
        #print(str("eout_{:.2f}_{:.1f}_{}.dat".format(energy, angles[i], job)))

        eout_arr = np.append(eout_arr, eout)
        #os.remove("eout_{:.2f}_{:.1f}_{}.dat".format(energy, angles[i], job)) 
        
        Pout = make_array(np.genfromtxt(str("Pout_{:.2f}_{:.1f}_{}.dat".format(energy, angles[i], job))))  
        #print(str("Pout_{:.2f}_{:.1f}_{}.dat".format(energy, angles[i], job)))
        #print("Pout elements = ", Pout)

        Pout_arr = np.append(Pout_arr, Pout)
        #os.remove("Pout_{:.2f}_{:.1f}_{}.dat".format(energy, angles[i], job))
       
        pexit = make_array(np.genfromtxt(str("pexit_{:.2f}_{:.1f}_{}.dat".format(energy, angles[i], job)))) 
        #print("pexit vals = ", pexit[2], pexit[3])

        pexit_noregen_arr = np.append(pexit_noregen_arr, pexit[2])
        pexit_regen_arr = np.append(pexit_regen_arr, pexit[3])
        #os.remove("pexit_{:.2f}_{:.1f}_{}.dat".format(energy, angles[i], job))

    #print(eout_arr)
    with open ("eout_{:.2f}_{:.1f}.dat".format(energy, angles[i]), "a") as eout_file:
        for e in eout_arr:
            eout_file.write("%.5e \n" % (e))
    
    if Pout_arr.size==0:
        Pfinal = -1
    else:
        Pfinal = -1*np.mean(Pout_arr)

    #print(Pout_arr)
    #print(Pfinal)
    with open ("Pout_{:.2f}_{:.1f}.dat".format(energy, angles[i]), "a") as Pout_file:                    
        Pout_file.write("%.5e \n" % (Pfinal)) 
   
    prob_no_regen = np.sum(pexit_noregen_arr)
    prob_regen = np.sum(pexit_regen_arr)
    with open ("pexit_{:.2f}_{:.1f}.dat".format(energy, angles[i]), "a") as pexit_file:                  
        pexit_file.write("%.5e\t%.5e\t%.5e\t%.5e\n" % (10**(energy), angles[i], prob_no_regen/float(stats), prob_regen/float(stats)) ) 
    
    #print(prob_regen/float(stats)) 

for i in range(len(angles)):
    shutil.move(str("eout_{:.2f}_{:.1f}.dat".format(energy, angles[i])), path + str("eout_{:.2f}_{:.1f}.dat".format(energy, angles[i])))
    shutil.move(files_path + str("pexit_{:.2f}_{:.1f}.dat".format(energy, angles[i])), files_path + path + str("pexit_{:.2f}_{:.1f}.dat".format(energy, angles[i])))
    shutil.move(files_path + str("Pout_{:.2f}_{:.1f}.dat".format(energy, angles[i])), files_path + path + str("Pout_{:.2f}_{:.1f}.dat".format(energy, angles[i])))
'''
nu_type = 'neutrino' # nu for neutrino & anu for anti-neutrino
ch_lepton = 'tau' # type of charged lepton
cross_section_model = 'ct18nlo' # neutrino cross-section model
pn_model = 'bdhm' # photonuclear energy loss model
idepth = 4 # depth of water layer in km
stats = 1e8 # no. of ingoing neutrinos 
prop_type = 'stochastic' # type of energy loss; can be stochastic or continuous

data.process_htc_out(files_path, nu_type, ch_lepton, energy, idepth, cross_section_model, pn_model, prop_type, stats, elep_mode=True, arg=None)

'''
for i in range(len(angles)):
        os.remove(path + "eout_{:.2f}_{:.1f}.dat".format(energy, angles[i]))
        os.remove(path + "Pout_{:.2f}_{:.1f}.dat".format(energy, angles[i]))
        os.remove(path + "pexit_{:.2f}_{:.1f}.dat".format(energy, angles[i]))
'''
