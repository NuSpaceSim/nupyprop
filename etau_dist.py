import numpy as np
import h5py 
import matplotlib.pyplot as plt


# Open the .h5 file in read mode
with h5py.File('output_nu_tau_4km_ct18nlo_allm_stochastic_1e5.h5', 'r') as f:
    for key in f.keys():
        print(key)
    dataset = f['CLep_out_energies']['10.0']['1.0'][:]
    dataset = np.array(dataset.astype(float))
    
plt.hist(dataset, bins=25, color='skyblue', edgecolor='black')
plt.xlabel('Lepton Energy')
plt.xscale('log')
plt.ylabel('Frequency')
plt.title('Tau Energy Distribution')
plt.grid(True)
plt.show()
