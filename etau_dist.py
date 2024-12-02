import numpy as np
import h5py 
import matplotlib.pyplot as plt

dataset = []
angles = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
angles_longctau = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,1.0]
print(angles[2])

# Open the .h5 file in read mode
with h5py.File('output_nu_tau_4km_ct18nlo_allm_stochastic_1e7.h5', 'r') as f:
    for key in f.keys():
        print(key)
    
    for angle in angles_longctau:
        grab = f['CLep_out_energies']['9.0'][str(angle)][:]
        dataset.append(np.array(grab.astype(float)))

#omin = min(min(data) for data in dataset)
#omax = max(max(data) for data in dataset)

#bin_num = 15
#bin_edges = np.linspace(omin,omax,bin_num)


for i in range(len(dataset)):
    data = np.array(dataset[i])
    data = 10**data
    minE = min(data)
    maxE = max(data)
    bins = np.logspace(np.log10(minE),np.log10(maxE), 50)
    print(bins)
    hist= np.histogram(data, bins = bins)
    #bin_centers = 0.5 * (bins[:-1] + bins[1:])
    #print(bin_centers)
    
    plt.hist(data,bins=bins, edgecolor= 'black', alpha = 0.75, label=f"{angles[i]}")
    #plt.scatter(bin_centers, hist, s=20,alpha = 0.75,label=f"{angles_longctau[i]}") 
    plt.grid(True, alpha=0.5)
    plt.xlabel('Tau Energy (GeV)')
    plt.xscale('log')
    #plt.yscale('log')
    plt.ylabel('frequency')
    plt.title(' Tau Energy Distribution')
    plt.legend()
    plt.show()
    plt.clf()
    #plt.savefig('EDTau' + str(angles_longctau[i]) + 'longctau_.pdf', format='pdf')
    #plt.clf()
    
#plt.legend()
#plt.show()

'''
plot_number = 0
plt.hist(dataset[plot_number], bins=25, color='skyblue', edgecolor='black',alpha=0.5,label=angles[plot_number])
plt.xlabel('Lepton Energy')
plt.xscale('log')
plt.ylabel('Frequency')
plt.title('Tau Energy Distribution')
plt.grid(True)
plt.legend()
plt.show()
'''