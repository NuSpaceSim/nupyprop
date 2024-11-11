"""
Created on Wed July 31 12:01:03 2024

@author: Luke Kupari

"""

from nupyprop import propagation
import numpy as np
import matplotlib.pyplot as plt


def single_stat(energy, angle, nu_xc, nu_ixc, depth, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,
                alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton, fac_nu,
                prop_type, e_file, p_file):
    """Propagates a single ingoing neutrino event

    Args:
        energy (_type_): Incoming neutrino energy, in GeV.
        angle (_type_): Earth emergence angle (beta), in degrees.
        nu_xc (_type): 2D array containing neutrino CC & NC cross-section values, in cm^2.
        nu_ixc (_type_): 3D array containing neutrino integrated cross-section CDF values.
        depth (_type_): Maximum column depth for neutrino propagation, in kmwe.
        depthE (_type_): Total column depth for a given Earth emergence angle, in kmwe.
        dwater (_type_): Column depth along the chord for a given Earth emergence angle, in kmwe.
        xc_water (_type_): 2D array containing N_A/A*charged lepton-nucleon cross-section values in water, in cm^2/g.
        xc_rock (_type_): 2D array containing N_A/A*charged lepton-nucleon cross-section values in rock, in cm^2/g.
        lep_ixc_water (_type_): 3D array containing charged lepton integrated cross-section CDF values in water.
        lep_ixc_rock (_type_): 3D array containing charged lepton integrated cross-section CDF values in rock.
        alpha_water (_type_): 1D array containing ionization energy loss values in water, in (GeV*cm^2)/g.
        alpha_rock (_type_): 1D array containing ionization energy loss values in rock, in (GeV*cm^2)/g.
        beta_water (_type_): 2D array of beta values in water, in cm^2/g.
        beta_rock (_type_): 2D array of beta values in rock, in cm^2/g.
        xalong (_type_): 1D array containing distance in water, in km.
        cdalong (_type_): 1D array containing column depth at xalong, in g/cm^2.
        ithird (_type_): Choice for neutrino -> charged lepton energy fraction selection.
        idepth (_type_): Depth of water layer in km.
        lepton (_type_): Type of charged lepton. 1=tau; 2=muon.
        fac_nu (_type_): Rescaling factor for SM cross-sections.
        prop_type (_type_): Type of energy loss propagation. 1=stochastic, 2=continuous.
        u (_type_): File object for energies
        w (_type_): File object for Pout
        
    Returns:
        no_regen_tot (integer): No. of outgoing leptons without regeneration
        regen_tot (integer): No. of outgoing leptons with regeneration
        Pout (float):
        e_out (float):
    """
    
    e_format = "{:5.2f}"
    p_format = "{:8.5f}"
    regen_tot = 0
    no_regen_tot = 0 
    
    depth0 = 0.0 #start with this each time
    
    #tnu goes until neutrino either goes to dtot, or converts to a tau
    #print('depth= ', depth)
    ip,dtr,ef = propagation.propagate_nu(energy, nu_xc, nu_ixc, depth, fac_nu)
    
    #how far did the neutrino go? dtr is how far traveled
    
    depth0 = depth0 + dtr #how far is the neutrino on trajectory?
    
    dleft = depth - depth0 # how far is left for the neutrino to travel?

    if ip == 0: # still a neutrino at the end of the road
        #print('first go is a neutrino')
        return no_regen_tot, regen_tot
    
    regen_cnt = 1 # tau out after first interaction
    
    etauin = ef
    #still need to propagate the tau, dolumn depth to go
    ipp, dfinal, etauf, Pi =  propagation.tau_thru_layers(angle, depth, dwater, depth0, etauin, xc_water, xc_rock, lep_ixc_water,
                                                         lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong,
                                                        idepth, lepton, prop_type)
    
    #print('just propagated tau, ipp =', ipp, ' ipp=1: Tau, ipp=0: neutrino')
    dleft = depth-dfinal

    if ipp ==1 and dleft <= 0.0:
        Pout = Pi
        no_regen_tot = no_regen_tot + 1
        regen_tot = regen_tot + 1 #update the regen tau array once
        #print('no regeneration, tau makes it out')
        #print('Pout w/ no regen = ', Pout)
        e_file.write(str(e_format.format(np.log10(etauf))) + '\n')
        p_file.write(str(p_format.format(Pout)) + '\n')
        return no_regen_tot,regen_tot # break outside stat; continue is correct here
    
    
    #must be a neutrino. Is there still column depth to propagate?
    counter = 0
    ipp3 = 99 #dummy variable
    while dfinal < depthE and ipp3 != 1 and regen_cnt <=6:
        #print('in regeneration loop' , counter)
        etauin = etauf # regen finds neutrino energy
        
        ipp3, dtau2, ef2, Pint = propagation.regen(angle, etauin, depth, dwater, dfinal, nu_xc , nu_ixc, ithird, xc_water, xc_rock,
                                                    lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock,
                                                    xalong, cdalong, idepth, lepton, fac_nu, prop_type ,Pi)
        #print('after regeneration, Pint = ', Pint)
        #print('ipp3 =', ipp3, ' ipp3 =1: tau; ipp3 = 0 neutrino')
        
        regen_cnt = regen_cnt + 1
        if ipp3 ==1: # then we are back to a tau at the end of the road
            regen_tot = regen_tot + 1
            Pout = Pint
            #print('back to tau after regenewrated neutrino, exiting, P = ', Pout)
            e_file.write(str(e_format.format(np.log10(etauf))) + '\n')
            p_file.write(str(p_format.format(Pout)) + '\n')
            return no_regen_tot, regen_tot
        if regen_cnt > 6:
            #print('regeneration exceeded 6')
            return no_regen_tot, regen_tot #only if regen > 6, break and go to run_stat for next iteraction
        
        etauf = ef2
        dfinal = dtau2
        Pi = Pint
        counter = counter + 1
    
    return no_regen_tot,regen_tot
            
def run_stat_single(energy, angle, nu_xc, nu_ixc, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,
                    alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton, fac_nu,
                    stats, prop_type):
    """run a loop for all ingoing neutrinos

    Args:
        energy (float): Incoming neutrino energy, in GeV.
        angle (float): Earth emergence angle (beta), in degrees.
        nu_xc (np.ndarray): 2D array containing neutrino CC & NC cross-section values, in cm^2.
        nu_ixc (np.ndarray): 3D array containing neutrino integrated cross-section CDF values.
        depthE (float): Total column depth for a given Earth emergence angle, in kmwe.
        dwater (float): Column depth along the chord for a given Earth emergence angle, in kmwe.
        xc_water (np.ndarray): 2D array containing N_A/A*charged lepton-nucleon cross-section values in water, in cm^2/g.
        xc_rock (np.ndarray): 2D array containing N_A/A*charged lepton-nucleon cross-section values in rock, in cm^2/g.
        lep_ixc_water (np.ndarray): 3D array containing charged lepton integrated cross-section CDF values in water.
        lep_ixc_rock (np.ndarray): 3D array containing lepton integrated cross-section CDF values in rock.
        alpha_water (np.ndarray): 1D array containing ionization energy loss values in water, in (GeV*cm^2)/g.
        alpha_rock (np.ndarray): 1D array containing ionization energy loss values in rock, in (GeV*cm^2)/g.
        beta_water (np.ndarray): 2D array of beta values in water, in cm^2/g.
        beta_rock (np.ndarray): 2D array of beta values in rock, in cm^2/g.
        xalong (np.ndarray): 1D array containing distance in water, in km.
        cdalong (np.ndarray): 1D array containing column depth at xalong, in g/cm^2.
        ithird (integer): Choice for neutrino -> charged lepton energy fraction selection.
        idepth (integer): Depth of water layer in km.
        lepton (integer): Type of charged lepton. 1=tau; 2=muon.
        fac_nu (float): Rescaling factor for SM neutrino cross-sections.
        stats (integer): Statistics or no. of ingoing neutrinos.
        prop_type (integer): Type of energy loss propagation. 1=stochastic, 2=continuous.
    
    Returns:
        no_regen_tot (integer): No. of outgoing charged leptons without regeneration.
        regen_tot (integer): No. of outgoing charged leptons with regeneration.
    """
    
    format = "{:.2f}"
    Efilename = 'eout_'+ str(format.format(np.log10(energy))) + '_' + str(angle) + '.dat'
    Pfilename = 'Pout_'+ str(format.format(np.log10(energy))) + '_' + str(angle) + '.dat'  
    no_regen_tot = 0
    regen_tot = 0
    '''
    with open(Efilename, 'a') as e_file, open(Pfilename, 'a') as p_file:
        depth = depthE
        
        for i in range(0,stats):
            tempnrt, temprt = single_stat(energy, angle, nu_xc, nu_ixc, depth, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,
            alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton, fac_nu, prop_type, 
            e_file, p_file)
            no_regen_tot = no_regen_tot + tempnrt
            regen_tot = regen_tot + temprt
            
    e_file.close()
    p_file.close()
    '''

    #PROPAGATE_NU DEBUGGING BLOCK
    part_type = []
    d_travel = []
    e_fin = []
    #print('depthE' , depthE)
    for i in range(0,stats):
        ##print(i)
        result = propagation.propagate_nu(energy,nu_xc,nu_ixc,depthE,fac_nu)
        part_type.append(result[0])
        d_travel.append(result[1])
        e_fin.append(result[2])
            
    
    no_regen_tot =0
    regen_tot = 0
    part_type = np.array(part_type)
    d_travel = np.array(d_travel)
    e_fin = np.array(e_fin)
    e_init = np.zeros_like(e_fin)
    e_init[:] = energy

    bin_edges = np.arange(0, d_travel.max() + 100, 100)
    
    count_ones = np.count_nonzero(part_type == 1)
    count_zeros = np.count_nonzero(part_type == 0)
    
    print('1s: ', count_ones,'0s: ', count_zeros)
    print(np.mean(d_travel))
    print(np.mean(e_fin),np.min(e_fin),np.max(e_fin))
    
    fig, axs = plt.subplots(2,2, figsize=(10,8))
    axs[0,0].hist(part_type, bins =[-0.5, 0.5, 1.5], edgecolor='black', rwidth=0.7)
    axs[0,0].set_title('propagate_nu particle type')
    axs[0,0].set_xlabel('0: neutrino or 1: charged lepton ')
    #axs[0,0].set_ylabel('Normalized Counts')

    axs[0,1].hist(d_travel, bins=bin_edges, edgecolor = 'black')
    axs[0,1].set_title('propagate_nu dtravel')
    axs[0,1].set_xlabel('dtravel value')
    axs[0,1].set_yscale('log')
    #axs[0,1].set_ylabel('Normalized counts')
    
    axs[1,0].hist(e_init, bins = 3, edgecolor = 'black')
    axs[1,0].set_title('propagate_nu initial energy')
    axs[1,0].set_xlabel('e_init')
    #axs[1,0].set_ylabel(' counts')
    
    axs[1,1].hist(e_fin, bins = 20, edgecolor = 'black')
    axs[1,1].set_title('propagate_nu e_fin')
    axs[1,1].set_xlabel('e_fin')
    #axs[1,1].set_ylabel(' counts')
    
    plt.suptitle('propagate_nu output for 10^10 GeV, 30 deg, 1e7 stats')
    plt.tight_layout()
    plt.show()

    '''
    ### PROPAGATE_LEP_WATER BLOCK
    no_regen_tot = 1
    regen_tot = 1
    part_id = []
    d_fin = []
    e_fin = []
    pcthf = []
    
    for i in range(0,stats):
        results = propagation.propagate_lep_water(energy,xc_water,lep_ixc_water,alpha_water,beta_water,dwater,lepton,                                       prop_type,cthi=.25,Pi=1)
        part_id.append(results[0])
        d_fin.append(results[1])
        e_fin.append(results[2])
        pcthf.append(results[3])
    
    count_ones = part_id.count(1)
    count_zeros = part_id.count(0)
    #print(part_id)
    #print('1s: ', count_ones,'0s: ', count_zeros)
    #print(np.mean(d_fin))
    #print(np.mean(e_fin),np.min(e_fin),np.max(e_fin))
    #print(np.mean(pcthf))
    ### PROPAGATE_LEP_ROCK BLOCK
    no_regen_tot = 1
    regen_tot = 1
    part_id = []
    d_fin = []
    e_fin = []
    cthf = []
    Pf = []
    depth_traj = depthE - idepth
    for i in range(0,stats):
        results = propagation.propagate_lep_rock(angle ,energy, xc_rock,lep_ixc_rock,
                                                 alpha_rock, beta_rock,depthE,depth_traj,
                                                 xalong,cdalong,idepth,lepton, prop_type)
    '''
    '''
    ###TAU_THRU_LAYERS BLOCK
    depth0 = 0.0 #start with this each time
    
    #tnu goes until neutrino either goes to dtot, or converts to a tau
    
    ip,dtr,ef = propagation.propagate_nu(energy, nu_xc, nu_ixc, depthE, fac_nu)
    
    #how far did the neutrino go? dtr is how far traveled
    
    depth0 = depth0 + dtr #how far is the neutrino on trajectory?
    
    dleft = depthE - depth0 # how far is left for the neutrino to travel?
    
        
    regen_cnt = 1 # tau out after first interaction
    
    etauin = ef
    #still need to propagate the tau, dolumn depth to go
    ipp, dfinal, etauf, Pi =  propagation.tau_thru_layers(angle, depthE, dwater, depth0, etauin, xc_water, xc_rock, lep_ixc_water,
                                                         lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong,
                                                         idepth, lepton, prop_type)
    
    dleft = depthE-dfinal
      
    no_regen_tot=1 
    regen_tot = 1
    '''
    print('no_regen', no_regen_tot)
    print('regen', regen_tot)
    return no_regen_tot,regen_tot
