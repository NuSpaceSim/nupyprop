"""
Created on Wed July 31 12:01:03 2024

@author: Luke Kupari

"""

from nupyprop import propagation
import numpy as np


def single_stat(energy, angle, nu_xc, nu_ixc, depth, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,
                alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton, fac_nu,
                prop_type, u, w):
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
        u (_type_): Filename for Etau_out character size.
        w (_type_): Filename for Polarization_out character size.
        
    Returns:
        no_regen_tot (integer): No. of outgoing leptons without regeneration
        regen_tot (integer): No. of outgoing leptons with regeneration
    """
    
    e_format = "{:5.2f}"
    p_format = "{:8.5f}"
    no_regen_tot = 0
    regen_tot = 0
    
    depth0 = 0.0 #start with this each time
    
    #tnu goes until neutrino either goes to dtot, or converts to a tau
    
    ip,dtr,ef = propagation.propagate_nu(energy, nu_xc, nu_ixc, depth, fac_nu)
    
    #how far did the neutrino go? dtr is how far traveled
    
    depth0 = depth0 + dtr #how far is the neutrino on trajectory?
    
    dleft = depth - depth0 # how far is left for the neutrino to travel?
    
    if ip == 0:
        return no_regen_tot, regen_tot
    
    regen_cnt = 1 # tau out after first interaction
    
    etauin = ef
    #still need to propagate the tau, dolumn depth to go
    
    ipp, dfinal, etauf, Pi =  propagation.tau_thru_layers(angle, depth, dwater, depth0, etauin, xc_water, xc_rock, lep_ixc_water,
                                                         lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong,
                                                         idepth, lepton, prop_type)
    
    dleft = depth-dfinal
    
    with open(u, 'w') as e_file , open(w, 'w') as p_file:
    
        if ipp ==1 and dleft <= 0.0:
            Pout = Pi
            no_regen_tot = no_regen_tot + 1
            regen_tot = regen_tot + 1 #update the regen tau array once
            
            e_file.write(str(format.e_format(np.log10(etauf))) + '\n')
            p_file.write(str(format.p_format(Pout)) + '\n')
            
            return no_regen_tot,regen_tot # break outside stat; continue is correct here
        
        
        #must be a neutrino. Is there still column depth to propagate?
        
        ipp3 = 99 #dummy variable
        while dfinal < depthE and ipp3 != 1 and regen_cnt <=6:
            etauin = etauf # regen finds neutrino energy
            
            ipp3, dtau2, ef2, Pint = propagation.regen(angle, etauin, depth, dwater, dfinal, nu_xc, ithird, xc_water, xc_rock,
                                                       lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock,
                                                       xalong, cdalong, idepth, lepton, fac_nu, Pi)
        
            regen_cnt = regen_cnt + 1
            
            if ipp3 ==1: # then we are back toa tau at the end of the road
                regen_tot = regen_tot + 1
                Pout = Pint
                e_file.write(str(format.e_format(np.log10(ef2))) + '\n')
                p_file.write(str(format.p_format(Pout)) + '\n')
                
                return no_regen_tot, regen_tot
            
            if regen_cnt > 6:
                return no_regen_tot, regen_tot #only if regen > 6, break and go to run_stat for next iteraction
            
            etauf = ef2
            dfinal = dtau2
            Pi = Pint
            
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
    Efilename = 'eout_'+ str(format.format(np.log10(energy))) + '_ ' + str(angle) + '.dat'
    Pfilename = 'Pout_'+ str(format.format(np.log10(energy))) + '_ ' + str(angle) + '.dat'
    with open(Efilename, 'w') as e_file, open(Pfilename, 'w') as p_file:
        
        depth = depthE
        regen_cnt = 0
        no_regen_tot = 0
        regen_tot = 0
        
        for i in range(0,stats):
            no_regen_tot,regen_tot = single_stat(energy, angle, nu_xc, nu_ixc, depth, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,
            alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton, fac_nu, prop_type, 
            Efilename, Pfilename)
        
        
        e_file.close()
        p_file.close()
        return no_regen_tot,regen_tot