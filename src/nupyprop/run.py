"""
Created on Wed July 31 12:01:03 2024

@author: Luke Kupari

"""

from nupyprop import propagation
import numpy as np
import matplotlib.pyplot as plt
import nupyprop.constants as const

batch_num = const.batch_num #batch size to divide the total stats

Emin = const.Emin #min threshold energy for leptons

rho_rock = const.rho_rock # rock density
rho_iron = const.rho_iron # iron density

E_nu = const.E_nu # Neutrino energy numpy array
E_lep = const.E_lep # Lepton energy numpy array
yvals = const.yvals # inelasticity from 1e-3 to 1

#variables for tau-lepton polarization
ypol, Pcthp, P = const.ypol, const.Pcthp, const.P

def single_stat(energy, angle, nu_xc, nu_ixc, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,
                alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton, fac_nu,
                prop_type, stats, earth_model): #, e_file, p_file):
    '''
    Propagates a single ingoing neutrino event

    Args:
        energy (_type_): Incoming neutrino energy, in GeV.
        angle (_type_): Earth emergence angle (beta), in degrees.
        nu_xc (_type): 2D array containing neutrino CC & NC cross-section values, in cm^2.
        nu_ixc (_type_): 3D array containing neutrino integrated cross-section CDF values.
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
        stats (integer): Statistics or no. of ingoing neutrinos.
        earth_model (str): Earth density model, prem or ak135.
        u (_type_): File object for energies
        w (_type_): File object for Pout

    Returns:
        no_regen_tot (integer): No. of outgoing leptons without regeneration
        regen_tot (integer): No. of outgoing leptons with regeneration
        Pout (float):
        e_out (float):
    '''

    e_format = "{:5.2f}"
    p_format = "{:8.5f}"
    regen_tot = 0
    no_regen_tot = 0

    depth0 = 0.0 #start with this each time
    
    # --- Create per-particle arrays ---
    energy_arr = np.full(stats, energy)       # all start with same initial energy
    depth0_arr = np.zeros(stats)              # start column depth for each neutrino 
    dwater_arr = np.full(stats, dwater)       # total water column depth for given idepth (kmwe) 
    depthE_arr = np.full(stats, depthE)       # total column depth for given angle (kmwe)
    Pi_arr = np.ones(stats)                   # all initial polarization = 1
    cthi_arr = np.cos(np.zeros(stats))        # all initial costheta = 1
    
    '''part_id, d_fin, e_fin, cthf, Pf = propagation.propagate_lep(energy_arr, angle, xc_water, lep_ixc_water, alpha_water, beta_water, depth0_arr, dwater_arr, lepton, 
                                                                prop_type,
                                                               'water', # either 'water' or 'rock'
                                                               cthi_arr, Pi_arr,
                                                               idepth, earth_model, Emin, E_nu, E_lep, yvals, ypol, Pcthp, P, xalong=None, cdalong=None)'''
    
    d_in = depthE - depth0 - dwater #propagate this far in rock
    print("d_in = ", d_in)
    din_arr = np.full(stats, d_in) 
    part_id, d_fin, e_fin, cthf, Pf = propagation.propagate_lep(energy_arr, angle, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, depth0_arr, din_arr, lepton, prop_type,
                                                               'rock', # either 'water' or 'rock'
                                                               cthi_arr, Pi_arr,
                                                               idepth, earth_model, Emin, E_nu, E_lep, yvals, ypol, Pcthp, P, xalong, cdalong)
    
    '''ip, dtr, ef = propagation.propagate_nu(energy_arr, nu_xc, nu_ixc, depthE_arr, fac_nu, stats, Emin, E_nu, E_lep, yvals)

    return ip, dtr, ef'''

    #tnu goes until neutrino either goes to dtot, or converts to a tau
    #print('depth= ', depth)
    '''ip, dtr, ef = propagation.propagate_nu(energy, nu_xc, nu_ixc, depthE, fac_nu, stats, Emin, E_nu, E_lep, yvals)

    #how far did the neutrino go? dtr is how far traveled

    depth0 = depth0 + dtr #how far is the neutrino on trajectory?

    dleft = depthE - depth0 # how far is left for the neutrino to travel?

    if ip == 0: # still a neutrino at the end of the road
        #print('first go is a neutrino')
        return no_regen_tot, regen_tot

    regen_cnt = 1 # tau out after first interaction

    etauin = ef
    #still need to propagate the tau, dolumn depth to go
    ipp, dfinal, etauf, Pi =  propagation.tau_thru_layers(angle, depthE, dwater, depth0, etauin, xc_water, xc_rock, lep_ixc_water,
                                                         lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong,
                                                        idepth, lepton, prop_type, Emin, E_nu, E_lep, yvals, ypol, Pcthp, P, earth_model)

    #print('just propagated tau, ipp =', ipp, ' ipp=1: Tau, ipp=0: neutrino')
    dleft = depthE-dfinal

    if ipp ==1 and dleft <= 0.0:
        Pout = Pi
        no_regen_tot = no_regen_tot + 1
        regen_tot = regen_tot + 1 #update the regen tau array once
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

        ipp3, dtau2, ef2, Pint = propagation.regen(angle, etauin, depthE, dwater, dfinal, nu_xc , nu_ixc, ithird, xc_water, xc_rock,
                                                    lep_ixc_water, lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock,
                                                    xalong, cdalong, idepth, lepton, fac_nu, prop_type ,Pi, Emin, E_nu, E_lep, yvals, ypol, Pcthp, P,
                                                    earth_model)

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
            print('regeneration exceeded 6')
            return no_regen_tot, regen_tot #only if regen > 6, break and go to run_stat for next iteraction

        etauf = ef2
        dfinal = dtau2
        Pi = Pint
        counter = counter + 1

    return no_regen_tot,regen_tot'''
    return part_id, d_fin, e_fin, cthf, Pf

def run_stat_single(energy, angle, nu_xc, nu_ixc, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,
                    alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton, fac_nu,
                    stats, prop_type, earth_model):
    '''
    Run a loop for all ingoing neutrinos

    Parameters
    ----------
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
    earth_model (str): Earth density model, prem or ak135.

    Returns
    --------
    no_regen_tot (integer): No. of outgoing charged leptons without regeneration.
    regen_tot (integer): No. of outgoing charged leptons with regeneration.
    '''

    format = "{:.2f}"
    Efilename = 'eout_'+ str(format.format(np.log10(energy))) + '_' + str(angle) + '.dat'
    Pfilename = 'Pout_'+ str(format.format(np.log10(energy))) + '_' + str(angle) + '.dat'
    no_regen_tot = 0
    regen_tot = 0
    with open(Efilename, 'a') as e_file, open(Pfilename, 'a') as p_file:
        #iparr, dtrarr, efarr = [], [], []
        '''for i in tqdm(range(batch_num)):
            tempnrt, temprt, ef = single_stat(energy, angle, nu_xc, nu_ixc, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,
            alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton, fac_nu, prop_type, int(stats/batch_num),
            e_file, p_file)
            # iparr.append(tempnrt)
            # dtrarr.append(temprt)
            # efarr.append(ef)
            #no_regen_tot = no_regen_tot + tempnrt
            #regen_tot = regen_tot + temprt'''

        part_id, d_fin, e_fin, cthf, Pf = single_stat(energy, angle, nu_xc, nu_ixc, depthE, dwater, 
                                                      xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,
                                                      alpha_water, alpha_rock, beta_water, beta_rock,
                                                      xalong, cdalong, ithird, idepth, lepton, fac_nu,
                                                      prop_type, stats, earth_model)

    e_file.close()
    p_file.close()

    plt.hist(part_id, 50)
    plt.xlabel("id")
    plt.yscale('log')
    plt.title(f"Angle={angle}")
    plt.savefig(f"1e{np.log10(energy)}GeV_{angle}deg_part_id.png")
    plt.show()

    plt.hist(d_fin, 50)
    plt.xlabel("dfinal")
    plt.yscale('log')
    plt.title(f"Angle={angle}")
    plt.savefig(f"1e{np.log10(energy)}GeV_{angle}deg_dfinal.png")
    plt.show()

    plt.hist(e_fin, 20)
    plt.xlabel("efinal")
    plt.title(f"Angle={angle}")
    plt.loglog()
    plt.savefig(f"1e{np.log10(energy)}GeV_{angle}deg_efinal.png")
    plt.show()

    plt.hist(Pf*cthf, 50)
    plt.xlabel("polarization")
    plt.title(f"Angle={angle}")
    plt.semilogy()
    plt.savefig(f"1e{np.log10(energy)}GeV_{angle}deg_polarization.png")
    plt.show()

    return no_regen_tot,regen_tot
