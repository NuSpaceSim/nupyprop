"""
Created on Wed July 31 12:01:03 2024

@author: Luke Kupari and Diksha Garg
Comments: This code calls functions from propagation routine to 
propagate neutrinos and charged leptons.
"""

from nupyprop import propagation
import numpy as np
#import matplotlib.pyplot as plt
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

def run_stat(energy, angle, nu_xc, nu_ixc, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,
                alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton, fac_nu,
                prop_type, stats, earth_model): 
    '''
    Propagates a single ingoing neutrino event

    Parameters
    ----------
    energy : float 
        Incoming neutrino energy, in GeV.
    angle : float
        Earth emergence angle (beta), in degrees.
    nu_xc : np.ndarray 
        2D array containing neutrino CC & NC cross-section values, in cm^2.
    nu_ixc : np.ndarray 
        3D array containing neutrino integrated cross-section CDF values.
    depthE : float
        Total column depth for a given Earth emergence angle, in kmwe.
    dwater : float
        Column depth along the chord for a given Earth emergence angle, in kmwe.
    xc_water : np.ndarray 
        2D array containing N_A/A*charged lepton-nucleon cross-section values in water, in cm^2/g.
    xc_rock : np.ndarray 
        2D array containing N_A/A*charged lepton-nucleon cross-section values in rock, in cm^2/g.
    lep_ixc_water : np.ndarray 
        3D array containing charged lepton integrated cross-section CDF values in water.
    lep_ixc_rock : np.ndarray 
        3D array containing charged lepton integrated cross-section CDF values in rock.
    alpha_water : np.ndarray 
        1D array containing ionization energy loss values in water, in (GeV*cm^2)/g.
    alpha_rock : np.ndarray 
        1D array containing ionization energy loss values in rock, in (GeV*cm^2)/g.
    beta_water : np.ndarray 
        2D array of beta values in water, in cm^2/g.
    beta_rock : np.ndarray 
        2D array of beta values in rock, in cm^2/g.
    xalong : np.ndarray 
        1D array containing distance in water, in km.
    cdalong : np.ndarray  
        1D array containing column depth at xalong, in g/cm^2.
    ithird : int
        Choice for neutrino -> charged lepton energy fraction selection.
    idepth : int
        Depth of water layer in km.
    lepton : int 
        Type of charged lepton. 1=tau; 2=muon.
    fac_nu : float
        Rescaling factor for SM cross-sections.
    prop_type : int
        Type of energy loss propagation. 1=stochastic, 2=continuous.
    stats : integer 
        Statistics or no. of ingoing neutrinos.
    earth_model : str
        Earth density model, prem or ak135.

    Returns (stats_cl : number of charged leptons produced after first neutrino interaction)
    -------
    no_regen_tot : np.ndarray (stats_cl,)
        1 if event produced an exiting lepton with NO regeneration, else 0.
    regen_tot    : np.ndarray (stats_cl,)
        1 if event produced an exiting lepton WITH regeneration, else 0.
    e_out        : np.ndarray (stats_cl,)
        Final exiting lepton energy (GeV); 0 for events with no exiting lepton.
    P_out        : np.ndarray (stats_cl,)
        Final exiting lepton polarization; 0 for events with no exiting lepton.
    '''

    # --- Create per-particle arrays ---
    energy_arr = np.full(stats, energy)       # all start with same initial energy 
    depthE_arr = np.full(stats, depthE)       # total column depth for given angle (kmwe)

    ############# propagate neutrinos #############
    _, dlep, elep = propagation.propagate_nu(energy_arr, nu_xc, nu_ixc, depthE_arr, fac_nu, stats, 
                                           Emin, E_nu, E_lep, yvals)

    stats_cl = len(dlep) # number of charged leptons
    if stats_cl == 0:
        # no charged leptons produced at all
        return (np.zeros(0, dtype=int),
                np.zeros(0, dtype=int),
                np.zeros(0, dtype=float),
                np.zeros(0, dtype=float))

    depthE_arr = np.full(stats_cl, depthE) # adjusting len of array, tot col depth
    dwater_arr = np.full(stats_cl, dwater)       # total water column depth for given idepth (kmwe) 

    # --- state for propagation --- #
    active = np.ones(stats_cl, dtype=bool)  # track which leptons still propagate
    regen_count = np.ones(stats_cl, dtype=int) # charged lepton out after first interaction

    ###### initializing the final output arrays ######
    e_out = np.zeros(stats_cl, dtype=float)
    P_out = np.zeros(stats_cl, dtype=float)
    no_regen_tot, regen_tot = np.zeros(stats_cl, dtype=int), np.zeros(stats_cl, dtype=int)

    ############# now propagate charged lepons #############
    ipp, dfinal, etauf, Pi =  propagation.tau_thru_layers(angle, depthE_arr, dwater_arr, dlep, 
                                                          elep, xc_water, xc_rock, lep_ixc_water,
                                                          lep_ixc_rock, alpha_water, alpha_rock, beta_water,
                                                          beta_rock, xalong, cdalong, idepth, lepton, 
                                                          prop_type, Emin, E_nu, E_lep, yvals, ypol, 
                                                          Pcthp, P, earth_model, stats_cl)
    
    dleft = depthE_arr - dfinal # col depth left for charged lepton to travel, kmwe
    no_regen_mask = (ipp == 1) & (dleft <= 0.0) # cl that exit without regenerating

    e_out[no_regen_mask]        = etauf[no_regen_mask]
    P_out[no_regen_mask]        = Pi[no_regen_mask]
    no_regen_tot[no_regen_mask] = 1
    regen_tot[no_regen_mask]    = 1   # same as your scalar code
    
    active[no_regen_mask] = False 

    if not np.any(active):
        return no_regen_tot, regen_tot, e_out, P_out

    # ----------------------- Regeneration Loop ----------------------- #
    max_regen = 6

    # Continue while *any* lepton is still active
    while np.any(active):
        # eligible for regeneration
        elig = active & (dfinal < depthE) & (regen_count <= max_regen)
        if not np.any(elig):
            break

        idx = np.where(elig)[0]

        ipp3, dtau2, ef2, Pint = propagation.regen(
            angle, etauf[idx], depthE, dwater, dfinal[idx], nu_xc, 
            nu_ixc, ithird, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,
            alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong,
            idepth, lepton, fac_nu, prop_type, Pi[idx], Emin, E_nu, E_lep,
            yvals, ypol, Pcthp, P, earth_model)

        # Update regeneration counter
        regen_count[idx] += 1

        # ------------------ Case 1: these regenerated into taus (ipp3 == 1) ------------------
        exit_mask = (ipp3 == 1)

        global_exit_idx = idx[exit_mask]     # map back to global indexing

        regen_tot[global_exit_idx] = 1
        e_out[global_exit_idx]     = ef2[exit_mask]
        P_out[global_exit_idx]     = Pint[exit_mask]

        # deactivate these leptons
        active[global_exit_idx] = False

        # ------------------ Case 2: continue propagation for neutrinos ------------------
        cont_mask = ~exit_mask 
        cont_idx = idx[cont_mask]

        etauf[cont_idx]   = ef2[cont_mask]
        dfinal[cont_idx] = dtau2[cont_mask]
        Pi[cont_idx]     = Pint[cont_mask]

        # Regen_count > 6 automatically deactivated next iteration

    return no_regen_tot, regen_tot, e_out, P_out

'''part_id, d_fin, e_fin, cthf, Pf = propagation.propagate_lep(energy_arr, angle, xc_water, lep_ixc_water, alpha_water, beta_water, depth0_arr, dwater_arr, lepton, 
                                                            prop_type,
                                                            'water', # either 'water' or 'rock'
                                                            cthi_arr, Pi_arr,
                                                            idepth, earth_model, Emin, E_nu, E_lep, yvals, ypol, Pcthp, P, xalong=None, cdalong=None)'''

'''d_in = depthE - depth0 - dwater #propagate this far in rock
print("d_in = ", d_in)
din_arr = np.full(stats, d_in) 
part_id, d_fin, e_fin, cthf, Pf = propagation.propagate_lep(energy_arr, angle, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, depth0_arr, din_arr, lepton, prop_type,
                                                            'rock', # either 'water' or 'rock'
                                                            cthi_arr, Pi_arr,
                                                            idepth, earth_model, Emin, E_nu, E_lep, yvals, ypol, Pcthp, P, xalong, cdalong)'''

# def run_stat_single(energy, angle, nu_xc, nu_ixc, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,
#                     alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton, fac_nu,
#                     stats, prop_type, earth_model):
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

    # format = "{:.2f}"
    # Efilename = 'eout_'+ str(format.format(np.log10(energy))) + '_' + str(angle) + '.dat'
    # Pfilename = 'Pout_'+ str(format.format(np.log10(energy))) + '_' + str(angle) + '.dat'
    # no_regen_tot = 0
    # regen_tot = 0
    # with open(Efilename, 'a') as e_file, open(Pfilename, 'a') as p_file:
    #     #iparr, dtrarr, efarr = [], [], []
    #     '''for i in tqdm(range(batch_num)):
    #         tempnrt, temprt, ef = single_stat(energy, angle, nu_xc, nu_ixc, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,
    #         alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton, fac_nu, prop_type, int(stats/batch_num),
    #         e_file, p_file)
    #         # iparr.append(tempnrt)
    #         # dtrarr.append(temprt)
    #         # efarr.append(ef)
    #         #no_regen_tot = no_regen_tot + tempnrt
    #         #regen_tot = regen_tot + temprt'''

    #     part_id, d_fin, e_fin, cthf, Pf = single_stat(energy, angle, nu_xc, nu_ixc, depthE, dwater, 
    #                                                   xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,
    #                                                   alpha_water, alpha_rock, beta_water, beta_rock,
    #                                                   xalong, cdalong, ithird, idepth, lepton, fac_nu,
    #                                                   prop_type, stats, earth_model)

    # e_file.close()
    # p_file.close()

    # plt.hist(part_id, 50)
    # plt.xlabel("id")
    # plt.yscale('log')
    # plt.title(f"Angle={angle}")
    # plt.savefig(f"1e{np.log10(energy)}GeV_{angle}deg_part_id.png")
    # plt.show()

    # plt.hist(d_fin, 50)
    # plt.xlabel("dfinal")
    # plt.yscale('log')
    # plt.title(f"Angle={angle}")
    # plt.savefig(f"1e{np.log10(energy)}GeV_{angle}deg_dfinal.png")
    # plt.show()

    # plt.hist(e_fin, 20)
    # plt.xlabel("efinal")
    # plt.title(f"Angle={angle}")
    # plt.loglog()
    # plt.savefig(f"1e{np.log10(energy)}GeV_{angle}deg_efinal.png")
    # plt.show()

    # plt.hist(Pf*cthf, 50)
    # plt.xlabel("polarization")
    # plt.title(f"Angle={angle}")
    # plt.semilogy()
    # plt.savefig(f"1e{np.log10(energy)}GeV_{angle}deg_polarization.png")
    # plt.show()

    # return no_regen_tot,regen_tot
