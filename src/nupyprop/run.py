"""
Created on Wed July 31 12:01:03 2024

@author: Luke Kupari and Diksha Garg
Comments: This code calls functions from propagation routine to 
propagate neutrinos and charged leptons.
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

def run_stat(energy, angle, nu_xc, nu_ixc, nu_bsm_xc, nu_bsm_ixc, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,
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
    nu_bsm_xc : np.ndarray 
        2D array containing neutrino BSM CC & NC cross-section values, in cm^2.
    nu_bsm_ixc : np.ndarray 
        3D array containing BSM neutrino integrated cross-section CDF values.
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
#    print('length of energy array before propagate_nu',len(energy_arr),len(depthE_arr))
    ############# propagate neutrinos #############
    id, dlep, elep = propagation.propagate_nu(energy_arr, nu_xc, nu_ixc, nu_bsm_xc, nu_bsm_ixc, depthE_arr, fac_nu, stats, 
                                           Emin, E_nu, E_lep, yvals)
    
    #print('length after propagate_nu', len(elep),len(dlep))
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
    
    active[no_regen_mask] = False 
    if not np.any(active):
        return no_regen_tot, regen_tot, e_out, P_out

    # ----------------------- Regeneration Loop ----------------------- #
    max_regen = 6
    
    while np.any(active):
        # eligible for regeneration: taus that did not exit, still have depth remaining,
        # and have not exceeded the regeneration limit
        elig = active & (ipp != 1) & (dfinal < depthE_arr) & (regen_count <= max_regen)
        if not np.any(elig):
            break

        cand = np.where(elig)[0]  # candidates (size m)
        regen_count[cand] += 1
        #print("regen: cand", cand.size)

        kept_local, ipp3, dtau2, ef2, Pint = propagation.regen_compact(
            angle,
            etauf[cand],
            depthE_arr[cand],
            dwater_arr[cand],
            dfinal[cand],
            nu_xc,
            nu_ixc,
            ithird,
            xc_water,
            xc_rock,
            lep_ixc_water,
            lep_ixc_rock,
            alpha_water,
            alpha_rock,
            beta_water,
            beta_rock,
            xalong,
            cdalong,
            idepth,
            lepton,
            fac_nu,
            prop_type,
            Pi[cand],
            Emin,
            E_nu,
            E_lep,
            yvals,
            ypol,
            Pcthp,
            P,
            earth_model,
        )
        #print("regen: kept_local", kept_local.size)
        #print("regen: exiting after regen", np.sum((ipp3 == 1) & ((depthE_arr[cand[kept_local]] - dtau2) <= 0.0)))
        # If nothing regenerated into a tau via CC, all candidates die as neutrinos
        if kept_local.size == 0:
            active[cand] = False
            continue

        surv = cand[kept_local]  # survivors mapped back to stats_cl indexing (size k)

        # Exiting taus after regen
        exit_mask = (ipp3 == 1) & ((depthE_arr[surv] - dtau2) <= 0.0)
        if np.any(exit_mask):
            exit_idx = surv[exit_mask]
            regen_tot[exit_idx] = 1
            e_out[exit_idx] = ef2[exit_mask]
            P_out[exit_idx] = Pint[exit_mask]
            active[exit_idx] = False

        # Continued taus (decayed again inside Earth) remain eligible for another regen iteration
        cont_mask = ~exit_mask
        if np.any(cont_mask):
            cont_idx = surv[cont_mask]
            etauf[cont_idx] = ef2[cont_mask]
            dfinal[cont_idx] = dtau2[cont_mask]
            Pi[cont_idx] = Pint[cont_mask]
            ipp[cont_idx] = 0

        # Candidates that did not survive the neutrino step are dead as neutrinos
        dead_local = np.ones(cand.size, dtype=bool)
        dead_local[kept_local] = False
        if np.any(dead_local):
            active[cand[dead_local]] = False

            # --- consistency checks (debug) ---
    exit_from_energy = (e_out > 0.0)
    exit_from_flags  = (no_regen_tot == 1) | (regen_tot == 1)

    n = e_out.size
    #print("DEBUG run_stat:")
    #rint("  stats_cl:", n)
    #print("  nonzero e_out:", np.count_nonzero(exit_from_energy))
    #print("  exit flags:",   np.count_nonzero(exit_from_flags))
    #print("  mismatches:",   np.count_nonzero(exit_from_energy != exit_from_flags))

    # If mismatches are nonzero, this tells you which events disagree
    # (only print if small enough to not spam)
    mm = np.flatnonzero(exit_from_energy != exit_from_flags)
    if mm.size and mm.size < 20:
        print("  mismatch idx:", mm)
        print("  e_out[mm]:", e_out[mm])
        print("  no_regen_tot[mm]:", no_regen_tot[mm])
        print("  regen_tot[mm]:", regen_tot[mm])
        
    return no_regen_tot, regen_tot, e_out, P_out