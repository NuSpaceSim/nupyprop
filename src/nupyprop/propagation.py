"""
Created on Wed June 19 15:21:03 2024

@author: Luke Kupari and Diksha Garg
Comments: This code is the main backbone, doing the propagation of neutrinos, 
charged leptons, and regeneration of taus.

"""

import numpy as np
import nupyprop.transport as transport
import nupyprop.geometry as geometry
import nupyprop.constants as const

rho_water = const.rho_water # density of water, in g/cm3
step_size = const.step_size # step size for continuous energy loss, in cm
m_tau, m_mu = const.m_tau, const.m_mu # mass of tau and muon, in GeV
ctau_tau, ctau_mu = const.ctau_tau, const.ctau_mu # ctau of tau and muon, in cm

def propagate_nu(e_init, nu_xc, nu_ixc, depth_max, fac_nu, stats, Emin, E_nu, E_lep, yvals, return_index: bool = False):
    '''
    Propagates a neutrino inside the Earth.

    Parameters
    ----------
    e_init : float
        Inital neutrino energy, in GeV
    nu_xc : np.ndarray
        2D array containing neutrino CC & NC cross-section values, in cm^2
    nu_ixc : np.ndarray
        3D array containing neutrino integrated cross-section CDF values
    depth_max : float
        Maximum column depth for neutrino propagation, in kmwe
    fac_nu : float
        Rescaling factor for SM neutrino cross-sections
    stats : int
        Number of iterations to evaluate
    Emin : float
        Minimum threshold energy for leptons, in GeV
    E_nu : np.ndarray
        Array of neutrino energiesm in GeV
    E_lep : np.ndarray
        Array of charged lepton energies, in GeV
    yvals : np.ndarray
        Array of min. y values from which the cross-section CDF is calculated

    Returns
    -------
    part_type : integer array
        Type of outgoing particle. 0 = neutrino; 1 = charged lepton
    d_travel : float array
        Distance traveled until converted to charged lepton or total distance traveled by neutrino 
        (if no conversion to charged lepton), in kmwe
    e_fin : float array
        Final neutrino energy, in GeV
    '''
    part_type = np.zeros(stats, dtype=int)  # 0 = neutrino, 1 = charged lepton
    e_fin = e_init.astype(float).copy()  # Initialize all with e_init
    x_0 = np.zeros(stats, dtype=np.float32)  # Depth in kmwe
    d_travel = depth_max.astype(float).copy()  # Default to depth_max in kmwe

    active = np.ones(stats, dtype=bool)  # Track active simulations

    while np.any(active):  # Continue until all neutrinos stop
        r = np.random.random(stats)
        int_depth = transport.int_depth_nu(e_fin, nu_xc, fac_nu, E_nu)
        step_size = -int_depth * np.log(r) * 1e-5 # in kmwe

        x_0[active] += step_size[active]

        # Check which simulations exceeded total column depth
        exceeded = (x_0 > depth_max) & active
        active[exceeded] = False  # Stop these simulations

        # Compute interaction type
        int_type = transport.interaction_type_nu(e_fin, nu_xc, fac_nu, E_nu)
        converted = (part_type == 0) & (int_type == 0) & active # Neutrino to Charged Lepton
        part_type[converted] = 1  # Convert neutrinos to charged leptons

        # Compute energy loss
        y_fraction = transport.find_y(e_fin, nu_ixc, int_type, E_nu, E_lep, yvals)
        e_fin[active] *= (1 - y_fraction[active])

        # Check for charged leptons
        new_leptons = (part_type == 1) & active
        d_travel[new_leptons] = x_0[new_leptons]  # Update travel distances
        active[new_leptons] = False  # Stop these simulations

        # Check which events should stop due to energy loss
        energy_depleted = (e_fin <= Emin) & active
        d_travel[energy_depleted] = x_0[energy_depleted]
        active[energy_depleted] = False  # Stop simulations

    final_mask = ~((e_fin <= Emin) | (part_type == 0))  # eliminates low-energy charged leptons and neutrino-ended events

    if return_index:
        kept_local = np.nonzero(final_mask)[0]
        return kept_local, part_type[final_mask], d_travel[final_mask], e_fin[final_mask]

    return part_type[final_mask], d_travel[final_mask], e_fin[final_mask]
    
def propagate_lep(e_init, angle, xc, lep_ixc, alpha, beta, d_entry, d_in, lepton, prop_type,
    medium, # either 'water' or 'rock'
    cthi, Pi,
    idepth, earth_model, Emin, E_nu, E_lep, yvals, ypol, Pcthp, P,
    xalong = None, cdalong = None): 
    '''
    Propagates a charged lepton through Earth 

    Parameters
    ----------    
    e_init : float 
        Initial energy of the charged lepton, in GeV
    angle : float 
        Earth emergence angle (beta)
    xc : np.ndarray
        2D array containing N_A/A*charged lepton-nucleon cross-section values in water/rock, in cm^2/g
    lep_ixc : np.ndarray
        3D array containing charged lepton integrated cross-section CDF values in water/rock
    alpha : np.ndarray
        1D array containing ionization energy loss values in water/rock, in (GeV*cm^2)/g
    beta : np.ndarray 
        2D array of beta values in water/rock, in cm^2/g
    d_entry : float
        Column depth along the chord for a given Earth emergence angle, in kmwe. Defaults to None  
    d_in : float
        Maximum distance for charged lepton to propagate in water/rock, in kmwe
    lepton : integer
        Type of charged lepton. 1=tau; 2=muon
    prop_type : integer 
        Type of energy loss propagation. 1=stochastic, 2=continuous
    medium : string 
        charged lepton propagation medum, can be rock or water
    cthi : float 
        costheta value obtained from tau EM interaction with rock
    Pi : float 
        Degree of Polarization obtained from tau EM interaction with rock
    idepth : integer 
        Depth of water layer in km.
    earth_model : string
        prem or ak135 Earth model 
    Emin : float
        Minimum threshold energy for leptons, in GeV
    E_nu : np.ndarray
        Array of neutrino energiesm in GeV
    E_lep : np.ndarray
        Array of charged lepton energies, in GeV
    yvals : np.ndarray
        Array of min. y values from which the cross-section CDF is calculated
    ypol : float array
        Predefined inelasticity array
    Pcthp :
        np.cos(theta_P), where theta_P is the polar angle of the spin vector in tau rest frame
    P : float array
        Magnitude of polarization vector, defining the degree of polarization
    xalong : float 
        1D array containing distance, in km. Defaults to None
    cdalong : float 
        1D array containing column depth at xalong, in g/cm^2. Defaults to None

    Returns
    -------
    part_id : integer
        Type of outgoing charged lepton. 0=decayed; 1=not decayed; 2=don't count
    d_fin : float
        Distance traveled before charged lepton decays or total distance traveled by charged lepton, in kmwe
    e_fin : float
        Final energy of the charged lepton, in GeV
    cthf : float
        Final costheta value for polarization after EM interaction of the tau lepton
    Pf : float 
        Final degree of polarization value after EM interaction of the tau lepton
    '''
    rng = np.random.default_rng() #seed for random number generation

    stats = e_init.size 

    # Output arrays
    part_id = np.ones(stats, dtype=int) # default 1 = not decayed
    d_fin = np.zeros(stats, dtype=float) # final distance (kmwe)
    e_fin = e_init.astype(float).copy() # final energies (GeV)
    cthf = np.ones(stats, dtype=float)
    Pf = np.ones(stats, dtype=float)      

    energy = e_init.astype(float).copy()
    x0 = np.zeros(stats, dtype=float) # distance along chord in cmwe
    finished = np.zeros(stats, dtype=bool) # finished particles

    # masses and ctau
    m_le = np.where(lepton == 1, m_tau, m_mu)
    ctau_le = np.where(lepton == 1, ctau_tau, ctau_mu)

    x_total = d_in * 1e5 # in cmwe
    col_depth = d_entry * 1e5  # only meaningful for rock initially, in cmwe

    # initial Pin and theta
    Pin = Pi.astype(float).copy()
    theta_in = np.arccos(cthi)  

    # -------------------------
    # STOCHASTIC PROPAGATION
    # -------------------------
    if prop_type == 1:
        alive = (energy > Emin)  # main mask to check for alive particles
        while np.any(alive):
            idx = np.nonzero(alive)[0] # indices of alive particles
            if idx.size == 0:
                break

            # fetch per-particle arrays
            E_alive = energy[idx]

            # IMPORTANT: use physical mass/ctau arrays computed above, not the lepton ID
            m_le_alive = m_le[idx]
            ctau_le_alive = ctau_le[idx]

            x0_alive = x0[idx]
            x_total_alive = x_total[idx] # cmwe
            
            if medium == 'water':
                rho_a = np.full(idx.size, rho_water, dtype=float)

            else: # for rock
                # current column depth for these particles:
                col_depth_alive = col_depth[idx] # cmwe
                x_interp_rock = transport.cd2distd(xalong, cdalong, col_depth_alive)
                _, rho_a = geometry.densityatx(x_interp_rock, angle, idepth, earth_model)
                                                  
            # --- interaction length and random step ---
            int_depth = transport.int_depth_lep(E_alive, xc, rho_a, m_le_alive, ctau_le_alive, E_lep)  # returns length in cmwe
            # sample exponential distances vectorized
            u = rng.random(int_depth.size)
            xs = -int_depth * np.log(u)  # cmwe

            x_f = x0_alive + xs
            escaped_mask = (x_f >= x_total_alive)   # those that reach/exit before interacting

            # --- handles particles that reach x_total this step ---
            if np.any(escaped_mask):
                esc_local = np.nonzero(escaped_mask)[0]
                esc_global = idx[esc_local]

                # remaining distance to boundary (cmwe)
                x_step = x_total_alive[esc_local] - x0_alive[esc_local]

                # continuous loss applied over x_step
                E_before = energy[esc_global]
                a_vals = transport.int_alpha(E_before, alpha, E_lep)
                b_vals = transport.int_beta(E_before, beta, rho_a[esc_local], E_lep)
                e_after = E_before - (E_before * b_vals + a_vals) * x_step
                e_after = np.maximum(e_after, Emin)

                # store outputs
                e_fin[esc_global] = e_after
                d_fin[esc_global] = d_in[esc_global]   # they reached the provided d_in (kmwe)
                part_id[esc_global] = np.where(e_after <= Emin, 2, 1) 
                cthf[esc_global] = np.cos(theta_in[esc_global])
                Pf[esc_global] = Pin[esc_global]
                finished[esc_global] = True
                alive[esc_global] = False

            # --- handle interacting particles (those that interact before reaching boundary) ---
            interact_local = np.nonzero(~escaped_mask)[0]
            if interact_local.size > 0:
                int_global = idx[interact_local]

                # advance x0 by xs for these particles (cmwe)
                x0[int_global] += xs[interact_local]

                # For rock particles, update their col_depth by xs (they moved through matter)
                # note: col_depth is per-global-particle
                if medium == 'rock':
                    col_depth[int_global] += xs[interact_local]

                # continuous losses across the step xs
                E_before = energy[int_global]
                # IMPORTANT: use physical mass/ctau arrays computed above, not the lepton ID
                m_le_step = m_le[int_global]
                ctau_le_step = ctau_le[int_global]
                a_vals = transport.int_alpha(E_before, alpha, E_lep)
                b_vals = transport.int_beta(E_before, beta, rho_a[interact_local], E_lep)
                e_int = E_before - (E_before * b_vals + a_vals) * xs[interact_local]
                e_int = np.maximum(e_int, Emin)

                e_avg = 10 ** ((np.log10(E_before) + np.log10(e_int)) / 2.0)
                a_avg = transport.int_alpha(e_avg, alpha, E_lep)
                b_avg = transport.int_beta(e_avg, beta, rho_a[interact_local], E_lep)
                e_int = transport.em_cont_part(E_before, a_avg, b_avg, xs[interact_local], m_le_step)

                # update energies in state
                energy[int_global] = e_int
                e_fin[int_global] = e_int

                # Those that fell to or below Emin after the step -> mark as "don't count" 
                low_mask = (e_int <= Emin)
                if np.any(low_mask):
                    low_global = int_global[low_mask]
                    part_id[low_global] = 2 
                    d_fin[low_global] = d_in[low_global]
                    e_fin[low_global] = Emin
                    cthf[low_global] = np.cos(theta_in[low_global])
                    Pf[low_global] = Pin[low_global]
                    finished[low_global] = True
                    alive[low_global] = False

                # remaining indices that still have e_int > e_min:
                cont_mask = ~low_mask
                if np.any(cont_mask):
                    cont_local = np.nonzero(cont_mask)[0]
                    cont_global = int_global[cont_local]

                    # compute interaction types vectorized
                    e_for_int = e_int[cont_local]
                    rho_for_int = rho_a[interact_local][cont_local]

                    m_le_for_int = m_le_step[cont_local]
                    ctau_le_for_int = ctau_le_step[cont_local]
                    int_types = transport.interaction_type_lep(e_for_int, xc, rho_for_int, m_le_for_int, ctau_le_for_int, E_lep)

                    # decays: int_type == 2
                    dec_mask = (int_types == 2)
                    if np.any(dec_mask):
                        dec_local = np.nonzero(dec_mask)[0]
                        dec_global = cont_global[dec_local]
                        part_id[dec_global] = 0
                        d_fin[dec_global] = x0[dec_global] / 1e5   # kmwe
                        e_fin[dec_global] = energy[dec_global] 
                        cthf[dec_global] = np.cos(theta_in[dec_global])
                        Pf[dec_global] = Pin[dec_global]
                        finished[dec_global] = True
                        alive[dec_global] = False

                    # for interactions that do not decay, sample 'y' and update energies
                    survive_mask = ~dec_mask
                    if np.any(survive_mask):
                        surv_local = np.nonzero(survive_mask)[0]
                        surv_global = cont_global[surv_local]

                        e_after_int = e_for_int[surv_local]
                        int_types_surv = int_types[survive_mask]

                        # sample y (vectorized)
                        y = transport.find_y(e_after_int, lep_ixc, int_types_surv, E_nu, E_lep, yvals)

                        E_new = e_after_int * (1.0 - y)
                        E_new = np.maximum(E_new, Emin)

                        # update energies
                        energy[surv_global] = E_new
                        e_fin[surv_global] = E_new

                        # polarization updates when int_type == 5
                        pol_mask = (int_types_surv == 5)
                        if np.any(pol_mask):
                            pol_local = np.nonzero(pol_mask)[0]
                            pol_global = surv_global[pol_local]
                            y_pol = y[pol_local]
                            Pout_vals, theta_out_vals = transport.polarization(y_pol, Pin[pol_global], theta_in[pol_global], ypol, Pcthp, P)
                            Pin[pol_global] = Pout_vals
                            theta_in[pol_global] = theta_out_vals

                        # those that drop to <= e_min after inelastic update => mark finished
                        died_mask = (E_new <= Emin)
                        if np.any(died_mask):
                            died_global = surv_global[died_mask]
                            part_id[died_global] = 2 
                            d_fin[died_global] = d_in[died_global]
                            e_fin[died_global] = Emin
                            cthf[died_global] = np.cos(theta_in[died_global])
                            Pf[died_global] = Pin[died_global]
                            finished[died_global] = True
                            alive[died_global] = False

            alive = (~finished) & (energy > Emin) # updating alive mask

        # --- after main loop: set outputs for any not finished (they didn't interact or decayed)
        not_finished = ~finished
        if np.any(not_finished):
            # final distance is x_total (cmwe) converted to kmwe
            part_id[not_finished] = 1 # still leptons
            d_fin[not_finished] = x_total[not_finished] / 1e5
            e_fin[not_finished] = np.maximum(energy[not_finished], Emin)
            cthf[not_finished] = np.cos(theta_in[not_finished])
            Pf[not_finished] = Pin[not_finished]

        return part_id, d_fin, e_fin, cthf, Pf
    
    # -------------------------
    # CONTINUOUS PROPAGATION
    # -------------------------
    else: 
        delta_x = step_size  # cm
        d_max = d_in * 1e5  # cmwe

        # --- initialize densities for rock ---
        if medium == "rock":
            x_interp = transport.cd2distd(xalong, cdalong, col_depth)
            _, rho_rock_vals = geometry.densityatx(x_interp, angle, idepth, earth_model)
            rho = rho_rock_vals
            j_max = (d_max / (rho * delta_x)).astype(int)
        else:
            rho = np.full(col_depth.size, rho_water, dtype=float)
            j_max = (x_total / (delta_x * rho)).astype(int)

        max_steps = np.max(j_max)

        # --- main propagation loop --- 
        for _ in range(max_steps): #---> is this loop conditions correct? 
            active = (~finished) & (energy > Emin)

            if not np.any(active):
                break

            # update rock density if needed
            if medium == 'rock':
                x_interp = transport.cd2distd(xalong, cdalong, col_depth[active])
                _, rho_rock_vals = geometry.densityatx(x_interp, angle, idepth, earth_model)
                rho[active] = rho_rock_vals

            # depth step
            delta_d = delta_x * rho[active]
            if medium == 'rock':
                col_depth[active] = d_entry[active] * 1e5 + x0[active]

            x0[active] += delta_d
        
            # decay check
            part_id_step = transport.idecay(energy[active], delta_x, m_le[active], ctau_le[active])
            decayed_local = (part_id_step == 0)

            if np.any(decayed_local):
                active_idx = np.where(active)[0]
                decayed = active_idx[decayed_local]   # map back to global indices
                finished[decayed] = True
                part_id[decayed] = 0
                e_fin[decayed] = energy[decayed]
                d_fin[decayed] = x0[decayed] / 1e5

            # energy loss
            alpha_val = transport.int_alpha(energy[active], alpha, E_lep)
            beta_val = transport.int_beta(energy[active], beta, rho[active], E_lep)
            energy[active] -= (energy[active] * beta_val + alpha_val) * delta_d

            # update final values for survivors
            e_fin[active] = energy[active]
            d_fin[active] = x0[active] / 1e5

            # stop conditions ---> FROM HERE WORK ON IT
            reached_max_local = (x0[active] >= x_total[active])
            below_local = (energy[active] <= Emin)

            active_idx = np.where(active)[0]  # map subset → global
            stop_global = active_idx[reached_max_local | below_local]
            below_global = active_idx[below_local]

            if stop_global.size > 0:
                finished[stop_global] = True
                d_fin[stop_global] = x0[stop_global] / 1e5

            if below_global.size > 0:
                part_id[below_global] = 2
                e_fin[below_global] = Emin

        # --- finalize all particles ---
        not_finished = ~finished
        if np.any(not_finished):
            d_fin[not_finished] = np.where(
                medium == 'rock', d_max[not_finished] / 1e5, x_total[not_finished] / 1e5
            )
            e_fin[not_finished] = np.maximum(energy[not_finished], Emin)
            Pf[not_finished] = Pin[not_finished]
            cthf[not_finished] = np.cos(theta_in[not_finished])

        # --- outputs ---
        return part_id, d_fin, e_fin, cthf, Pf

def tau_thru_layers(
    angle, depth, d_water, depth_traj, e_lep_in, xc_water, xc_rock, lep_ixc_water,
    lep_ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong,
    idepth, lepton, prop_type, Emin, E_nu, E_lep, yvals, ypol, Pcthp, P, earth_model, stats
):
    """
    Paramters
    ---------
        angle : float
            Earth emergence angle (beta)
        depth : float
            Total column depth of the chord, in kmwe.
        d_water : float
            Column depth of the final layer of water (or full distance in water if only water layer), in kmwe.
        depth_traj : float
            Column depth along the chord for a given Earth emergence angle, in kmwe.
        e_lep_in : float
            Ingoing charged lepton energy, in GeV.
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
        idepth : Integer 
            Depth of water layer in km.
        lepton : 
            Integer Type of charged lepton. 1=tau; 2=muon.
        prop_type : 
            Integer Type of energy loss propagation. 1=stochastic, 2=continuous.
        Emin : float
            Minimum threshold energy for leptons, in GeV
        E_nu : np.ndarray
            Array of neutrino energies in GeV
        E_lep : np.ndarray
            Array of charged lepton energies, in GeV
        yvals : np.ndarray
            Array of min. y values from which the cross-section CDF is calculated
        ypol : float array
            Predefined inelasticity array
        Pcthp :
            np.cos(theta_P), where theta_P is the polar angle of the spin vector in tau rest frame
        P : float array
            Magnitude of polarization vector, defining the degree of polarization
        earth_model : str
            Earth density model, prem or ak135.

    Returns
    -------
        part_type : np.ndarray
            Type of outgoing charged lepton. 0=decayed, 1=not decayed.
        d_fin : np.ndarray
            Distance traveled before charged lepton decays or total distance traveled by charged lepton, in kmwe.
        e_fin : np.ndarray
            Outgoing charged lepton energy, in GeV.
        pcthf : np.ndarray
            Final polarization after EM interaction of the tau lepton
    """
    # Normalize inputs to arrays
    depth = np.asarray(depth, dtype=float)
    d_water = np.asarray(d_water, dtype=float)
    depth_traj = np.asarray(depth_traj, dtype=float)
    stats = e_lep_in.size

    lepton = np.asarray(lepton, dtype=int)
    # Broadcast scalar lepton id (e.g. 1 or 2) to per-event array
    if lepton.ndim == 0 or lepton.size == 1:
        lepton = np.full(stats, int(lepton.item()), dtype=int)
    elif lepton.size != stats:
        raise ValueError(f"lepton must be scalar or length {stats}, got length {lepton.size}")

    part_type = np.ones(stats, dtype=int)
    d_fin = depth_traj.copy()
    e_fin = e_lep_in.copy()
    pcthf = np.ones(stats, dtype=float)

    # Low-energy exit
    low_energy = e_lep_in < 1e3
    if np.any(low_energy):
        part_type[low_energy] = 0
        d_fin[low_energy] = depth[low_energy]
        e_fin[low_energy] = e_lep_in[low_energy]
    active = ~low_energy

    if not np.any(active):
        return part_type, d_fin, e_fin, pcthf
    
    ones = np.ones_like(e_lep_in)
    
    remaining = depth - depth_traj
    
    #change to not be built off stats
    rho_now = np.full(stats, rho_water, dtype=float)

    need_rho = active & (remaining >= d_water)
    if np.any(need_rho):
        col_depths = depth_traj[need_rho] * 1e5
        x_interp = transport.cd2distd(xalong, cdalong, col_depths)
        r, rho_vals = geometry.densityatx(x_interp, angle, idepth, earth_model)

        rho_now[need_rho] = rho_vals
            
    in_rock_now = active & (rho_now > 1.5) 
    in_water_now = active & ~in_rock_now

    # --- Rock -> Water branch ---
    if np.any(in_rock_now):
        rock_idx = np.where(in_rock_now)[0]

        d_in_rock = depth[in_rock_now] - depth_traj[in_rock_now] - d_water[in_rock_now]
        d_in_rock = np.maximum(d_in_rock, 0.0)

        r_part, r_d, r_e, r_cth, r_P = propagate_lep(
            e_lep_in[in_rock_now], angle, xc_rock, lep_ixc_rock,
            alpha_rock, beta_rock, depth_traj[in_rock_now], d_in_rock,
            lepton[in_rock_now], prop_type, 'rock', ones[in_rock_now], ones[in_rock_now],
            idepth, earth_model, Emin, E_nu, E_lep, yvals, ypol, Pcthp, P,
            xalong=xalong, cdalong=cdalong
        )
        
        depth_after_rock = depth_traj[in_rock_now] + r_d

        # Survivors propagate through water if a water layer exists
        survivors = (r_part == 1) & (idepth != 0)
        if np.any(survivors):
            surv_idx = rock_idx[survivors]
            w_part, w_d, w_e, w_cth, w_P = propagate_lep(
                r_e[survivors], angle, xc_water, lep_ixc_water,
                alpha_water, beta_water, depth_after_rock[survivors], d_water[in_rock_now][survivors],
                lepton[surv_idx], prop_type, 'water', r_cth[survivors], r_P[survivors],
                idepth, earth_model, Emin, E_nu, E_lep, yvals, ypol, Pcthp, P
            )
            part_type[surv_idx] = w_part
            d_fin[surv_idx] = depth_after_rock[survivors] + w_d
            e_fin[surv_idx] = w_e
            pcthf[surv_idx] = w_P * w_cth
        # stop/decay in rock or no water layer
        non_surv = ~survivors
        if np.any(non_surv):
            ns_idx = rock_idx[non_surv]
            part_type[ns_idx] = r_part[non_surv]
            d_fin[ns_idx] = depth_after_rock[non_surv]
            e_fin[ns_idx] = r_e[non_surv]
            pcthf[ns_idx] = r_P[non_surv] * r_cth[non_surv]

        active[in_rock_now] = False

    # --- Water-only branch (full remaining distance) ---
    if np.any(in_water_now):
        w_part, w_d, w_e, w_cth, w_P = propagate_lep(
            e_lep_in[in_water_now], angle, xc_water, lep_ixc_water,
            alpha_water, beta_water, depth_traj[in_water_now], remaining[in_water_now],
            lepton[in_water_now], prop_type, 'water', ones[in_water_now], ones[in_water_now],
            idepth, earth_model, Emin, E_nu, E_lep, yvals, ypol, Pcthp, P,
            xalong=xalong, cdalong=cdalong
        )
        part_type[in_water_now] = w_part
        d_fin[in_water_now] = depth_traj[in_water_now] + w_d
        e_fin[in_water_now] = w_e
        pcthf[in_water_now] = w_P * w_cth

        active[in_water_now] = False

    return part_type, d_fin, e_fin, pcthf

def distnu_arr(r, ithird, Pin, tol=1e-3, max_iter=25):
    '''
    Determines the neutrino energy from tau decay. The energy fraction is determined by tau energy CDF
    approximated by leptonic decay channel. Approximate (good enough) for taus; exact for muons.
    
    Vectorized sampler for the distnu distribution using batched bisection.

    Parameters
    ----------
    r : array-like
        Random numbers in [0, 1).
    ithird : array-like or int
        Selector for neutrino-to-lepton energy fraction (1 -> fixed 1/3).
    Pin : array-like
        Tau polarization values.
    tol : float, optional
        Bisection stopping tolerance on interval width.
    max_iter : int, optional
        Maximum number of bisection iterations.

    Returns
    -------
    np.ndarray
        Sampled energy fractions y, y = E_nu_tau/E_tau
    '''
    r = np.asarray(r, dtype=float)
    ithird = np.asarray(ithird)
    Pin = np.asarray(Pin, dtype=float)

    out = np.empty_like(r, dtype=float)

    # ithird == 1 -> fixed 1/3
    mask_third = ithird == 1
    out[mask_third] = 1.0 / 3.0

    mask_calc = ~mask_third
    if np.any(mask_calc):
        rr = r[mask_calc]
        pin = Pin[mask_calc]

        y_lo = np.zeros_like(rr)
        y_hi = np.ones_like(rr)

        for _ in range(max_iter):
            mid = 0.5 * (y_lo + y_hi)
            term1 = (mid / 3.0) * (5.0 - 3.0 * mid**2 + mid**3)
            term2 = (-1.0 * pin) * (mid / 3.0 * (1.0 - 3.0 * mid**2 + 2.0 * mid**3))
            fm = term1 + term2

            go_high = fm < rr
            y_lo = np.where(go_high, mid, y_lo)
            y_hi = np.where(go_high, y_hi, mid)

            if np.max(y_hi - y_lo) < tol:
                break

        out[mask_calc] = 0.5 * (y_lo + y_hi)

    return out

def regen_compact(angle, e_lep, depth, d_water, d_lep, nu_xc, nu_ixc, ithird, xc_water, xc_rock,
                  ixc_water, ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, idepth,
                  lepton, fac_nu, prop_type, Pin, Emin, E_nu, E_lep, yvals, ypol, Pcthp, P, earth_model):
    '''
    Regeneration loop, should take a pin and also throw out pout. pin will be used by distnu and the pout
       it throws will be used by the regen again as input in the single_stat() function.

    This mirrors the scalar/Fortran logic but allows compaction by returning the mapping from inputs
    (the arrays passed to this function) to the compacted survivor arrays.

    Parameters
    ----------
    angle : float
        Earth emergence angle (beta), in degrees.
    e_lep : np.ndarray
        Incoming charged lepton energies, in GeV.
    depth : np.ndarray
        Total column depths of the chord, in kmwe.
    d_water : np.ndarray
        Column depths of the terminal water layer, in kmwe.
    d_lep : np.ndarray
        Column depths along the chord at entry, in kmwe.
    nu_xc : np.ndarray
        Neutrino CC & NC cross-sections, in cm^2.
    nu_ixc : np.ndarray
        Neutrino integrated cross-section CDF values.
    ithird : np.ndarray or int
        Selector for neutrino-to-lepton energy fraction.
    xc_water : np.ndarray
        Lepton-nucleon cross-sections in water, N_A/A * sigma, in cm^2/g.
    xc_rock : np.ndarray
        Lepton-nucleon cross-sections in rock, N_A/A * sigma, in cm^2/g.
    ixc_water : np.ndarray
        Lepton integrated cross-section CDF values in water.
    ixc_rock : np.ndarray
        Lepton integrated cross-section CDF values in rock.
    alpha_water : np.ndarray
        Ionization energy loss in water, in (GeV*cm^2)/g.
    alpha_rock : np.ndarray
        Ionization energy loss in rock, in (GeV*cm^2)/g.
    beta_water : np.ndarray
        Beta values in water, in cm^2/g.
    beta_rock : np.ndarray
        Beta values in rock, in cm^2/g.
    xalong : np.ndarray
        Distances in water, in km.
    cdalong : np.ndarray
        Column depths at xalong, in g/cm^2.
    idepth : int
        Depth of water layer, in km.
    lepton : np.ndarray or int
        1=tau; 2=muon.
    fac_nu : float
        Rescaling factor for SM neutrino cross-sections.
    prop_type : int
        1=stochastic; 2=continuous.
    Pin : np.ndarray
        Input polarization from prior tau propagation.
    Emin : float
        Minimum lepton energy threshold, in GeV.
    E_nu : np.ndarray
        Neutrino energy grid, in GeV.
    E_lep : np.ndarray
        Lepton energy grid, in GeV.
    yvals : np.ndarray
        Minimum y values for the cross-section CDF.
    ypol : np.ndarray
        Predefined inelasticity array.
    Pcthp : np.ndarray
        cos(theta_P) lookup for polarization.
    P : np.ndarray
        Polarization magnitude lookup.
    earth_model : str
        Earth density model, "prem" or "ak135".

    Returns
    -------
    ipp3 : integer array (k,)
        Type of outgoing particle. 0=neutrino; 3=exit.
    dtau2 : float array (k,)
        Distance traveled before charged lepton decays or total distance traveled by charged lepton, in kmwe.
    ef2 : float array (k,)
        Final particle energy, in GeV.
    Pout : float array (k,)
        Tau's polarization output from regen function
    '''
    e_lep = np.asarray(e_lep, dtype=float)
    n = e_lep.size

    def _as_float_array(x):
        arr = np.asarray(x, dtype=float)
        if arr.ndim == 0:
            arr = np.full(n, float(arr.item()), dtype=float)
        return arr

    depth = _as_float_array(depth)
    d_water = _as_float_array(d_water)
    d_lep = _as_float_array(d_lep)
    Pin = _as_float_array(Pin)

    lepton_arr = np.asarray(lepton, dtype=int)
    if lepton_arr.ndim == 0 or lepton_arr.size == 1:
        lepton_arr = np.full(n, int(lepton_arr.item()), dtype=int)
    elif lepton_arr.size != n:
        raise ValueError(f"lepton must be scalar or length {n}, got length {lepton_arr.size}")

    # 1) tau decay -> neutrino energy
    r = np.random.random(n)
    frac = distnu_arr(r, ithird, Pin)
    e_nu = frac * e_lep

    # 2) available depth for neutrino propagation
    d_left0 = depth - d_lep
    room0 = (d_left0 > 0.0)
    if not np.any(room0):
        return (np.zeros(0, dtype=int),
                np.zeros(0, dtype=int),
                np.zeros(0, dtype=float),
                np.zeros(0, dtype=float),
                np.zeros(0, dtype=float))

    idx0 = np.flatnonzero(room0)

    # 3) propagate neutrinos; this compacts to CC->tau survivors above Emin
    kept0, _ptype, dtr, etau2 = propagate_nu(
        e_nu[idx0], nu_xc, nu_ixc, d_left0[idx0], fac_nu,
        stats=idx0.size,
        Emin=Emin, E_nu=E_nu, E_lep=E_lep, yvals=yvals,
        return_index=True,
    )

    if kept0.size == 0:
        return (np.zeros(0, dtype=int),
                np.zeros(0, dtype=int),
                np.zeros(0, dtype=float),
                np.zeros(0, dtype=float),
                np.zeros(0, dtype=float))

    kept_local = idx0[kept0]

    # 4) tau entry point and remaining depth (Fortran: if no room, drop)
    d_lep1 = d_lep[kept_local] + dtr
    d_left1 = depth[kept_local] - d_lep1
    room1 = (d_left1 > 0.0)
    if not np.any(room1):
        return (np.zeros(0, dtype=int),
                np.zeros(0, dtype=int),
                np.zeros(0, dtype=float),
                np.zeros(0, dtype=float),
                np.zeros(0, dtype=float))

    kept_local = kept_local[room1]
    etau2 = etau2[room1]
    d_lep1 = d_lep1[room1]

    # 5) propagate tau through layers
    ipp3, dtau2, ef2, Pout = tau_thru_layers(
        angle,
        depth[kept_local],
        d_water[kept_local],
        d_lep1,
        etau2,
        xc_water,
        xc_rock,
        ixc_water,
        ixc_rock,
        alpha_water,
        alpha_rock,
        beta_water,
        beta_rock,
        xalong,
        cdalong,
        idepth,
        lepton_arr[kept_local],
        prop_type,
        Emin,
        E_nu,
        E_lep,
        yvals,
        ypol,
        Pcthp,
        P,
        earth_model,
        stats=kept_local.size,
    )

    return kept_local, ipp3, dtau2, ef2, Pout
