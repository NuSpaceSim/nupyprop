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

def propagate_nu(e_init, nu_xc, nu_ixc, depth_max, fac_nu, stats, Emin, E_nu, E_lep, yvals):
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

    final_mask = ~((e_fin <= Emin)) # | (part_type == 0)) # this will eliminate charged leptons 
                                                       # with energy < Emin and neutrinos

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
            int_depth = transport.int_depth_lep(E_alive, xc, rho_a, m_le, ctau_le, E_lep)  # returns length in cmwe
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
                a_vals = transport.int_alpha(E_before, alpha)
                b_vals = transport.int_beta(E_before, beta, rho_a[esc_local])
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
                a_vals = transport.int_alpha(E_before, alpha, E_lep)
                b_vals = transport.int_beta(E_before, beta, rho_a[interact_local], E_lep)
                e_int = E_before - (E_before * b_vals + a_vals) * xs[interact_local]
                e_int = np.maximum(e_int, Emin)

                e_avg = 10 ** ((np.log10(E_before) + np.log10(e_int)) / 2.0)
                a_avg = transport.int_alpha(e_avg, alpha, E_lep)
                b_avg = transport.int_beta(e_avg, beta, rho_a[interact_local], E_lep)
                e_int = transport.em_cont_part(E_before, a_avg, b_avg, xs[interact_local], m_le)

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

                    int_types = transport.interaction_type_lep(e_for_int, xc, rho_for_int, m_le, ctau_le, E_lep)

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
            part_id_step = transport.idecay(energy[active], delta_x, m_le, ctau_le)
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

def tau_thru_layers(angle, depth,d_water, depth_traj, e_lep_in, xc_water, xc_rock, lep_ixc_water,
                    lep_ixc_rock, alpha_water, alpha_rock, beta_water,
                    beta_rock, xalong, cdalong, idepth, lepton, 
                    prop_type, Emin, E_nu, E_lep, yvals, ypol, 
                    Pcthp, P, earth_model, stats):
    """
    Vectorized version - arrays of input parameters
    """
    n_events = stats
    
    # Initialize output arrays
    part_types = np.ones(n_events, dtype=int)
    d_fins = depth_traj.copy()
    e_fins = e_lep_in.copy()
    pcthfs = np.ones(n_events)
    
    # Vectorized energy threshold check
    low_energy_mask = (e_lep_in < Emin)
    part_types[low_energy_mask] = 0
    d_fins[low_energy_mask] = depth[low_energy_mask]
    
    active_mask = ~low_energy_mask
    col_depths = depth_traj[active_mask]*1e5    
    
    if not np.any(active_mask):
        return ( part_types[low_energy_mask],
                 d_fins[low_energy_mask],
                 e_fins[low_energy_mask],
                 pcthfs[low_energy_mask])
    
    active_indices = np.where(active_mask)[0] # get indices for all events we want to propagate

    water_condition = (depth[active_mask] - depth_traj[active_mask]) < d_water[active_mask]
    
    rhos = np.full(np.sum(active_mask), rho_water)
    
    rock_mask = ~water_condition
    if np.any(rock_mask):
        #rock_indices = active_indices[rock_mask]
        active_rock_indices = np.where(rock_mask)[0] # get active subset
        for j in active_rock_indices:
            orig_index = active_indices[j]
            x = transport.cd2distd(xalong, cdalong, col_depths[j])
            _, rho = geometry.densityatx(x, angle, idepth, model_name='prem')
            rhos[j] = rho
            
        # split into processing groups
        in_rock_mask_active = (rhos > 1.5)
        in_water_mask_active = ~in_rock_mask_active
        
        if np.any(in_rock_mask_active):
            rock_active_subix = np.where(in_rock_mask_active)[0]
            rock_indices = active_indices[rock_active_subix]
            
            for local_index, orig_index in enumerate(rock_indices):
                e_lep = e_fins[orig_index]
                depth_traj = depth_traj[orig_index]
                d_in_rock = depth[orig_index] - depth_traj - d_water[orig_index]
                
                #propagate through rock
                part_type_r, df_r, efin_r, cthf, Pf = propagate_lep_rock(
                    angle, e_lep, xc_rock, lep_ixc_rock, alpha_rock, beta_rock,
                    depth_traj, d_in_rock, xalong, cdalong, idepth, lepton,
                    prop_type, Emin, E_nu, E_lep, yvals, ypol, Pcthp, P, earth_model
                )
                
                depth_traj_after_rock = depth_traj + df_r
                
                if part_type_r == 1 and idepth != 0:
                    #water propagation after rock
                    e_lep_water = efin_r
                    d_in_water = d_water[orig_index]
                    Pi = Pf
                    cthi = cthf
                    
                    part_type_w, df_w, e_fin_w, ptchf_w = propagate_lep_water(
                        e_lep_water, xc_water, lep_ixc_water, alpha_water,
                        beta_water, d_in_water, lepton, prop_type,
                        cthi, Pi, Emin, E_nu, E_lep, yvals, ypol, Pcthp, P
                    )
    return part_types, d_fins, e_fins, pcthfs

def distnu(r, ithird, Pin):
    '''
    Determines the neutrino energy from tau decay. The energy fraction is determined by tau energy CDF
    approximated by leptonic decay channel. Approximate (good enough) for taus; exact for muons

    Args:
        r (Float): Random number
        ithird (IntegerL): Choice for neutrino -> charged lepton energy fraction selection
        Pin (Float): Tau's polarization

    Returns:
        dist_val (Float): Energy fraction , y = E_nu_tau/E_tau
    '''
    #P = -1.0 * Pin, -1 = fully LH polarized; 0 = fully unpolarized tau
    def fnu(y):

        return y/3.0 * (5.0 - 3.0*y**2 + y**3) + (-1.0 *Pin) * (y/3.0 *(1.0 -3.0 * y**2 + 2.0 *y**3))

    if ithird !=1:
        fm = 1.0 #max value of distribution
        ff = r*fm
        y0 = 0.0
        y1 = fm

        while abs(y1-y0) > 0.1e-2:
            y = (y0+y1)/2.0
            if fnu(y) < ff:
                y0 = y
            else:
                y1 = y
            dist_val = (y0+y1)/2.0

    else: # if ithird ==1; use 1/3 of energy of 3 body decay.
        dist_val = 1.0/3.0

    return dist_val

def regen(angle, e_lep, depth, d_water, d_lep, nu_xc, nu_ixc, ithird, xc_water, xc_rock,
            ixc_water, ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, idepth,
            lepton, fac_nu, prop_type, Pin, Emin, E_nu, E_lep, yvals, ypol, Pcthp, P, earth_model, stats):


    n = int(stats)

    # Outputs
    part_type = np.full(n, 3, dtype=int)     # default 'exit' sentinel like original
    d_exit    = d_lep.copy()
    e_fin     = np.empty(n, dtype=float)
    Pout      = Pin.copy()

    # 1) Sample y = E_nu/E_tau from distnu
    r = np.random.random(n)
    distnu_vec = np.vectorize(lambda rr, it, pin: distnu(float(rr), int(it), float(pin)), otypes=[float])
    frac = distnu_vec(r, ithird, Pin)
    e_nu = frac * e_lep
    e_fin[:] = e_nu

    # 2) Geometry mask – remaining depth the neutrino can travel
    d_left = depth - d_lep
    no_room = (d_left <= 0.0)
    if np.any(no_room):
        d_exit[no_room]    = depth[no_room]
        part_type[no_room] = 3
        Pout[no_room]      = Pin[no_room]

    # 3) Propagate neutrinos for the rest
    nu_active = ~no_room
    if np.any(nu_active):
        act_idx = np.flatnonzero(nu_active)

        # propagate_nu expects a batch; assume it returns arrays aligned to inputs
        ip, dtr, etau2 = propagate_nu(
            e_nu[nu_active], nu_xc, nu_ixc, d_left[nu_active], fac_nu,
            np.count_nonzero(nu_active), Emin, E_nu, E_lep, yvals
        )

        # ip != 1 -> still a neutrino at the end
        nu_ends_mask = (ip != 1)
        if np.any(nu_ends_mask):
            end_idx = act_idx[nu_ends_mask]
            part_type[end_idx] = 0
            d_exit[end_idx]    = depth[end_idx]
            e_fin[end_idx]     = etau2[nu_ends_mask]
            Pout[end_idx]      = -2.0  # sentinel used by upstream code

        # The rest produced taus
        tau_mask_global = nu_active.copy()
        if np.any(nu_ends_mask):
            tau_mask_global[act_idx[nu_ends_mask]] = False

        if np.any(tau_mask_global):
            tau_idx = np.flatnonzero(tau_mask_global)

            # Set up per-event tau inputs
            d_lep_tau  = d_lep[tau_idx] + dtr[~nu_ends_mask]
            d_left_tau = d_left[tau_idx] - dtr[~nu_ends_mask]
            etau2_tau  = etau2[~nu_ends_mask]

            # If no room after CC, finalize as neutrino-like outcome
            tau_no_room = (d_left_tau <= 0.0)
            if np.any(tau_no_room):
                nr_idx = tau_idx[tau_no_room]
                part_type[nr_idx] = 0
                d_exit[nr_idx]    = depth[nr_idx]
                e_fin[nr_idx]     = etau2_tau[tau_no_room]
                Pout[nr_idx]      = Pin[nr_idx]

            # Go for tau through layers where there is room
            tau_go_mask = ~tau_no_room
            if np.any(tau_go_mask):
                go_idx = tau_idx[tau_go_mask]

                ptype_sub, dexit_sub, efin_tau_sub, Pi_sub = tau_thru_layers(
                    angle[go_idx], depth[go_idx], d_water[go_idx], d_lep_tau[tau_go_mask],
                    etau2_tau[tau_go_mask],
                    xc_water, xc_rock, ixc_water, ixc_rock, alpha_water, alpha_rock,
                    beta_water, beta_rock, xalong, cdalong, idepth, lepton, prop_type,
                    Emin, E_nu, E_lep, yvals, ypol, Pcthp, P, earth_model,
                    stats=len(go_idx)
                )

                part_type[go_idx] = ptype_sub
                d_exit[go_idx]    = dexit_sub
                e_fin[go_idx]     = efin_tau_sub
                Pout[go_idx]      = Pi_sub

    return part_type, d_exit, e_fin, Pout