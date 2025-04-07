"""
Created on Wed June 19 15:21:03 2024

@author: Luke Kupari

"""

import numpy as np
import nupyprop.constants as const
import nupyprop.transport as transport
import nupyprop.geometry as geometry

#rho_rock = const.rho_rock # rock density
#rho_iron = const.rho_iron # iron density

#E_nu = const.E_nu # Neutrino energy numpy array
#E_lep = const.E_lep # Lepton energy numpy array
#yvals = const.yvals # inelasticity from 1e-3 to 1


def propagate_nu(e_init, nu_xc, nu_ixc, depth_max, fac_nu):

    """Propagates a neutrino inside the Earth.

    Args:
        e_init (float): Inital neutrino energy, in GeV
        nu_xc (np.ndarray): 2D array containing neutrino CC & NC cross-section values, in cm^2
        nu_ixc (np.ndarray): 3D array containing neutrino integrated cross-section CDF values
        depth_max (float): Maximum column depth for neutrino propagation, in kmwe
        fac_nu (float): Rescaling factor for SM neutrino cross-sections

    Returns:
        part_type (integer): Type of outgoing particle. 0=neutrino; 1=charged lepton
        d_travel (float): Distance traveled until converted to charged lepton or total distance traveled by neutrino (if no conversion to charged lepton), in kmwe
        e_fin (float): Final neutrino energy, in GeV
    """
    part_type = 0 #int_typ=0CC, 1:NC
    col_depth_total = depth_max*1e5 #kmwe to cmwe
    e_nu = e_init
    e_fin = e_init
    x_0 = 0.0 #starting depth in g/cm^2
    d_travel = depth_max #added this in case there is aprogram and needed for breaking out when E<1e3
    #print('e_nu = ', e_nu)
    while e_nu > 1e3:

        r = np.random.random()
        int_depth = 1*transport.int_depth_nu(e_nu,nu_xc,fac_nu)
        x= -1*int_depth*np.log(r) # prob of interaction=exp(-x/int_depth)
        x_f = x_0 + x   #xf is tracking total column depth traveled

        if (x_f > col_depth_total): #total col depth exceeded
            return part_type,depth_max,e_fin

        x_0 = x_f
        int_type = transport.interaction_type_nu(e_nu,nu_xc,fac_nu) #CC or NC interaction

        if part_type == 0 and int_type == 1:
            part_type = 0
        elif part_type ==0 and int_type ==0:
            part_type = 1

        y = transport.find_y(e_nu,nu_ixc,int_type)

        e_fin = e_nu*(1-y) #Energy Transfer y = (E_i-E_f)/E_i

        if (part_type==1):
            d_travel = x_0*1e-5 #cmwe to kmwe
            return part_type,d_travel,e_fin

        e_nu = e_fin
        d_travel = x_0*1e-5 #cmwe to kmwe

        if e_nu <= 1e3:
            return part_type,d_travel,e_fin

    return part_type,d_travel,e_fin


def propagate_lep_water(e_init, xc_water, lep_ixc, alpha_water, beta_water, d_in, lepton, prop_type, cthi, Pi):
    """Propagates a charged lepton in water inside the Earth

    Args:
        e_init (float): Initial energy of the charged lepton, in GeV
        xc_water (np.ndarray): 2D array containing N_A/A*charged lepton-nucleon cross-section values in water, in cm^2/g
        lep_ixc (np.ndarray): 3D array containing charged lepton integrated cross-section CDF values in water
        alpha_water (np.ndarray): 1D array containing ionization energy loss values in water, in (GeV*cm^2)/g
        beta_water (np.ndarray): 2D array of beta values in water, in cm^2/g
        d_in (float): Maximum distance for charged lepton to propagate in water, in kmwe
        lepton (Integer): Type of charged lepton. 1=tau; 2=muon
        prop_type (integer): Type of energy loss propagation. 1=stochastic, 2=continuous
        cthi (float): costheta value obtained from tau EM interaction with rock
        Pi (float): Degree of Polarization obtained from tau EM interaction with rock

    Returns:
        part_id (integer): Type of outgoing charged lepton. 0=decayed; 1=not decayed; 2=don't count
        d_fin (float): Distance traveled before charged lepton decays or total distance traveled by charged lepton, in kmwe
        e_fin (float): Final energy of the charged lepton, in GeV
        pcthf (float): Final polarization after EM interaction of the tau lepton
    """
    e_min = 380*1e6 #Minimum tau energy, in GeV
    part_id = 1 #Start with tau that's obviously not decayed

    x_total = d_in*1e5 # kmwe to cmwe
    e_lep = e_init
    e_fin = e_init #in case the first interaction is too far
    x_0 = 0.0 #haven't gone anywhere yet; initiate tracker
    pcthf =0

    d_fin = d_in #temp flag
    if lepton == 1:
        m_le = 1.77682 #m_tau in GeV
        # change c_tau to whatever it is times 10^6
        c_tau = 8.793e-3 # c*lifetime, in cm, for taus (taken from PDB)
    else:
        m_le = 0.10565837550000001 #m-mu in GeV
        c_tau = 6.586384e4 #c*lifetime, in cm, for muons (taken from PDB 2020)

    if prop_type ==1:
        Pin = Pi
        theta_in = np.arccos(cthi)

        while (e_lep > e_min):

            if e_lep <= e_min:
                break #taken care of outside the while loop

            r = np.random.random()
            int_depth = transport.int_depth_lep(e_lep,xc_water,const.rho_water,m_le,c_tau)
            x = -1*int_depth*np.log(r) #basically the step size
            # prob of interaction = exp(-x/int_depth)

            x_f = x_0 + x #how far have we traveled here
            d_fin = x_f/1e5 #make sure it is not past the old number, in kmwe

            if (x_f >= x_total): #already past maximum depth but still a tau
                x_step = x_total - x_0 #backtrack one step

                alpha = transport.int_alpha(e_lep,alpha_water)
                beta = transport.int_beta(e_lep,beta_water,const.rho_water)

                e_fin = e_lep - (e_lep*beta + alpha)*x_step
                d_fin=d_in

                if (e_fin <= e_min): #tau has decayed
                    d_fin = d_in
                    e_fin = e_min
                    part_id = 2 # don't count

                pcthf = Pin*np.cos(theta_in)
                return part_id, d_fin, e_fin, pcthf

            x_0 = x_f #update x_0 and keep going
            alpha = transport.int_alpha(e_lep, alpha_water)
            beta = transport.int_beta(e_lep,beta_water,const.rho_water)

            e_int = e_lep - (e_lep*beta + alpha)*x #find some intermediate energy to get reasonable values of energy between inital and final energy, a la MUSIC

            if (e_int <= e_min):
                e_int = e_min

            e_avg = 10.0**((np.log10(e_lep)+np.log10(e_int))/2.0) #avg of 10^7 & 10^8 is 10^7.5

            alpha = transport.int_alpha(e_avg,alpha_water)
            beta = transport.int_beta(e_avg,beta_water,const.rho_water)

            e_int = transport.em_cont_part(e_lep,alpha,beta,x,m_le)

            if (e_int <= e_min): # below min energy
                break #taken care of outside while loop

            int_type = transport.interaction_type_lep(e_int,xc_water,const.rho_water,m_le,c_tau)

            if (int_type == 2): # tau decayed
                part_id = 0
                e_fin = e_int
                d_fin = x_f/1e5
                pcthf = Pin*np.cos(theta_in)
                return part_id,d_fin,e_fin,pcthf

            #tau didn't decay. Now how much energy does it have after interaction?

            y = transport.find_y(e_int,lep_ixc,int_type) #stochastic energy loss sampling

            #outgoing tau energy is old e_lep*(1-y)
            e_lep = e_int*(1.0 - y) #this is the energy for the next interaction
            e_fin = e_lep

            if (int_type ==5):
                Pout, theta_out = transport.polarization(y,Pin,theta_in)
                Pin = Pout
                theta_in = theta_out

        #Now outside the while loop, e_lep has to be <= e_min
        if (e_lep <= e_min):
            #d_fin = d_in # max distance in water
            d_fin = x_f/1e5 #temp flag
            e_fin = e_min
            part_id = 2 # don't count this
            pcthf=Pin*np.cos(theta_in)
            return part_id,d_fin,e_fin,pcthf

        else: # continuous energy loss
            #print('in else statement continuous energy loss')

            j_max = int(x_total/(const.step_size*const.rho_water)) # we will get close to exit

            #not sure if this should be j_max +1 or j_max +2 bc of fortran not 0 indexing
            for i in range(0,j_max +2): #takes care of the integer truncation issue
                if (e_lep < e_min): #taken care of outside the do while loop
                    break

                delta_x = const.step_size * const.rho_water #distance goes into decay
                x_f = x_0 + delta_x
                #does the particle decay over this distance?
                part_id = transport.idecay(e_lep,const.step_size,m_le,c_tau)

                if part_id ==0 : #we are all done
                    e_fin = e_lep
                    #d_fin = d_in
                    d_fin = x_f/1e5 #temp flag
                else: # find the new energy; assume alpha and beta are total values, not cut values
                    alpha = transport.int_alpha(e_lep,alpha_water)
                    beta = transport.int_beta(e_lep,beta_water,const.rho_water)

                    e_fin = e_lep - (e_lep*beta + alpha)*delta_x
                    d_fin = x_f/1e5
                    x_0 = x_f
                    e_lep = e_fin

                    if (i >= j_max):
                        break


            if i >= j_max:
                x_step = x_total - x_f #backtrack a little
                if x_step > 0.0: # last little energy loss
                    e_fin = e_lep - (e_lep * beta + alpha)*x_step #take care of that last little dx
                else:
                    if (e_fin <= e_min):
                        e_fin = e_min
                        d_fin = d_in
                        part_id = 2
                    elif e_fin > e_init:
                        e_fin = e_init #sanity check
            #outside the while e_lep has to be < e_min
            if e_lep <= e_min:
                d_fin = x_f/1e5
                e_fin = e_min
                part_id = 2 #don't count this
                return part_id,d_fin,e_fin,pcthf

#    return part_id, d_fin, e_fin, pcthf

def propagate_lep_rock(angle, e_init, xc_rock, lep_ixc, alpha_rock, beta_rock, d_entry,
                       d_in, xalong, cdalong, idepth, lepton, prop_type, earth_model):
    """Propagates a charged lepton in rock & iron inside the Earth

    Args:
        angle (Float): Earth emergence angle (beta)
        e_init (Float): Initial energy of charged lepton, in GeV.
        xc_rock (np.ndarray): 2D array containing N_A/A*charged lepton-nucleon cross-section values in rock, in cm^2/g.
        lep_ixc (np.ndarray): 3D array containing charged lepton integrated cross-section CDF values in rock.
        alpha_rock (np.ndarray): 1D array containing ionization energy loss values in rock, in (GeV*cm^2)/g.
        beta_rock (np.ndarray): 2D array of beta values in rock, in cm^2/g.
        d_entry (Float): Column depth along the chord for a given Earth emergence angle, in kmwe.
        d_in (Float): How much distance in rock/iron the charged lepton is supposed to travel, in kmwe.
        xalong (Float): 1D array containing distance in water, in km.
        cdalong (Float): 1D array containing column depth at xalong, in g/cm^2.
        idepth (Integer): Depth of water layer in km.
        lepton (Integer): Type of charged lepton. 1=tau; 2=muon.
        prop_type (Integer): Type of energy loss propagation. 1=stochastic, 2=continuous.
        earth_model (str): Earth density model, prem or ak135.

    Returns:
        part_id (Integer): Type of outgoing charged lepton. 0=decayed; 1=not decayed; 2=don't count.
        d_fin (Float): Distance traveled before charged lepton decays or total distance traveled by charged lepton, in kmwe
        e_fin (Float): Final charged lepton energy, in GeV.
        cthf (Float): Final costheta after EM interaction of the tau lepton.
        Pf (Float): Final degree of polarization after EM interaction of the tau lepton.

    """
    e_min = 1e3 #min tau energy, in GeV
    part_id = 1 #start with tau
    col_depth = d_entry*1e5 #how far in
    d_max = d_in*1e5 # how much to go, in cmwe
    e_lep = e_init
    x_0 = 0.0
    cnt = 0
    # x_total = d_in*1e5
    # rho = rho_rock #g/cm^3 USE FOR TESTING P_SURV FOR ROCK ONLY !
    if lepton ==1:
        m_le = 1.77682
        c_tau = 8.703e-3
    else:
        m_le = 0.10565837550000001
        c_tau = 6.586384e4

    if prop_type == 1 : #stochastic energy loss
        Pin = 1.0
        theta_in = 0.0

        while e_lep > e_min:
            cnt = cnt + 1
            x_interp = transport.cd2distd(xalong,cdalong,col_depth)
            r,rho = geometry.densityatx(x_interp,angle,idepth,earth_model)

            dy = np.random.random()
            int_len = transport.int_depth_lep(e_lep, xc_rock, rho, m_le, c_tau)

            x = -int_len * np.log(dy)
            col_depth = col_depth + x #update along trajectory, from the start of the chord
            x_f = x_0 + x #going 1D
            d_fin = x_f/1e5

            if x_f > d_max: #already past max depth
                #go to 30
                d_rem = d_max - x_0
                alpha = transport.int_alpha(e_lep, alpha_rock)
                beta = transport.int_beta(e_lep, beta_rock, rho)
                e_fin = e_lep - (e_lep*beta + alpha)*d_rem
                d_fin = d_max/1e5
                if e_fin > e_init:
                    e_fin = e_init
                #print('Passed Earth already')
                Pf = Pin
                cthf = np.cos(theta_in)

                return part_id,d_fin,e_fin,cthf,Pf

            x_0 = x_f #update x_0 and keep going
            alpha = transport.int_alpha(e_lep, alpha_rock)
            beta = transport.int_beta(e_lep, beta_rock, rho)
            e_int = e_lep - (e_lep*beta + alpha)*x #find some intermediate energy to get reasonable values of energy between inital and final energy, a la MUSIC

            if e_int <= e_min:
                e_int = e_min

            e_avg = 10.0**((np.log10(e_lep) + np.log10(e_int))/2) # does this work ?

            alpha = transport.int_alpha(e_avg, alpha_rock)
            beta = transport.int_beta(e_avg, beta_rock, rho)

            e_int = transport.em_cont_part(e_lep, alpha, beta, x, m_le) # get the continuous energy

            if e_int <= e_min:# is it below minimum energy now ?
                #print('below min energy')
                e_fin = e_int
                d_fin = d_max/1e5
                e_fin = e_min
                part_id = 0
                #print('too small to handle')
                #print('pcthf = ', Pin*np.cos(theta_in))
                Pf = Pin
                cthf = np.cos(theta_in)

                return part_id,d_fin,e_fin,cthf,Pf

            int_type = transport.interaction_type_lep(e_int, xc_rock, rho, m_le, c_tau)

            if int_type == 2: #tau has decayed
                #print('tau decay')
                part_id = 0
                e_fin = e_int
                Pf = Pin
                cthf = np.cos(theta_in)

                return part_id,d_fin,e_fin,cthf,Pf

            #tau didn't decay. Now how much energy does it have after interaction?
            y = transport.find_y(e_int, lep_ixc, int_type)

            #outgoing tau energy is old e_lep*(1-y)
            e_lep = e_int*(1.0-y) #this is the energy for the next ineraction
            e_fin = e_lep

            # This routine will run for only PN interaction
            if (int_type ==5):
                Pout, theta_out = transport.polarization(y,Pin,theta_in)
                Pin = Pout
                theta_in = theta_out
        #outside the while loop, e_lep has to be < e_min
        if e_lep <= e_min: #only continuous energy loss
            #print('ENERGY TOO LOW')
            d_fin = d_max/1e5
            e_fin = e_min
            part_id = 0 #dayed or no_count?? should be decayed
            Pf = Pin
            cthf = np.cos(theta_in)
            return part_id,d_fin,e_fin,cthf,Pf

    else: #continuous energy loss

        d_0 = 0.0
        delta_x = const.step_size # for now, not adaptive, distance into decay, cm; works for taus

        x_interp = transport.cd2distd(xalong, cdalong, col_depth) # find how far we are along the chord for a given beta
        r, rho = geometry.densityatx(x_interp, angle, idepth, earth_model) # find rho at x

        j_max = int(d_max/(rho*delta_x))

        if j_max == 0:
            e_fin = e_init
            d_fin = d_in
            return part_id,d_fin,e_fin,cthf,Pf

        for i in range(j_max +2):

            if e_lep < e_min:
                break

            x_interp = transport.cd2distd(xalong, cdalong, col_depth) # find how far we are along the chord for a given beta
            r, rho = geometry.densityatx(x_interp, angle, idepth, earth_model) # find rho at x

            delta_d = delta_x * rho
            x_0 = x_0 + delta_x
            d_0 = d_0 + delta_d

            #does the particle decay over this distance?
            part_id = transport.idecay(e_lep, delta_x, m_le, c_tau)

            if part_id == 0:
                e_fin = e_lep
                d_fin = d_0/1e5
                return part_id,d_fin,e_fin,cthf,Pf
            else: #find the new energy; assume alpha and beta are total values, not cut values
                alpha = transport.int_alpha(e_lep, alpha_rock)
                beta = transport.int_beta(e_lep, beta_rock, rho)
                e_lep = e_lep - (e_lep*beta + alpha)*delta_d
                d_fin = d_0/1e5 #updating the d_final
                e_fin = e_lep
                col_depth = d_entry*1e5 + d_0 # in order to update rho

                if (i >= j_max):
                    break # if so, break out of the loop
        if cnt >= j_max:
            if delta_x > 0.0:
                d_fin = d_in
                return part_id,d_fin,e_fin,cthf,Pf
            else:
                if e_fin <= e_min:
                    e_fin = e_min
                    d_fin = d_in
                    part_id = 2
                elif e_fin > e_init:
                    e_fin = e_init #sanity check

                return part_id,d_fin,e_fin,cthf,Pf

        if e_lep <= e_min: #only continuous energy loss
            d_fin = d_in
            e_fin = e_min
            part_id =2 #don't count this
            return part_id,d_fin,e_fin,cthf,Pf
        elif cnt >= j_max:
            d_fin = d_in
            return part_id,d_fin,e_fin,cthf,Pf
   # return part_id, d_fin, e_fin, cthf, Pf



def tau_thru_layers(angle,depth,d_water,depth_traj,e_lep_in,xc_water,xc_rock,lep_ixc_water,lep_ixc_rock,alpha_water,
                       alpha_rock,beta_water,beta_rock,xalong,cdalong,idepth,lepton,prop_type, earth_model):
    """

    Args:
        angle (Float): Earth emergence angle (beta), in degrees.
        depth (Float): Total column depth of the chord, in kmwe.
        d_water (Float): Column depth of the final layer of water (or full distance in water if only water layer), in kmwe.
        depth_traj (Float): Column depth along the chord for a given Earth emergence angle, in kmwe.
        e_lep_in (Float): Ingoing charged lepton energy, in GeV.
        xc_water (np.ndarray): 2D array containing N_A/A*charged lepton-nucleon cross-section values in water, in cm^2/g.
        xc_rock (np.ndarray): 2D array containing N_A/A*charged lepton-nucleon cross-section values in rock, in cm^2/g.
        lep_ixc_water (np.ndarray): 3D array containing charged lepton integrated cross-section CDF values in water.
        lep_ixc_rock (np.ndarray): 3D array containing charged lepton integrated cross-section CDF values in rock.
        alpha_water (np.ndarray): 1D array containing ionization energy loss values in water, in (GeV*cm^2)/g.
        alpha_rock (np.ndarray): 1D array containing ionization energy loss values in rock, in (GeV*cm^2)/g.
        beta_water (np.ndarray): 2D array of beta values in water, in cm^2/g.
        beta_rock (np.ndarray): 2D array of beta values in rock, in cm^2/g.
        xalong (np.ndarray): 1D array containing distance in water, in km.
        cdalong (np.ndarray): 1D array containing column depth at xalong, in g/cm^2.
        idepth (Integer): Depth of water layer in km.
        lepton (Integer): Type of charged lepton. 1=tau; 2=muon.
        prop_type (Integer): Type of energy loss propagation. 1=stochastic, 2=continuous.
        earth_model (str): Earth density model, prem or ak135.

    Returns:
        part_type (Integer): Type of outgoing charged lepton. 0=decayed, 1=not decayed.
        d_fin (Float): Distance traveled before charged lepton decays or total distance traveled by charged lepton, in kmwe.
        e_fin (Float): Outgoing charged lepton energy, in GeV.
        pcthf (FLoat): Final polarization after EM interaction of the tau lepton
    """
    d_fin = depth_traj
    col_depth = depth_traj*1e5 #g/cm^2
    e_lep = e_lep_in # so e_lep_in doesn't change
    e_fin = e_lep
    part_type = 1 #tau going in
    if e_lep < 1e3: #just in case
        part_type = 0
        d_fin = depth
        return part_type, d_fin,e_fin,None

    if (depth-depth_traj) < d_water:
        rho = const.rho_water
    else:
        x = transport.cd2distd(xalong, cdalong, col_depth)
        r, rho = geometry.densityatx(x, angle, idepth, earth_model) # find rho at x

        if rho <= 0.0:
            print('rho is 0')
        if rho <= 1.5 and r < 6365.0:
            print('rho too small')

    if rho > 1.5: # we aren't in water yet
        #print('WERE NOT IN WATER HELP')
        d_in = depth - depth_traj - d_water #propagate this far in rock

        part_type, d_f, e_fin, cthf, Pf = propagate_lep_rock(angle, e_lep, xc_rock, lep_ixc_rock, alpha_rock,
                                                             beta_rock, depth_traj, d_in, xalong, cdalong,
                                                             idepth, lepton,prop_type, earth_model)
        depth_traj = depth_traj + d_f

        if part_type == 1 and idepth != 0: #still a tau; added and clause on 3/18

            e_lep = e_fin
            d_in = d_water
            #depth_traj = depth_traj + d_f
            #d_fin = depth_traj
            Pi = Pf
            cthi = cthf

            part_type, d_f, e_fin, pcthf = propagate_lep_water(e_lep, xc_water, lep_ixc_water, alpha_water,
                                                               beta_water, d_in, lepton, prop_type, cthi, Pi)

            depth_traj = depth_traj + d_f
            d_fin = depth_traj

        else:
            pcthf = Pf*cthf
            d_fin = depth_traj
            return part_type,d_fin,e_fin,pcthf

    else:
        d_in = depth - depth_traj
        Pi = 1.0
        cthi = np.cos(0.0)

        part_type,d_f,e_fin,pcthf = propagate_lep_water(e_lep, xc_water, lep_ixc_water, alpha_water, beta_water, d_in, lepton,
                                                        prop_type, cthi, Pi)
        depth_traj = depth_traj + d_f
        d_fin = depth_traj
        return part_type,d_fin,e_fin,pcthf


    return part_type,d_fin,e_fin,pcthf

def distnu(r, ithird, Pin):
    """Determines the neutrino energy from tau decay. The energy fraction is determined by tau energy CDF
    approximated by leptonic decay channel. Approximate (good enough) for taus; exact for muons

    Args:
        r (Float): Random number
        ithird (IntegerL): Choice for neutrino -> charged lepton energy fraction selection
        Pin (Float): Tau's polarization

    Returns:
        dist_val (Float): Energy fraction , y = E_nu_tau/E_tau
    """
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
          lepton, fac_nu, prop_type, earth_model, Pin):
    """Regeneration loop, should take a pin and also throw out pout. pin will be used by distnu and the pout
       it throws will be used by the regen again as input in the single_stat() function.

    Args:
        angle (float): Earth emergence angle (beta), in degrees.
        e_lep (float): Incoming charged lepton energy, in GeV.
        depth (float): Total column depth of the chord, in kmwe.
        d_water (float): Column depth of the final layer of water (or full distance in water if only water layer), in kmwe.
        d_lep (float): Column depth along the chord for a given Earth emergence angle, in kmwe.
        nu_xc (np.ndarray): 2D array containing neutrino CC & NC cross-section values, in cm^2.
        nu_ixc (np.ndarray): 3D array containing neutrino integrated cross-section CDF values.
        ithird (integer): Choice for neutrino -> charged lepton energy fraction selection.
        xc_water (np.ndarray): 2D array containing N_A/A*charged lepton-nucleon cross-section values in water, in cm^2/g.
        xc_rock (np.ndarray): 2D array containing N_A/A*charged lepton-nucleon cross-section values in rock, in cm^2/g.
        ixc_water (np.ndarray): 3D array containing charged lepton integrated cross-section CDF values in water.
        ixc_rock (np.ndarray): 3D array containing charged lepton integrated cross-section CDF values in rock.
        alpha_water (np.ndarray): 1D array containing ionization energy loss values in water, in (GeV*cm^2)/g.
        alpha_rock (np.ndarray): 1D array containing ionization energy loss values in rock, in (GeV*cm^2)/g.
        beta_water (np.ndarray): 2D array of beta values in water, in cm^2/g.
        beta_rock (np.ndarray): 2D array of beta values in rock, in cm^2/g.
        xalong (np.ndarray): 1D array containing distance in water, in km.
        cdalong (np.ndarray): 1D array containing column depth at xalong, in g/cm^2.
        idepth (integer): Depth of water layer in km.
        lepton (integer): Type of charged lepton. 1=tau; 2=muon.
        fac_nu (float): Rescaling factor for SM neutrino cross-sections.
        prop_type (integer): Type of energy loss propagation. 1=stochastic, 2=continuous.
        earth_model (str): Earth density model, prem or ak135.
        Pin (float): tau's polarization input from tau passing thru layers

    Returns:
        part_type (integer): Type of outgoing particle. 0=neutrino; 3=exit.
        d_exit (float): Distance traveled before charged lepton decays or total distance traveled by charged lepton, in kmwe.
        e_fin (float): Final particle energy, in GeV.
        Pout (float): Tau's polarization output from regen function
    """
    r = np.random.random()
    frac = distnu(r,ithird,Pin) # that pola input goes here as input and cal. the neutrino energy
    e_nu = frac*e_lep

    d_left = depth - d_lep # this is how are the neutrino can go
    e_fin = e_nu #in case we need to exit
    part_type = 3 #in case we need to exit

    if d_left <= 0.0: # past the point of interactions allowed
        d_exit = depth
        Pout = Pin
        return part_type, d_exit, e_fin, Pout

    d_exit = d_lep #we are starting this far into the Earth with a neutrino
    int_part = 0 # starting with a neutrino with energy e_nu; change later to string

    int_part, dtr, etau2 = propagate_nu(e_nu, nu_xc, nu_ixc, d_left, fac_nu)
    #print('int part after propagate_nu', int_part)
    if int_part != 1: # neutrinos at the end
        d_exit = depth
        part_type = 0 #HLS =0
        e_fin = etau2 #final neutrino energy
        Pout = -2.0 #fake polarization value which will help us filter out the neutrino events
        return part_type, d_exit, e_fin, Pout

    #otherwise we have a tau
    d_lep = d_lep + dtr
    d_left = d_left - dtr
    #print('dleft', d_left)
    if d_left <= 0.0:
        d_exit = depth
        e_fin = etau2
        part_type = 0
        Pout = Pin
        return part_type, d_exit, e_fin, Pout

    #we have a tau with room to travel for tauthrulayers

    part_type, d_exit, e_fin, Pi = tau_thru_layers(angle, depth, d_water, d_lep, etau2, xc_water, xc_rock, ixc_water, ixc_rock,
         alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, idepth, lepton, prop_type, earth_model)

    Pout = Pi

    return part_type, d_exit, e_fin, Pout



