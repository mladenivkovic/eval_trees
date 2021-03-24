#!/usr/bin/python

#------------------------------
# Cosmology related functions
#------------------------------

import numpy as np



#=========================================
def compute_cosmo_quantities(p, c, sd):
#=========================================
    """
    Compute times and Hubble parameter given the expansion factor. 
    Then compute rho_crit at given aexp in M_Sol/Mpc^3 and the 
    cosmological time at that point.

        p:  params object
        c:  constants object
        sd: snapshotdata object
    """

    print("Computing cosmo stuff")

    sd.times = np.zeros(sd.aexp.shape[0], dtype='float')
    H = np.zeros(sd.aexp.shape[0], dtype='float')

    a_out, t_out, H_out = friedman(sd.aexp.min(), c)

    i = 1
    for j,a in enumerate(sd.aexp):
        while (a_out[i] > a) and (i < a_out.shape[0]-1):
            i += 1
        sd.times[j] = t_out[i]*(a-a_out[i-1])/(a_out[i]-a_out[i-1]) + \
            t_out[i-1]*(a-a_out[i])/(a_out[i-1]-a_out[i])
        sd.H[j] = H_out[i]*(a-a_out[i-1])/(a_out[i]-a_out[i-1]) + \
            H_out[i-1]*(a-a_out[i])/(a_out[i-1]-a_out[i])


    sd.redshift = 1.0/sd.aexp - 1

    #----------------------
    # find where to start
    #----------------------
    # find sd.redshift value closest to 0
    # (you need values of later sd.times (sd.times z<0) to check for jumpers correctly)
    p.z0 = np.argmin(np.absolute(sd.redshift))
    print( "z = 0 snapshot is", p.outputnrs[p.z0])


    #-----------------------
    # get physical units
    #-----------------------
    alpha = c.H0 * 1e5/c.Mpc*c.Gyr # 1e5: cm->km

    sd.times *= 1.0/alpha # get sd.times in Gyrs: sd.times were calculated in units of H_0

    # turn sd.time around: t=0 at start to get cosmic sd.time
    timestart = -sd.times[-1]+sd.times[p.z0]
    sd.times += timestart


    sd.H *= alpha # get H in Gyrs^-1
    # sanity check: H0 in units [Gyrs^-1] = 0.07195/Gyr
    sd.rho_crit = 3 * sd.H**2 / (8 * np.pi * c.G)    
    
    
    return





#==============================
def friedman(axp_min, c):
#==============================
    """
    Integrate friedman equation to get table of
    expansion factors and times.
    Gives a in units of H0.
    See ramses/utils/f90/utils.f90/subroutine friedman

        axp_min: smallest aexp
        c:       constants object
    """

    def dadtau(axp_tau):
        dadtau = axp_tau**3 * (c.omega_m + c.omega_b + c.omega_l*axp_tau**3 + c.omega_k * axp_tau)
        return np.sqrt(dadtau)

    def dadt(axp_t):
        dadt = 1/axp_t * (c.omega_m + c.omega_b + c.omega_l*axp_t**3 + c.omega_k*axp_t)
        return np.sqrt(dadt)


    epsilon = 1e-4 # tolerance

    axp_tau = 1.0 # expansion factor
    axp_t = 1.0   # expansion factor
    tau = 0       # conformal sd.time
    t = 0         # look-back sd.time

    a_out = np.zeros(1000000, dtype='float')
    t_out = np.zeros(1000000, dtype='float')
    #  tau_out = np.zeros(1000000, dtype='float')
    H_out = np.zeros(1000000, dtype='float')

    i = 0
    while axp_tau >= axp_min or axp_t >= axp_min:
        dtau = epsilon * axp_tau / dadtau(axp_tau)
        axp_tau_pre = axp_tau - dadtau(axp_tau)*dtau/2
        axp_tau = axp_tau - dadtau(axp_tau_pre)*dtau
        tau = tau - dtau

        dt = epsilon * axp_t / dadt(axp_t)
        axp_t_pre = axp_t - dadt(axp_t)*dt/2
        axp_t = axp_t - dadt(axp_t_pre)*dt
        t = t - dt

        if (abs((t - t_out[i])/t) >= 1e-5):
            a_out[i] = axp_tau
            H_out[i] = dadtau(axp_tau)/axp_tau
            t_out[i] = t
            #  tau_out[i] = tau

            i += 1

    a_out[i] = axp_tau
    t_out[i] = t
    #  tau_out[i] = tau
    H_out[i] = dadtau(axp_tau)/axp_tau
    i+=1

    a_out = a_out[:i]
    t_out = t_out[:i]
    #  tau_out = tau_out[:i]
    H_out = H_out[:i]


    return a_out, t_out, H_out




def compute_R200(M, rhoC):
    """
    Compute R200, radius that contains 200 times the
    critical density, for given halo mass M
    """

    return (0.75 / np.pi * M/ rhoC)**(1./3)
