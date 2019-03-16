#!/usr/bin/env python3


#===================================
# CLASSES
# data containers for various things
#===================================



import numpy as np



#===================
class params():
#===================
    """
    Global parameters to be stored
    """

    def __init__(self, mthresh_main, mthresh_sub):
       self.workdir = ""             # current work directory
       self.lastdir = ""             # last output_XXXXX directory
       self.lastdirnr = -1           # XXXX from lastdir
       self.ncpu = 1            
       self.noutput = 1              # how many output_XXXXX directories exist
       self.nout = 1                 # how many outputs we're gonna deal with. (Some might not have merger tree data)
       self.outputnrs = None         # numpy array of output numbers
       self.output_start = 0         # lowest snapshot number that we're dealing with (>= 1)
       self.z0 = 0                   # index of z=0 snapshot
       self.mth_main = mthresh_main  # mass threshold for main haloes
       self.mth_sub  = mthresh_sub   # mass threshold for sub haloes
       return





#======================
class snapshotdata():
#======================
    """
    Snapshot specific data
    """

    def __init__(self, par):
        """
        par: params object
        """
        self.aexp        = np.zeros(par.noutput)
        self.H           = np.zeros(par.noutput)
        self.rho_crit    = np.zeros(par.noutput)
        self.times       = np.zeros(par.noutput)
        self.unit_l      = np.zeros(par.noutput)
        self.unit_m      = np.zeros(par.noutput)
        self.unit_t      = np.zeros(par.noutput)
        self.unit_dens   = np.zeros(par.noutput)
        return
 



#======================
class mtreedata():
#======================
    """ 
    mergertree data lists 
    """

    def __init__(self, par):
        """
        par: params object
        """
        self.progenitors             = [np.zeros(1) for i in range(par.noutput)]
        self.descendants             = [np.zeros(1) for i in range(par.noutput)]
        self.progenitor_outputnrs    = [np.zeros(1) for i in range(par.noutput)]
        self.mass                    = [np.zeros(1) for i in range(par.noutput)]
        self.npart                   = [np.zeros(1) for i in range(par.noutput)]
        self.x                       = [np.zeros(1) for i in range(par.noutput)]
        self.y                       = [np.zeros(1) for i in range(par.noutput)]
        self.z                       = [np.zeros(1) for i in range(par.noutput)]
        self.vx                      = [np.zeros(1) for i in range(par.noutput)]
        self.vy                      = [np.zeros(1) for i in range(par.noutput)]
        self.vz                      = [np.zeros(1) for i in range(par.noutput)]


        self.nhalosmax = 0
        return





#======================
class halodata():
#======================
    """
    Data from halo_XXXXX.txtYYYYY 
    """
    def __init__(self, par):
        """
        par: params object
        """

        self.haloid = [np.zeros(1) for i in range(par.noutput)]
        return






#======================
class constants():
#======================
    """
    Class holding constants.
    """

    def __init__(self):
        self.Mpc     = 3.086e24          # cm
        self.M_Sol   = 1.98855e33        # g
        self.Gyr     = (24*3600*365*1e9) # s
        self.G       = 4.492e-15         # Mpc^3/(M_sol Gyr^2)

        self.H0      = 70.4
        self.omega_m = 0.272
        self.omega_l = 0.728
        self.omega_k = 0.0
        self.omega_b = 0.0

        return





#======================
class results():
#======================
    """
    Store results
    """

    def __init__(self):
        self.hmg = None     # halo mass growth
        self.hmf = None     # halo mass fluctuation
        self.shmg = None    # subhalo mass growth
        self.shmf = None    # subhalo mass fluctuation

        self.mbl = None     # main branch length
        self.nbr = None     # number of branches

        self.njumpers = 0   # total number of jumpers
        self.npruned = 0    # number of pruned trees
        return






#========================
def init_lists(par):
#========================
    """
    Initialize list for snapshot and halo data
    """

    sd = snapshotdata(par)
    mtd = mtreedata(par)

    return sd, mtd





