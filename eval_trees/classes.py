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
        self.progenitors             = [np.zeros(1) for i in range(par.noutput)] # progenitor IDs
        self.descendants             = [np.zeros(1) for i in range(par.noutput)] # descendant IDs
        self.progenitor_outputnrs    = [np.zeros(1) for i in range(par.noutput)] # snapshot number of progenitor
        self.mass                    = [np.zeros(1) for i in range(par.noutput)] # descendant mass
        self.npart                   = [np.zeros(1) for i in range(par.noutput)] # descendant npart exclusive
        self.x                       = [np.zeros(1) for i in range(par.noutput)] # descendant x
        self.y                       = [np.zeros(1) for i in range(par.noutput)] # descendant y
        self.z                       = [np.zeros(1) for i in range(par.noutput)] # descendant z
        self.vx                      = [np.zeros(1) for i in range(par.noutput)] # descendant vel x
        self.vy                      = [np.zeros(1) for i in range(par.noutput)] # descendant vel y
        self.vz                      = [np.zeros(1) for i in range(par.noutput)] # descendant vel z

        self.is_halo                 = [np.zeros(1) for i in range(par.noutput)] # descendant is halo or not
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
        self.hmg = np.zeros(1)      # halo mass growth
        self.hmg_free = 0           # last used index for hmg
        self.hmf = np.zeros(1)      # halo mass fluctuation
        self.hmf_free = 0           # last used index for hmf
        self.shmg = np.zeros(1)     # subhalo mass growth
        self.shmg_free = 0          # last used index for shmg
        self.shmf = np.zeros(1)     # subhalo mass fluctuation
        self.shmf_free = 0          # last used index for shmf
        self.mg = np.zeros(1)       # all mass growth
        self.mg_free = 0            # last used index for mg
        self.mf = np.zeros(1)       # all mass fluctuation
        self.mf_free = 0            # last used index for mf

        self.mbl = None     # main branch length
        self.nbr = None     # number of branches

        self.njumpers = 0   # total number of jumpers
        self.npruned = 0    # number of pruned trees

        self.branch_bins = [100, 500, 1000]                             # particle numbers for bins of main branch lengths
        npartbins = len(self.branch_bins) + 1
        self.branchlengths = [np.zeros(1) for b in range(npartbins)]    # lengths of main branches, divided into particle bins
        self.branchlen_free = [0 for b in range(npartbins)]             # first free index for every bin
        self.nbranches = [np.zeros(1) for b in range(npartbins)]        # number of branches, divided into particle bins
        self.nbranch_free = [0 for b in range(npartbins)]               # first free index for every bin
        return



    #---------------------------------------
    def add_halo_growth(self, val):
    #---------------------------------------
        if self.hmg_free == self.hmg.shape[0]:
            self.hmg.resize(self.hmg.shape[0]+10000)

        self.hmg[self.hmg_free] = val
        self.hmg_free += 1
        return

    #---------------------------------------
    def add_subhalo_growth(self, val):
    #---------------------------------------
        if self.shmg_free == self.shmg.shape[0]:
            self.shmg.resize(self.shmg.shape[0]+10000)

        self.shmg[self.shmg_free] = val
        self.shmg_free += 1
        return

    #---------------------------------------
    def add_any_growth(self, val):
    #---------------------------------------
        if self.mg_free == self.mg.shape[0]:
            self.mg.resize(self.mg.shape[0]+10000)

        self.mg[self.mg_free] = val
        self.mg_free += 1
        return



    #---------------------------------------
    def add_halo_fluct(self, val):
    #---------------------------------------
        if self.hmf_free == self.hmf.shape[0]:
            self.hmf.resize(self.hmf.shape[0]+10000)

        self.hmf[self.hmf_free] = val
        self.hmf_free += 1
        return

    #---------------------------------------
    def add_subhalo_fluct(self, val):
    #---------------------------------------
        if self.shmf_free == self.shmf.shape[0]:
            self.shmf.resize(self.shmf.shape[0]+10000)

        self.shmf[self.shmf_free] = val
        self.shmf_free += 1
        return

    #---------------------------------------
    def add_any_fluct(self, val):
    #---------------------------------------
        if self.mf_free == self.mf.shape[0]:
            self.mf.resize(self.mf.shape[0]+10000)

        self.mf[self.mf_free] = val
        self.mf_free += 1
        return


    #---------------------------------------------
    def add_branch_length(self, val, npart):
    #---------------------------------------------
        """
        Add new branch length, put it in the right bin based on npart
        """

        # first find the right bin
        b = 0
        while npart > self.branch_bins[b]:
            b+=1
            if b == len(self.branch_bins):
                break


        if self.branchlen_free[b] == self.branchlengths[b].shape[0]:
            self.branchlengths[b].resize(self.branchlengths[b].shape[0]+10000)

        self.branchlengths[b][self.branchlen_free[b]] = val
        self.branchlen_free[b] += 1

        return



    #---------------------------------------------
    def add_nr_branches(self, val, npart):
    #---------------------------------------------
        """
        Add new branch length, put it in the right bin based on npart
        """

        # first find the right bin
        b = 0
        while npart > self.branch_bins[b]:
            b+=1
            if b == len(self.branch_bins):
                break


        if self.nbranch_free[b] == self.nbranches[b].shape[0]:
            self.nbranches[b].resize(self.nbranches[b].shape[0]+10000)

        self.nbranches[b][self.nbranch_free[b]] = val
        self.nbranch_free[b] += 1

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





