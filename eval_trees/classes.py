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

    def __init__(self):
        self.workdir = ""             # current work directory
        self.lastdir = ""             # last output_XXXXX directory
        self.lastdirnr = -1           # XXXX from lastdir
        self.ncpu = 1            
        self.noutput = 1              # how many output_XXXXX directories exist
        self.nout = 1                 # how many outputs we're gonna deal with. (Some might not have merger tree data)
        self.outputnrs = None         # numpy array of output numbers
        self.output_start = 0         # lowest snapshot number that we're dealing with (>= 1)
        self.z0 = 0                   # index of z=0 snapshot
        self.mth_main = 0             # mass threshold for main haloes
        self.mth_sub  = 0             # mass threshold for sub haloes
        self.sussing = False          # use sussing criteria
        self.use_npart_threshold = True # use a particle number threshold
        self.npart_thresh_sub = 200     # particle threshold if in use
        self.npart_thresh_main = 200    # particle threshold if in use

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
        # read in
        self.aexp        = np.zeros(par.noutput)
        self.unit_l      = np.zeros(par.noutput)
        self.unit_m      = np.zeros(par.noutput)
        self.unit_t      = np.zeros(par.noutput)
        self.unit_dens   = np.zeros(par.noutput)
        # to be computed
        self.redshift    = np.zeros(par.noutput) # z
        self.H           = np.zeros(par.noutput) # H(z)
        self.rho_crit    = np.zeros(par.noutput)
        self.times       = np.zeros(par.noutput) # t(z)
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
        #  self.y                       = [np.zeros(1) for i in range(par.noutput)] # descendant y
        #  self.z                       = [np.zeros(1) for i in range(par.noutput)] # descendant z
        self.v                      = [np.zeros(1) for i in range(par.noutput)] # descendant vel x
        #  self.vy                      = [np.zeros(1) for i in range(par.noutput)] # descendant vel y
        #  self.vz                      = [np.zeros(1) for i in range(par.noutput)] # descendant vel z

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
class jumper_data():
#======================
    """
    A class to store jumper data needed for flexible
    analysis later
    """

    def __init__(self, sd, sp, zd, zp, md, mp):
        
        self.snapshot_desc = sd
        self.snapshot_prog = sp
        self.z_desc = zd
        self.z_prog = zp
        self.mass_desc = md
        self.mass_prog = mp
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

        self.hlogM = np.zeros(1)      # halo dlogM/dlogt
        self.hlogM_free = 0           # last used index for hlogM
        self.shlogM = np.zeros(1)     # subhalo dlogM/dlogt
        self.shlogM_free = 0          # last used index for shlogM
        self.logM = np.zeros(1)       # all dlogM/dlogt
        self.logM_free = 0            # last used index for logM

        self.mbl = None     # main branch length
        self.nbr = None     # number of branches

        self.njumpers = 0   # total number of jumpers
        self.npruned = 0    # number of pruned trees


        self.branch_bins = [100, 500, 1000]                                             # particle numbers for bins of main branch lengths
        npartbins = len(self.branch_bins) + 1

        self.branchlengths = [np.zeros(1, dtype=np.int) for b in range(npartbins)]          # lengths of main branches, divided into particle bins
        self.branchlen_free = [0 for b in range(npartbins)]                                 # first free index for every bin
        self.branchlengths_sub = [np.zeros(1, dtype=np.int) for b in range(npartbins)]      # lengths of main branches, divided into particle bins
        self.branchlen_free_sub = [0 for b in range(npartbins)]                             # first free index for every bin
        self.branchlengths_main = [np.zeros(1, dtype=np.int) for b in range(npartbins)]     # lengths of main branches, divided into particle bins
        self.branchlen_free_main = [0 for b in range(npartbins)]                            # first free index for every bin

        self.nbranches = [np.zeros(1, dtype=np.int) for b in range(npartbins)]              # number of branches, divided into particle bins
        self.nbranch_free = [0 for b in range(npartbins)]                                   # first free index for every bin
        self.nbranches_sub = [np.zeros(1, dtype=np.int) for b in range(npartbins)]          # number of branches, divided into particle bins
        self.nbranch_free_sub = [0 for b in range(npartbins)]                               # first free index for every bin
        self.nbranches_main = [np.zeros(1, dtype=np.int) for b in range(npartbins)]         # number of branches, divided into particle bins
        self.nbranch_free_main = [0 for b in range(npartbins)]                              # first free index for every bin

        self.displacements = np.zeros(1, dtype=np.float)        # all displacement statistics
        self.displacement_free = 0                              # first free index

        self.dirprogs = np.zeros(1)         # number of direct progenitors
        self.dirprog_free = 0               # last used index for direct progenitors

        self.jumper_results = []            # store jumper_data objects here

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
    def add_branch_length_sub(self, val, npart):
    #---------------------------------------------
        """
        Add new branch length, put it in the right bin based on npart
        for subhaloes only!
        """


        # first find the right bin
        b = 0
        while npart > self.branch_bins[b]:
            b+=1
            if b == len(self.branch_bins):
                break

        if self.branchlen_free_sub[b] == self.branchlengths_sub[b].shape[0]:
            self.branchlengths_sub[b].resize(self.branchlengths_sub[b].shape[0]+10000)

        self.branchlengths_sub[b][self.branchlen_free_sub[b]] = val
        self.branchlen_free_sub[b] += 1

        return


    #---------------------------------------------
    def add_branch_length_main(self, val, npart):
    #---------------------------------------------
        """
        Add new branch length, put it in the right bin based on npart
        for mainhaloes only!
        """

        # first find the right bin
        b = 0
        while npart > self.branch_bins[b]:
            b+=1
            if b == len(self.branch_bins):
                break


        if self.branchlen_free_main[b] == self.branchlengths_main[b].shape[0]:
            self.branchlengths_main[b].resize(self.branchlengths_main[b].shape[0]+10000)

        self.branchlengths_main[b][self.branchlen_free_main[b]] = val
        self.branchlen_free_main[b] += 1

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


    #---------------------------------------------
    def add_nr_branches_sub(self, val, npart):
    #---------------------------------------------
        """
        Add new branch length, put it in the right bin based on npart
        only for subhaloes!
        """

        # first find the right bin
        b = 0
        while npart > self.branch_bins[b]:
            b+=1
            if b == len(self.branch_bins):
                break


        if self.nbranch_free_sub[b] == self.nbranches_sub[b].shape[0]:
            self.nbranches_sub[b].resize(self.nbranches_sub[b].shape[0]+10000)

        self.nbranches_sub[b][self.nbranch_free_sub[b]] = val
        self.nbranch_free_sub[b] += 1

        return


    #---------------------------------------------
    def add_nr_branches_main(self, val, npart):
    #---------------------------------------------
        """
        Add new branch length, put it in the right bin based on npart
        only for main haloes!
        """

        # first find the right bin
        b = 0
        while npart > self.branch_bins[b]:
            b+=1
            if b == len(self.branch_bins):
                break


        if self.nbranch_free_main[b] == self.nbranches_main[b].shape[0]:
            self.nbranches_main[b].resize(self.nbranches_main[b].shape[0]+10000)

        self.nbranches_main[b][self.nbranch_free_main[b]] = val
        self.nbranch_free_main[b] += 1

        return




    #---------------------------------------
    def add_halo_displacement(self, val):
    #---------------------------------------
        if self.displacement_free == self.displacements.shape[0]:
            self.displacements.resize(self.displacements.shape[0]+10000)

        self.displacements[self.displacement_free] = val
        self.displacement_free += 1
        return







    #---------------------------------------------
    def add_dirprog(self, val):
    #---------------------------------------------
        """
        Add number of direct progenitors of a halo in this snapshot
        """

        if self.dirprog_free == self.dirprogs.shape[0]:
            self.dirprogs.resize(self.dirprogs.shape[0]+10000)

        self.dirprogs[self.dirprog_free] = val
        self.dirprog_free += 1

        return





    #---------------------------------------
    def add_halo_logM(self, val):
    #---------------------------------------
        if self.hlogM_free == self.hlogM.shape[0]:
            self.hlogM.resize(self.hlogM.shape[0]+10000)

        self.hlogM[self.hlogM_free] = val
        self.hlogM_free += 1
        return

    #---------------------------------------
    def add_subhalo_logM(self, val):
    #---------------------------------------
        if self.shlogM_free == self.shlogM.shape[0]:
            self.shlogM.resize(self.shlogM.shape[0]+10000)

        self.shlogM[self.shlogM_free] = val
        self.shlogM_free += 1
        return

    #---------------------------------------
    def add_any_logM(self, val):
    #---------------------------------------
        if self.logM_free == self.logM.shape[0]:
            self.logM.resize(self.logM.shape[0]+10000)

        self.logM[self.logM_free] = val
        self.logM_free += 1
        return

    #----------------------------------------------------------------------
    def add_jumper_data(self, snapshot_desc, snapshot_prog, z_desc, z_prog, 
                        mass_desc, mass_prog):
    #----------------------------------------------------------------------
        self.jumper_results.append(jumper_data(snapshot_desc, snapshot_prog, z_desc, z_prog, mass_desc, mass_prog))
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





