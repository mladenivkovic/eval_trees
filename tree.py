#!/usr/bin/env python3


import numpy as np
from cosmo import compute_R200


#=========================================
def count_statistics_at_z0(p, r, mtd):
#=========================================
    """
    Count number of clumps and their median particle
    number at z = 0

        p:      params object
        r:      results object
        mtd:    mtreedata object
    """

    print("Gathering statistics at z = 0")

    proper_clumps = mtd.descendants[p.z0] > 0
    r.clumps_at_z0 = mtd.descendants[p.z0][proper_clumps].shape[0]
    r.median_clump_particlecount_at_z0 = np.median(mtd.npart[p.z0][proper_clumps])

    return





#=========================================
def count_pruned_trees(p, r, mtd, sd):
#=========================================
    """
    Count how many trees have been pruned.

        p:      params object
        r:      results object
        mtd:    mtreedata object
        sd:     snapshotdata object
    """

    print("Counting pruned trees")

    for out in range(p.z0+1, p.nout):
        absdescs = np.abs(mtd.descendants[out])
        absprogs = np.abs(mtd.progenitors[out-1])
        for i, de in enumerate(absdescs):
            # if it's a jumper, ignore
            if de == 0:
                continue
            if de not in absprogs:
                r.add_pruned_data(mtd.npart[out][i], mtd.is_halo[out][i])

    return



#=========================================
def count_total_descendants(p, mtd):
#=========================================
    """
    Count how many descendants are in simulation

        p:      params object
        r:      results object
        mtd:    mtreedata object
        sd:     snapshotdata object
    """

    total_descendants = 0

    for out in range(p.nout):
        
        total_descendants += mtd.descendants[out][mtd.descendants[out] > 0].shape[0]

    print("Total unique descendants in simulation:", total_descendants)

    return






#================================================
def clean_jumpers_only_after_z0(p, r, mtd, sd):
#================================================
    """
    Clean up jumpers:
    Search for jumpers in snapshots after z=0,
    then remove them in the snapshots with z > 0

        p:      params object
        r:      results object
        htd:    mtreedata object
        sd:     snapshotdata object
    """

    if p.z0==0:
        print('No snapshots with z < 0 found. Skipping cleaning jumpers')
        return


    nreplaced = 0

    for out in range(0, p.z0):
        for i, pr in enumerate(mtd.progenitors[out]):
            if pr < 0:
                # Subtract 1 here from snapind:
                # progenitor_outputnrs gives the snapshot number where the
                # jumper was a descendant for the last time
                # so you need to overwrite the merging one snapshot later,
                # where the clump is the progenitor
                snapind = _get_snap_ind(p, mtd.progenitor_outputnrs[out][i]) - 1

                #  print('found jumper ', pr, 'in snapshot', p.outputnrs[out], 'from snapshot', mtd.progenitor_outputnrs[out][i])
                #  print('searching for progenitor in', p.outputnrs[snapind])
                jumpind = np.where(mtd.progenitors[snapind]==-pr)[0] # where returns tuple

                if (jumpind.shape[0] != 1):
                    print('Got jumpind.shape=', jumpind.shape)
                    print('Something went wrong, exiting')
                    raise IndexError

                # replace jumper's descendant with 0
                if (mtd.descendants[snapind][jumpind] > 0):
                    print("Descendant > 0, not merger, what happened?")
                    raise IndexError
                else:
                    mtd.descendants[snapind][jumpind] = 0
                    nreplaced += 1

    print("Cleaned out", nreplaced, "jumpers")

    return








#=========================================
def clean_jumpers(p, r, mtd, sd):
#=========================================
    """
    Clean up jumpers:
    Search for jumpers everywhere, then remove them in the 
    snapshots with z > 0 where they appear to have merged

        p:      params object
        r:      results object
        mtd:    mtreedata object
        sd:     snapshotdata object
    """

    nreplaced = 0

    for out in range(p.nout):
        for i, pr in enumerate(mtd.progenitors[out]):
            if pr < 0:
                # Subtract 1 here from snapind:
                # progenitor_outputnrs gives the snapshot number where the
                # jumper was a descendant for the last time
                # so you need to overwrite the merging one snapshot later,
                # where the clump is the progenitor
                snapind = _get_snap_ind(p, mtd.progenitor_outputnrs[out][i]) - 1

                jumpind = np.where(mtd.progenitors[snapind]==-pr)[0] # where returns tuple

                if (jumpind.shape[0] != 1):
                    print('Got jumpind.shape=', jumpind.shape)
                    print('Something went wrong, exiting')
                    raise IndexError

                # replace jumper's descendant with 0
                if (mtd.descendants[snapind][jumpind] > 0):
                    print("Descendant > 0, not merger, what happened?")
                    raise IndexError
                else:

                    # store jumper data for analysis
                    snap_desc = p.outputnrs[out]
                    snap_prog = mtd.progenitor_outputnrs[out][i]
                    z_desc = sd.redshift[out]
                    z_prog = sd.redshift[snapind]
                    m_desc = mtd.mass[out][i]
                    desc_is_halo = mtd.is_halo[out][i]

                    snapind_where_prog_was_last_desc = snapind + 1
                    progind = np.where(mtd.descendants[snapind_where_prog_was_last_desc] == -pr)[0]
                    m_prog = np.asscalar(mtd.mass[snapind_where_prog_was_last_desc][progind])
                    prog_is_halo = np.asscalar(mtd.is_halo[snapind_where_prog_was_last_desc][progind])
                    r.add_jumper_data(snap_desc, snap_prog, z_desc, z_prog, m_desc, m_prog, desc_is_halo, prog_is_halo)

                    # now replace jumper's descendant with 0
                    mtd.descendants[snapind][jumpind] = 0

                    # count deletion
                    nreplaced += 1

    print("Cleaned out", nreplaced, "jumpers")

    return









#=============================================
def determine_if_halo(p, mtd, hd):
#=============================================
    """
    For all descendants, find whether they are haloes or
    subhaloes.
    
        p:      params object
        mtd:    mergertree data
        hd:     halodata
    """

    print("Checking which descendants are main haloes")

    for i,desc in enumerate(mtd.descendants):
        is_halo = np.zeros(desc.shape, dtype=np.bool) # zeroes = False
        d_inds = np.argsort(np.abs(desc))

        halos_sorted = np.sort(hd.haloid[i])

        id = 0
        ih = 0
        while id < desc.shape[0]:
            if abs(desc[d_inds[id]]) == halos_sorted[ih]:
                is_halo[d_inds[id]] = 1
                id += 1
                #  ih += 1 # don't do this: multiple progenitors may have same descendant, albeit a negative one
            elif abs(desc[d_inds[id]]) > halos_sorted[ih]:
                ih += 1
            else:
                id += 1

        mtd.is_halo[i] = is_halo.copy()


    return





#====================================
def get_mass_thresholds(p, mtd, hd):
#====================================
    """
    Compute mass thresholds:
    include 1000 main haloes and 200 subhaloes at z=0 for sussing

        p:      params object
        mtd:    mergertree data
        hd:     halodata
    """

    # use same as in Avila paper
    #------------------------------------

    if p.sussing:
        hids = hd.haloid[p.z0]
        ids = mtd.descendants[p.z0]
        is_halo = mtd.is_halo[p.z0]

        hmasses = mtd.mass[p.z0][is_halo]
        shmasses = mtd.mass[p.z0][np.logical_not(is_halo)]

        hmasses = np.sort(hmasses)
        shmasses = np.sort(shmasses)

        try:
            p.mth_main = hmasses[-1000]
        except IndexError:
            print("There aren't 1000 haloes in z=0. Setting no mass threshold.")
            p.mth_main = 0
        try:
            p.mth_sub = shmasses[-200]
        except IndexError:
            print("There aren't 200 subhaloes in z=0. Setting no mass threshold.")
            p.mth_sub = 0

    else:

        # use fixed number of particles
        #------------------------------------

        if p.use_npart_threshold:

            mp  = mtd.mass[p.z0][0]/mtd.npart[p.z0][0]

            p.mth_sub = p.npart_thresh_sub * mp
            p.mth_main = p.npart_thresh_main * mp


        # use fixed mass
        #------------------------------------
        else:
            p.mth_sub = 3e11
            p.mth_main = 3e11


    print("Main halo mass threshold is: {0:.3e}".format(p.mth_main))
    print("Subhalo mass threshold is: {0:.3e}".format(p.mth_sub))
    print("Particle mass is: {0:.3e}".format(mtd.mass[p.z0][0]/mtd.npart[p.z0][0]))

    return







#=====================================================
def get_main_branch_lengths(p, r, mtd, hd, sd):
#=====================================================
    """
    Compute the lenghts of the main branches for the 
    haloes and subhaloes that satisfy the mass threshold

        p:      params object
        r:      results object
        mtd:    mtreedata object
        hd:     clumpdata object
        sd:     snapshotdata object
    """

    print("Computing Main Branch Lengths")

    for c, clump in enumerate(mtd.descendants[p.z0]):

        d_npart = mtd.npart[p.z0][c]

        if clump > 0:
            # clump = 0 is removed jumper;
            # clump < 0 is merger; we're only taking main branch here

            dind = c
            desc_snap_ind = p.z0
            desc = mtd.descendants[p.z0][dind]
            prog = mtd.progenitors[p.z0][dind]
            desc_is_halo = mtd.is_halo[p.z0][dind]


            if p.sussing:
                # skip if threshold is not satisfied
                if desc_is_halo:
                    if mtd.mass[p.z0][dind] < p.mth_main:
                        continue
                else:
                    if mtd.mass[p.z0][dind] < p.mth_sub:
                        continue


            last_snap = p.z0 # the last snapshot this clump appeared in

            # now descend down the main branch
            while prog != 0:

                if prog > 0:    # if progenitor isn't jumper
                    # find progenitor's index in previous snapshot
                    p_snap_ind = desc_snap_ind + 1
                    pind = np.where(mtd.descendants[p_snap_ind]==prog)[0]

                elif prog < 0:
                    p_snap_ind = _get_snap_ind(p, mtd.progenitor_outputnrs[desc_snap_ind][dind])
                    pind = np.where(mtd.descendants[p_snap_ind]==-prog)[0]

                if pind.shape[0] != 1:
                    print('got pind.shape=', pind.shape)
                    print('something went wrong, exiting')
                    raise IndexError

                # prepare for next round
                desc_snap_ind = p_snap_ind
                dind = pind
                last_snap = p_snap_ind
                prog = mtd.progenitors[p_snap_ind][pind]

            # now put the result in the appropriate bin
            if (not p.sussing) or mtd.mass[p.z0][c] > p.mth_main:
                r.add_branch_length(last_snap - p.z0, mtd.npart[p.z0][c])

            if mtd.is_halo[p.z0][c]:
                r.add_branch_length_main(last_snap - p.z0, mtd.npart[p.z0][c])
            else:
                r.add_branch_length_sub(last_snap - p.z0, mtd.npart[p.z0][c])

            
            if p.do_short_branch_investigation:

                if d_npart > 1000 and last_snap - p.z0 < 6:
                    print(
                            "--- Short branches: clumpID", clump, 
                            "npart", d_npart, 
                            "halo?", mtd.is_halo[p.z0][c], 
                            "len", last_snap - p.z0, 
                            "last desc", mtd.descendants[desc_snap_ind][dind], 
                            "@", p.outputnrs[last_snap]
                            )

                    # repeat tree walk and print out each entry
                    dind = c
                    desc_snap_ind = p.z0
                    desc = mtd.descendants[p.z0][dind]
                    prog = mtd.progenitors[p.z0][dind]
                    desc_is_halo = mtd.is_halo[p.z0][dind]


                    if p.sussing:
                        # skip if threshold is not satisfied
                        if desc_is_halo:
                            if mtd.mass[p.z0][dind] < p.mth_main:
                                continue
                        else:
                            if mtd.mass[p.z0][dind] < p.mth_sub:
                                continue


                    last_snap = p.z0 # the last snapshot this clump appeared in

                    # now descend down the main branch
                    while prog != 0:

                        if prog > 0:    # if progenitor isn't jumper
                            # find progenitor's index in previous snapshot
                            p_snap_ind = desc_snap_ind + 1
                            pind = np.where(mtd.descendants[p_snap_ind]==prog)[0]

                        elif prog < 0:
                            p_snap_ind = _get_snap_ind(p, mtd.progenitor_outputnrs[desc_snap_ind][dind])
                            pind = np.where(mtd.descendants[p_snap_ind]==-prog)[0]

                        if pind.shape[0] != 1:
                            print('got pind.shape=', pind.shape)
                            print('something went wrong, exiting')
                            raise IndexError

                        # prepare for next round
                        desc_snap_ind = p_snap_ind
                        dind = pind
                        last_snap = p_snap_ind
                        prog = mtd.progenitors[p_snap_ind][pind]

                        print(
                            "----- snap", p.outputnrs[desc_snap_ind], 
                            "desc:", mtd.descendants[desc_snap_ind][dind], 
                            "npart", mtd.npart[desc_snap_ind][dind], 
                            "halo?", mtd.is_halo[desc_snap_ind][dind], 
                            )




    return







#===============================================
def get_nr_of_branches(p, r, mtd, hd):
#===============================================
    """
    Get the number of branches.
    Adapted from the old eval_tree script, since the new recursive version is too slow.
        p:      params object
        r:      results object
        mtd:    mtreedata object
        hd:     clumpdata object
    """


    print("Computing number of branches")

    # store nr of branches here by halo index 
    nr_of_branches = np.zeros(mtd.descendants[p.z0].shape, dtype=np.int)

    # store the index of the root of this tree
    # initialize to -1 for checks during treebuilding
    root = [ -np.ones(mtd.descendants[out].shape, dtype=np.int) for out in range(p.nout) ]


    # initialize root
    for i, desc in enumerate(mtd.descendants[p.z0]):
        if desc > 0:
            root[p.z0][i] = i


    # TODO: for debugging
    # debind = np.where(mtd.descendants[p.z0]==6096)[0]

    for out in range(p.z0, p.nout):
        for i, desc in enumerate(mtd.descendants[out]):

            # first propagate the root to mergers in this snapshot
            # you don't know the order a priori, so need two loops for no bugs
            if desc < 0:
                dind = np.where(mtd.descendants[out] == -desc)[0]
                root[out][i] = root[out][dind]

        for i, desc in enumerate(mtd.descendants[out]):
            prog = mtd.progenitors[out][i]

            if desc != 0:
                # propagate root
                if prog > 0:
                    # if not jumper
                    p_snap_ind = out+1
                    pind = np.where(mtd.descendants[p_snap_ind]==prog)[0]
                    root[p_snap_ind][pind] = root[out][i]
                elif prog < 0:
                    # if jumper:
                    p_snap_ind = _get_snap_ind(p, mtd.progenitor_outputnrs[out][i])
                    pind = np.where(mtd.descendants[p_snap_ind]==-prog)[0]
                    root[p_snap_ind][pind] = root[out][i]

            if desc < 0: # found merger
                # TODO: for debugging
                #  if root[out][i] == root[p.z0][debind]:
                #      print("Found merger for root", root[out][i], "root desc:", mtd.descendants[p.z0][debind])
                #      print("Desc:", desc, "prog", prog)
                if root[out][i] >= 0: # check that you're not working for a jumper at z = 0
                    nr_of_branches[root[out][i]] += 1


    # now write the results down properly
    for i, n in enumerate(mtd.npart[p.z0]):
        r.add_nr_branches(nr_of_branches[i]+1, n) # add + 1 for main branch
        # only store halos
        if mtd.is_halo[p.z0][i]:
            # TODO: mass threshold here?
            if (not p.sussing) or mtd.mass[p.z0][i] >= p.mth_main:
                r.add_nr_branches_main(nr_of_branches[i]+1, n) # add + 1 for main branch
        else:
            # TODO: mass threshold here?
            if (not p.sussing) or mtd.mass[p.z0][i] >= p.mth_sub:
                r.add_nr_branches_sub(nr_of_branches[i]+1, n) # add + 1 for main branch


    return








#=====================================================
def get_number_direct_progenitors(p, r, mtd, sd):
#=====================================================
    """
    Get the number of direct progenitors for every main halo
    i.e. how many progenitors from the past snapshot have been 
    merged into this descendant? Then add the main progenitor
    on top.

        p:      params object
        r:      results object
        mtd:    mtreedata object
        sd:     snapshotdata object
    """

    print("Computing number of direct progenitors")


    for snap in range(p.z0, p.nout):

        # only include up to redshift 2
        if sd.redshift[snap] > 2.0:
            continue

        for d, desc in enumerate(mtd.descendants[snap]):

            #  if not mtd.is_halo[snap][d]:
            #      continue

            if desc <= 0:
                continue

            #  if mtd.mass[snap][d] >= p.mth_main:

            if mtd.progenitors[snap][d] == 0:
                # this is a new clump. There can't be any mergers.
                dirprogs = 0
            else:
                # all mergers + 1 for main progenitor
                dirprogs = len(np.where(mtd.descendants[snap] == -desc)[0]) + 1 # + 1 for main progenitor

            r.add_dirprog(dirprogs)


    return








#  #=====================================================
#  def get_displacements(p, r, mtd, sd):
#  #=====================================================
#      """
#      Get the displacements for every main halo
#      progenitor-descendant pair
#
#          p:      params object
#          r:      results object
#          mtd:    mtreedata object
#          sd:     snapshotdata object
#      """
#
#      print("Computing displacements")
#
#      # TODO: something's not quite right here...
#
#
#      for snap in range(p.z0, p.nout):
#
#          for d, desc in enumerate(mtd.descendants[snap]):
#
#              if desc > 0 and mtd.is_halo[snap][d]:
#                  # skip mergers, pruned trees, and subhaloes
#
#                  if mtd.mass[snap][d] >= p.mth_main:
#                  #  if mtd.mass[snap][d] >= 1e12 / 0.703: # always use Srisawat threshold here
#                      continue
#
#                  prog = mtd.progenitors[snap][d]
#
#                  if prog <= 0:
#                      # no jumper or new descendant
#                      continue
#
#                  p_snap_ind = snap + 1
#                  pind = np.where(mtd.descendants[p_snap_ind]==prog)[0]
#
#                  if pind.shape[0] != 1:
#                      print('got pind.shape=', pind.shape)
#                      print('something went wrong, exiting')
#                      raise IndexError
#                  pind = np.asscalar(pind)
#
#                  if mtd.is_halo[p_snap_ind][pind]:
#                      if mtd.mass[p_snap_ind][pind] >= p.mth_main:
#                      #  if mtd.mass[p_snap_ind][pind] >= 1e12 / 0.703: # always use Srisawat threshold here
#
#                          delta_r = _calc_displacement(snap, d, p_snap_ind, pind, mtd, sd)
#                          r.add_halo_displacement(delta_r)
#
#      return













#=====================================================
def get_nr_of_branches_recursive(p, r, mtd, hd):
#=====================================================
    """
    Get the number of branches for all haloes at z0
    This method works, is elegant, but unfortunately too slow to
    be applicable.

        p:      params object
        r:      results object
        mtd:    mtreedata object
        hd:     clumpdata object
    """

    print("Counting number of branches")

    #-------------------------------------------------------
    def walk_tree_branches(dind, d_snap_ind, nbranches):
    #-------------------------------------------------------
        """
        Recursive function to go down the tree and count branches
        
            dind:       index of descendant in array at snapshot d_snap_ind
            d_snap_ind: snapshot index of descendant
            nbranches:  current number of branches

        returns nbranches:
            current number of branches
        """

        # go down main branch
        prog = mtd.progenitors[d_snap_ind][dind]

        if prog > 0:    # if progenitor isn't jumper
            # find progenitor's index in previous snapshot
            p_snap_ind = d_snap_ind + 1
            pind = np.where(mtd.descendants[p_snap_ind]==prog)[0]

        elif prog < 0:
            p_snap_ind = _get_snap_ind(p, mtd.progenitor_outputnrs[d_snap_ind][dind])
            pind = np.where(mtd.descendants[p_snap_ind]==-prog)[0]

        else:
            # progenitor = 0, you reached a leaf
            return nbranches


        if pind.shape[0] != 1:
            print('got pind.shape=', pind.shape)
            print("pind:", pind, "prog:", prog, "desc:", mtd.descendants[d_snap_ind][dind])
            print('something went wrong in case 1, exiting')
            raise IndexError


        nbranches = walk_tree_branches(pind, p_snap_ind, nbranches)

        # now go down branches
        desc = mtd.descendants[d_snap_ind][dind]

        for b, d in enumerate(mtd.descendants[d_snap_ind]):
            if d == - desc:
                prog = mtd.progenitors[d_snap_ind][b]

                if prog > 0:    # if progenitor isn't jumper
                    # find progenitor's index in previous snapshot
                    p_snap_ind = d_snap_ind + 1
                    pind = np.where(mtd.descendants[p_snap_ind]==prog)[0]
                elif prog < 0:
                    p_snap_ind = _get_snap_ind(p, mtd.progenitor_outputnrs[d_snap_ind][dind])
                    pind = np.where(mtd.descendants[p_snap_ind]==-prog)[0]
                else:
                    # progenitor = 0, you reached a leaf
                    # skip this loop iteration over branches
                    continue

                print("-- Adding a branch for clump", prog, "at", p.outputnrs[p_snap_ind], "nbranches=", nbranches+1)
                nbranches = walk_tree_branches(pind, p_snap_ind, nbranches+1)

        return nbranches




    #----------------------------------
    # Main Loop
    #----------------------------------

    for c, clump in enumerate(mtd.descendants[p.z0]):
    #  for c, clump in [(1, 2)]:

        if clump > 0:
            # clump = 0 is removed jumper;
            # clump < 0 is merger; we're only taking main branch here.
            # the mergers will be taken care of in the recursive function below.

            desc_is_halo = mtd.is_halo[p.z0][c]

            if desc_is_halo:
                print("-----------------------------------------------------") # in case you're printing which branch you're adding
                print("-- Root:", clump)
                nbranches = walk_tree_branches(c, p.z0, 0)
                r.add_nr_branches(nbranches, mtd.npart[p.z0][c])


    return









#=============================================
def get_mass_evolution(p, r, mtd, hd, sd, cd):
#=============================================
    """
    Compute mass fluctuations and growth along 
    the main branches of z=0 haloes.

        p:      params object
        r:      results object
        mtd:    mtreedata object
        hd:     halodata object
        sd:     snapshotdata object
        cd:     clumpdata object
    """


    print("Computing mass evolution")


    #----------------------------
    class masscompdata:
    #----------------------------
        """
        Objects to encapsulate only the necessary data
        for the tree walk to compute the mass growths and flucts
        """

        def __init__(self, ind=None, snap_ind=None):
            self.ind = ind           # index of clump
            self.snap_ind = snap_ind # index of snapshot of clump
            self.massgrowth = None   # Store mass growth to k+1, if available
            self.main = False        # computation for main haloes done?
            self.sub = False         # computation for subhaloes done?
            self.any = False         # computation for any case done?

            self.is_halo = None
            if ind is not None and snap_ind is not None:
                self.is_halo = mtd.is_halo[snap_ind][ind]
            else:
                self.is_halo = False

            if p.do_subhalo_investigation:
                self.clump_level = None
                if ind is not None and snap_ind is not None:
                    descID = mtd.descendants[snap_ind][ind]
                    desc_ind = np.where(cd.clumpid[snap_ind]==descID)[0]
                    try:
                        self.clump_level = np.asscalar(cd.clump_level[snap_ind][desc_ind])
                    except ValueError:
                        if mtd.npart[snap_ind][ind] > 20:
                            print(
                                "Error: Missing clump in clump_data with > 20 particles ID", 
                                descID, 
                                "snapshot", p.outputnrs[snap_ind], 
                                "nparts", mtd.npart[snap_ind][ind])
            return

    #-------------------------------------------
    def calc_halo_displacement(kplusone, kzero):
    #-------------------------------------------
        """
        Whether to calculate the displacement for the halo population
        """

        #  M > 10^12 M_Sol / h
        mthresh_srisawat = 1e12 / 0.703

        if kplusone.is_halo and kzero.is_halo:
            if kzero.snap_ind - kplusone.snap_ind == 1: # only do for non-jumpers!
                md = mtd.mass[kplusone.snap_ind][kplusone.ind]
                mp = mtd.mass[kzero.snap_ind][kzero.ind]
                if md >= p.mth_main and mp >= p.mth_main:
                #  if md >= mthresh_srisawat and mp >= mthresh_srisawat:
                    return True
        return False


    #-------------------------------------------
    def calc_dlogM(kplusone, kzero):
    #-------------------------------------------
        """
        Whether to calculate the log mass growth for the halo population
        """


        if (kzero.snap_ind - kplusone.snap_ind == 1): # only non-jumpers

            md = mtd.mass[kplusone.snap_ind][kplusone.ind]
            mp = mtd.mass[kzero.snap_ind][kzero.ind]

            if p.sussing:
                #  M > 10^12 M_Sol / h
                mthresh_srisawat = 1e12 / 0.703
                if (md > mthresh_srisawat) and (mp > mthresh_srisawat):
                    return True
            else:
                if mp >= p.mth_main and md >= p.mth_main:
                    return True

        return False


    #------------------------------------
    def calc_halo(kplusone, kzero):
    #------------------------------------
        """
        Whether to calculate the mass growth for the halo population
        """
        if kplusone.is_halo and kzero.is_halo:
            if kzero.snap_ind - kplusone.snap_ind == 1: # only do for non-jumpers!
                md = mtd.mass[kplusone.snap_ind][kplusone.ind]
                mp = mtd.mass[kzero.snap_ind][kzero.ind]
                if md >= p.mth_main and mp >= p.mth_main:
                    return True
        return False


    #------------------------------------
    def calc_subhalo(kplusone, kzero):
    #------------------------------------
        """
        Whether to calculate the mass growth for the subhalo population
        """
        if not kplusone.is_halo and not kzero.is_halo:
            if kzero.snap_ind - kplusone.snap_ind == 1: # only do for non-jumpers!
                md = mtd.mass[kplusone.snap_ind][kplusone.ind]
                mp = mtd.mass[kzero.snap_ind][kzero.ind]
                if md >= p.mth_sub and mp >= p.mth_sub:
                    return True
        return False

    #------------------------------------
    def calc_any(kplusone, kzero):
    #------------------------------------
        """
        Whether to calculate the mass growth for any population
        """
        if kzero.snap_ind - kplusone.snap_ind == 1: # only do for non-jumpers!
            md = mtd.mass[kplusone.snap_ind][kplusone.ind]
            mp = mtd.mass[kzero.snap_ind][kzero.ind]
            if mp >= p.mth_main and md >= p.mth_main:
                return True
        return False



    extreme_growth_desc_ind = []
    extreme_growth_desc_snap_ind = []
    extreme_growth_prog_ind = []
    extreme_growth_prog_snap_ind = []

    class extra_data_to_gather():
        def __init__(self):
            self.mass_and_level_decrease = 0
            self.mass_and_level_increase = 0
            self.total_subhaloes_in_growth = 0

    #-----------------------------------------
    def get_mass_growth(kplusone, kzero, extra):
    #-----------------------------------------
        """
        Compute the mass growth between k+1 and k_0, if applicable
        store results directly, and write them in the masscompdata
        objects where needed

        extra: extra data to gather
        """
        md = mtd.mass[kplusone.snap_ind][kplusone.ind]
        mp = mtd.mass[kzero.snap_ind][kzero.ind]
        td = sd.times[kplusone.snap_ind]
        tp = sd.times[kzero.snap_ind]

        if calc_any(kplusone, kzero):
            mgrowth = _calc_mass_growth(md, mp, td, tp)
            r.add_any_growth(mgrowth)
            kzero.any = True
            kzero.massgrowth = mgrowth

            if abs(mgrowth) > 0.95 and p.do_extra_mass_evolution_checks:
                #  print("--- Got growth {0:.3f} for descendant {1:d} snap {2:d} npart {3:.1f}".format(
                #          mgrowth, mtd.descendants[kplusone.snap_ind][kplusone.ind],
                #          p.outputnrs[kplusone.snap_ind], mtd.npart[kplusone.snap_ind][kplusone.ind])
                #          )
                extreme_growth_desc_ind.append(kplusone.ind)
                extreme_growth_desc_snap_ind.append(kplusone.snap_ind)
                extreme_growth_prog_ind.append(kzero.ind)
                extreme_growth_prog_snap_ind.append(kzero.snap_ind)

        if calc_halo(kplusone, kzero):
            if kzero.massgrowth is not None: # 'any' caught it
                r.add_halo_growth(kzero.massgrowth)
            else:
                mgrowth = _calc_mass_growth(md, mp, td, tp)
                r.add_halo_growth(mgrowth)
                kzero.massgrowth = mgrowth
            kzero.main = True
            #  dlogM = _calc_dlogMdlogt(md, mp, td, tp)
            #  r.add_halo_logM(dlogM)

        elif calc_subhalo(kplusone, kzero):
            if kzero.massgrowth is not None: # 'any' caught it
                r.add_subhalo_growth(kzero.massgrowth)
            else:
                mgrowth = _calc_mass_growth(md, mp, td, tp)
                r.add_subhalo_growth(mgrowth)
                kzero.massgrowth = mgrowth
            #  dlogM = _calc_dlogMdlogt(md, mp, td, tp)
            #  r.add_subhalo_logM(dlogM)

            kzero.sub = True
            
            #  if p.do_subhalo_investigation:
            if p.do_subhalo_investigation and abs(kzero.massgrowth) > 0.75:
                if kzero.clump_level is not None:
                    #  print("Subhalo Growth {0:.3f} levels {1:d} {2:d}".format(
                    #          kzero.massgrowth, kzero.clump_level, kplusone.clump_level))
                    if kzero.massgrowth > 0 and kzero.clump_level < kplusone.clump_level:
                        extra.mass_and_level_increase += 1
                    elif kzero.massgrowth < 0 and kzero.clump_level > kplusone.clump_level:
                        extra.mass_and_level_decrease += 1
                    extra.total_subhaloes_in_growth += 1
                else:
                    print("Using missing clump...")



        if calc_dlogM(kplusone, kzero):
            dlogM = _calc_dlogMdlogt(md, mp, td, tp)
            r.add_any_logM(dlogM)

        return extra




    #-------------------------------------------------
    def get_mass_fluct(kzero, kminusone):
    #-------------------------------------------------
        """
        Compute the mass growth fluctuation around k_0, if applicable
        store results directly
        """

        if kzero.massgrowth is None:
            return
        if kminusone.massgrowth is None:
            return
        
        fluct = _calc_mass_fluct(kminusone.massgrowth, kzero.massgrowth)

        if kzero.main and kminusone.main:
            r.add_halo_fluct(fluct)
        if kzero.sub and kminusone.sub:
            r.add_subhalo_fluct(fluct)
        if kzero.any and kminusone.any:
            # if main is true for one, the any flag will be set anyhow
            r.add_any_fluct(fluct)

        return


    #-----------------------------------------
    def get_displacement(kplusone, kzero):
    #-----------------------------------------
        """
        Compute the mass growth between k+1 and k_0, if applicable
        store results directly
        """

        mB = mtd.mass[kplusone.snap_ind][kplusone.ind]
        mA = mtd.mass[kzero.snap_ind][kzero.ind]
        xB = mtd.x[kplusone.snap_ind][kplusone.ind] / sd.aexp[kplusone.snap_ind]
        xA = mtd.x[kzero.snap_ind][kzero.ind] / sd.aexp[kzero.snap_ind]
        vB = mtd.v[kplusone.snap_ind][kplusone.ind] / sd.aexp[kplusone.snap_ind]
        vA = mtd.v[kzero.snap_ind][kzero.ind] / sd.aexp[kzero.snap_ind]
        tB = sd.times[kplusone.snap_ind]
        tA = sd.times[kzero.snap_ind]
        rhoCB = sd.rho_crit[kplusone.snap_ind]
        rhoCA = sd.rho_crit[kzero.snap_ind]
        r200B = compute_R200(mB, rhoCB) / sd.aexp[kplusone.snap_ind]
        r200A = compute_R200(mA, rhoCA) / sd.aexp[kzero.snap_ind]
        unit_lB = sd.unit_l[kplusone.snap_ind] / sd.aexp[kplusone.snap_ind]

        # Follow only  main branch
        if calc_halo_displacement(kplusone, kzero):
            delta_r = _calc_displacement(mA, mB, xA, xB, vA, vB, tA, tB, r200A, r200B, unit_lB)
            r.add_halo_displacement(delta_r)

        return







    #-----------------------------------------------------
    def walk_tree_main_branch(kplusone, kzero, extra):
    #-----------------------------------------------------
        """
        Walk the tree down the main branch only, compute the
        mass growths and mass growth fluctuations

            kplusone:   k+1: descendant of kzero
            kzero:      k_0: clump under investigation
            extra:      extra data to gather

        kplusone is needed for computation over non-adjacent snapshots
        """

        # first fill in the missing data in kminusone

        dind = kzero.ind
        d_snap_ind = kzero.snap_ind

        prog = mtd.progenitors[d_snap_ind][dind]

        if prog > 0:    # if progenitor isn't jumper
            # find progenitor's index in previous snapshot
            p_snap_ind = d_snap_ind + 1
            pind = np.where(mtd.descendants[p_snap_ind]==prog)[0]
        elif prog < 0:
            p_snap_ind = _get_snap_ind(p, mtd.progenitor_outputnrs[d_snap_ind][dind])
            pind = np.where(mtd.descendants[p_snap_ind]==-prog)[0]
        else:
            # you've reached the leaf
            return

        pind = np.asscalar(pind)


        kminusone = masscompdata(pind, p_snap_ind)


        # compute mass growth between k0 and k-1 if applicable
        extra = get_mass_growth(kzero, kminusone, extra)

        # compute mass growth fluctuation of k0 and if applicable
        get_mass_fluct(kzero, kminusone)

        # compute mass growth fluctuation of k0 and if applicable
        get_displacement(kzero, kminusone)

        # recurse
        extra = walk_tree_main_branch(kzero, kminusone, extra)
        return extra




    #--------------------
    # Main Loop
    #--------------------

    extra = extra_data_to_gather()

    for c, clump in enumerate(mtd.descendants[p.z0]):

        if clump > 0:
            # clump = 0 is removed jumper;
            # clump < 0 is merger; we're only taking main branch here

            kzero = masscompdata(c, p.z0) # initialize k_zero
            walk_tree_main_branch(masscompdata(), kzero, extra)

    

    # This part is to identify stong fluctuations
    # between different simulations.
    # Write stdout to a file, and keep only this printed
    # output. Then use find_clumps_not_in_file.py to
    # extract descendants and progenitors, which you then
    # should copypaste into run_extreme_growth_plots.sh
    # file to run the comparisons

    if p.do_extra_mass_evolution_checks:
        extreme_growth_desc_ind = np.array(extreme_growth_desc_ind)
        extreme_growth_prog_ind = np.array(extreme_growth_prog_ind)
        extreme_growth_desc_snap_ind = np.array(extreme_growth_desc_snap_ind)
        extreme_growth_prog_snap_ind = np.array(extreme_growth_prog_snap_ind)
        snapind_min = extreme_growth_desc_snap_ind.min()
        snapind_max = extreme_growth_desc_snap_ind.max()

        for s in range(snapind_min, snapind_max + 1):
            mask = extreme_growth_desc_snap_ind == s
            dinds = extreme_growth_desc_ind[mask]
            pinds = extreme_growth_prog_ind[mask]
            dsnapinds = extreme_growth_desc_snap_ind[mask]
            psnapinds = extreme_growth_prog_snap_ind[mask]
            dnparts = mtd.npart[s][dinds]
            sortinds = np.argsort(dnparts)[::-1]
            for i in range(sortinds.shape[0]):
                ind = sortinds[i]
                print(
                        mtd.descendants[s][dinds[ind]], 
                        p.outputnrs[s], 
                        mtd.mass[s][dinds[ind]], 
                        mtd.x[s][dinds[ind]][0], 
                        mtd.x[s][dinds[ind]][1], 
                        mtd.x[s][dinds[ind]][2], 
                        mtd.descendants[psnapinds[ind]][pinds[ind]], 
                        p.outputnrs[psnapinds[ind]], 
                        mtd.mass[psnapinds[ind]][pinds[ind]], 
                        mtd.x[psnapinds[ind]][pinds[ind]][0], 
                        mtd.x[psnapinds[ind]][pinds[ind]][1], 
                        mtd.x[psnapinds[ind]][pinds[ind]][2], 
                        dnparts[ind]
                    )



    if p.do_subhalo_investigation:
        print(
                "Total subhaloes:", extra.total_subhaloes_in_growth, 
                "Mass+level inc:", extra.mass_and_level_increase, 
                "Mass+level dec:", extra.mass_and_level_decrease)
    return







#=========================================================
def _calc_mass_fluct(current, previous):
#=========================================================
    """
    caclucate mass growth fluctuation from current
    mass growth in the loop (i.e. at current snapshot time)
    and previous loop (i.e. at later snapshot time)
    """

    return 0.5*(previous - current)



#==========================================================
def _calc_mass_growth(md, mp, td, tp):
#==========================================================
    """
    Calculate the logarithmic mass growth
    md, td: descendant mass/time
    mp, tp: progenitor mass/time
    """

    mass_growth = 2.0/np.pi*np.arctan((md - mp)*(tp + td)/(mp + md)/(td-tp))
    #  mass_growth = (md - mp)*(tp + td)/(mp + md)/(td-tp)
    #  mass_growth = (md - mp)/(td-tp)

    return mass_growth

#==========================================================
def _calc_dlogMdlogt(md, mp, td, tp):
#==========================================================
    """
    Calculate the logarithmic mass growth dlog M/dlog t
    md, td: descendant mass/time
    mp, tp: progenitor mass/time
    """

    mass_growth = (md - mp)*(tp + td)/(mp + md)/(td-tp)

    return mass_growth




#==================================================================================
def _calc_displacement(mA, mB, xA, xB, vA, vB, tA, tB, r200A, r200B, unit_lB):
#==================================================================================
    """
    calculate the displacement for halo A and B 

    mA, mB: masses
    xA, xB: locations
    vA, vB: velocities
    tA, tB: cosmic times
    r200A, r200B: radii of 200 critical densities containing mA, mB
    unit_lB: unit length at snapshot B
    """

    dt = tB - tA
    dx = xB - xA
    for i in range(3):
        # ASSUMING BOXSIZE = 1 IN INTERNAL UNITS
        if dx[i] > 0.5 * unit_lB:
            #  print("Correcting periodicity", dx[i], dx[i] - unit_lB, unit_lB)
            dx[i] -= unit_lB
        if dx[i] < -0.5 * unit_lB:
            #  print("Correcting periodicity", dx[i], dx[i] + unit_lB, unit_lB)
            dx[i] += unit_lB

    vector_quantity = dx - 0.5 * (vB + vA) * dt
    abs_quantity = vector_quantity[0]**2 + vector_quantity[1]**2 + vector_quantity[2]**2
    abs_quantity = np.sqrt(abs_quantity)


    vec_vel = (vA + vB)*dt
    vec_abs = vec_vel[0]**2 + vec_vel[1]**2 + vec_vel[2]**2
    vec_abs = np.sqrt(vec_abs)

    denominator = 0.5 * (r200A + r200B + vec_abs)

    delta = abs_quantity / denominator

    return delta


    

#============================
def _get_snap_ind(p, snap):
#============================
    """
    Computes the snapshot index in mtreedata/halodata/snapshotdata
    arrays for a given snapshot number snap

    p:      params object
    snap:   snapshot number (e.g. 70)
    """
    return np.asscalar(p.noutput - snap)

