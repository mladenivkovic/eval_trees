#!/usr/bin/env python3


import numpy as np



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

    pass






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
                    mtd.descendants[snapind][jumpind] = 0
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

        mtd.is_halo[i] = is_halo


    return





#====================================
def get_mass_thresholds(p, mtd, hd):
#====================================
    """
    Compute mass thresholds:
    include 1000 main haloes and 200 subhaloes at z=0

        p:      params object
        mtd:    mergertree data
        hd:     halodata
    """

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

        if clump > 0:
            # clump = 0 is removed jumper;
            # clump < 0 is merger; we're only taking main branch here

            dind = c
            desc_snap_ind = p.z0
            desc = mtd.descendants[p.z0][dind]
            prog = mtd.progenitors[p.z0][dind]
            desc_is_halo = mtd.is_halo[p.z0][dind]


            # skip if threshold is not satisfied
            if desc_is_halo:
                # apply mass threshold:
                #  if mtd.mass[p.z0][dind] < p.mth_main:
                #      continue

                # don't apply mass threshold:
                pass
            else:
                # only do it for main haloes
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
            r.add_branch_length(last_snap - p.z0, mtd.npart[p.z0][c])


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


    debind = np.where(mtd.descendants[p.z0]==6100)[0]

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
                if root[out][i] == root[p.z0][debind]:
                    print("Found merger for root", root[out][i], "root desc:", mtd.descendants[p.z0][debind])
                    print("Desc:", desc, "prog", prog)
                if root[out][i] >= 0: # check that you're not working for a jumper at z = 0
                    nr_of_branches[root[out][i]] += 1


    # now write the results down properly
    for i, n in enumerate(mtd.npart[p.z0]):
        # only store halos
        if mtd.is_halo[p.z0][i]:
            r.add_nr_branches(nr_of_branches[i]+1, n) # add + 1 for main branch


    return
















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









#=========================================
def get_mass_evolution(p, r, mtd, cd, sd):
#=========================================
    """
    Compute mass fluctuations and growth along 
    the main branches of z=0 haloes.

        p:      params object
        r:      results object
        mtd:    mtreedata object
        cd:     clumpdata object
        sd:     snapshotdata object
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
            return



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
            if md >= p.mth_sub and mp >= p.mth_sub:
                return True
        return False




    #-----------------------------------------
    def get_mass_growth(kplusone, kzero):
    #-----------------------------------------
        """
        Compute the mass growth between k+1 and k_0, if applicable
        store results directly, and write them in the masscompdata
        objects where needed
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

        if calc_halo(kplusone, kzero):
            if kzero.massgrowth is not None: # 'any' caught it
                r.add_halo_growth(kzero.massgrowth)
            else:
                mgrowth = _calc_mass_growth(md, mp, td, tp)
                r.add_halo_growth(mgrowth)
                kzero.massgrowth = mgrowth
            kzero.main = True

        elif calc_subhalo(kplusone, kzero):
            if kzero.massgrowth is not None: # 'any' caught it
                r.add_subhalo_growth(kzero.massgrowth)
            else:
                mgrowth = _calc_mass_growth(md, mp, td, tp)
                r.add_subhalo_growth(mgrowth)
                kzero.massgrowth = mgrowth

            kzero.sub = True

        return





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
        elif kzero.sub and kminusone.sub:
            r.add_subhalo_fluct(fluct)
        if kzero.any and kminusone.any:
            # if main is true for one, the any flag will be set anyhow
            r.add_any_fluct(fluct)

        return









    #-----------------------------------------------------
    def walk_tree_main_branch(kplusone, kzero):
    #-----------------------------------------------------
        """
        Walk the tree down the main branch only, compute the
        mass growths and mass growth fluctuations

            kplusone:   k+1: descendant of kzero
            kzero:      k_0: clump under investigation
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


        kminusone = masscompdata(pind, p_snap_ind)


        # compute mass growth between k0 and k-1 if applicable
        get_mass_growth(kzero, kminusone)

        # compute mass growth fluctuation of k0 and if applicable
        get_mass_fluct(kzero, kminusone)

        # recurse
        walk_tree_main_branch(kzero, kminusone)
        return




    #--------------------
    # Main Loop
    #--------------------

    for c, clump in enumerate(mtd.descendants[p.z0]):

        if clump > 0:
            # clump = 0 is removed jumper;
            # clump < 0 is merger; we're only taking main branch here

            is_subhalo = not mtd.is_halo[p.z0][c]
            m = mtd.mass[p.z0][c]

            if m > p.mth_main or (is_subhalo and m > p.mth_sub):
                kzero = masscompdata(c, p.z0) # initialize k_zero
                walk_tree_main_branch(masscompdata(), kzero)


    return







#===============================================
def get_mass_evolution_old(p, r, mtd, cd, sd):
#===============================================
    """
    Compute mass fluctuations and growth along 
    the main branches of z=0 haloes.

        p:      params object
        r:      results object
        mtd:    mtreedata object
        cd:     clumpdata object
        sd:     snapshotdata object
    """


    print("Computing mass evolution")




    # loop over clumps at z=0 and walk down their main branch

    debug = False
    extremecount = 0

    for c, clump in enumerate(mtd.descendants[p.z0]):
        prevhalomgrowth = None
        prevsubhalomgrowth = None
        prevbothmgrowth = None

        if clump > 0:
            # clump = 0 is removed jumper;
            # clump < 0 is merger; we're only taking main branch here

            desc_snap_ind = p.z0
            dind = c
            desc = mtd.descendants[p.z0][c]
            prog = mtd.progenitors[p.z0][c]

            loop = True

            #  while prog > 0: # this doesn't allow for jumpers
            while loop:

                desc_is_halo = desc in cd.haloid[desc_snap_ind]

                if debug: 
                    if desc_snap_ind == p.z0:
                        print("--------------------------------------")
                    print("Clump", desc, "is halo?", desc_is_halo , "snap", p.noutput-desc_snap_ind)

                if prog > 0:    # if progenitor isn't jumper
                    # find progenitor's index in previous snapshot
                    p_snap_ind = desc_snap_ind + 1
                    pind = np.where(mtd.descendants[p_snap_ind]==prog)[0]
                elif prog < 0:
                    p_snap_ind = _get_snap_ind(p, mtd.progenitor_outputnrs[desc_snap_ind][dind])
                    try:
                        pind = np.where(mtd.descendants[p_snap_ind]==-prog)[0]
                    except TypeError:
                        print(p_snap_ind)
                        print(mtd.descendants[p_snap_ind]==-prog)
                        print(pind)
                        print(type(pind))
                        quit()
                else:
                    break


                if pind.shape[0] != 1:
                    print('got pind.shape=', pind.shape)
                    print('something went wrong, exiting')
                    raise IndexError


                prog_is_halo = prog in cd.haloid[p_snap_ind]
                if debug:
                    print("Prog", prog, "is halo?", prog_is_halo , "snap", p.noutput-p_snap_ind)

                md = mtd.mass[desc_snap_ind][dind]
                mp = mtd.mass[p_snap_ind][pind]
                td = sd.times[desc_snap_ind]
                tp = sd.times[p_snap_ind]


                calc_desc_halo = desc_is_halo and (md >= p.mth_main)
                calc_desc_subhalo = (not desc_is_halo) and (md >= p.mth_sub)
                calc_desc_both = calc_desc_halo or (md >= p.mth_main) # where main and subhaloes are counted together

                calc_prog_halo = prog_is_halo and (mp >= p.mth_main)
                calc_prog_subhalo = (not prog_is_halo) and (mp >= p.mth_sub)
                calc_prog_both = calc_prog_halo or (mp >= p.mth_main)


                halomgrowth = None
                subhalomgrowth = None
                bothmgrowth = None

                if calc_desc_halo and calc_prog_halo:
                    halomgrowth = _calc_mass_growth(md, mp, td, tp)
                    r.add_halo_growth(halomgrowth)
                if calc_desc_subhalo and calc_prog_subhalo:
                    subhalomgrowth = _calc_mass_growth(md, mp, td, tp)
                    r.add_subhalo_growth(subhalomgrowth)

                if calc_desc_both and calc_prog_both:
                    bothmgrowth = _calc_mass_growth(md, mp, td, tp)
                    r.add_any_growth(bothmgrowth)
                    
                if prevhalomgrowth is not None and halomgrowth is not None: 
                    halomassfluct = _calc_mass_fluct(halomgrowth, prevhalomgrowth)
                    r.add_halo_fluct(halomassfluct)
                if prevsubhalomgrowth is not None and subhalomgrowth is not None: 
                    subhalomassfluct = _calc_mass_fluct(subhalomgrowth, prevsubhalomgrowth)
                    r.add_subhalo_fluct(subhalomassfluct)
                if prevbothmgrowth is not None and bothmgrowth is not None: 
                    bothmassfluct = _calc_mass_fluct(bothmgrowth, prevbothmgrowth)
                    r.add_any_fluct(bothmassfluct)

                
                # prepare next loop
                
                prevhalomgrowth = halomgrowth
                prevsubhalomgrowth = subhalomgrowth
                prevbothmgrowth = bothmgrowth

                desc_snap_ind = p_snap_ind
                dind = pind
                desc = prog
                prog = mtd.progenitors[p_snap_ind][pind]



                if halomgrowth is None and subhalomgrowth is None and bothmgrowth is None:
                    newgrowths = False
                else:
                    newgrowths = True
                if prevbothmgrowth is None and prevsubhalomgrowth is None and prevhalomgrowth is None:
                    prevgrowths = False
                else:
                    prevgrowths = True

                if (not newgrowths) and (not prevgrowths):
                    break


    print("Finished eval")




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

    return mass_growth

    

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
 

