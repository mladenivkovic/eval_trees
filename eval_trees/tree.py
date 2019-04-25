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





#=========================================
def clean_jumpers(p, r, mtd, sd):
#=========================================
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
    
    is_halo = np.empty(ids.shape, dtype=np.bool)
    for i, id in enumerate(ids):
        is_halo[i] = id in hids

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

    return
    





#=========================================
def get_mass_flucts(p, r, mtd, cd, sd):
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

            while prog > 0:

                desc_is_halo = desc in cd.haloid[desc_snap_ind]

                if debug: 
                    if desc_snap_ind == p.z0:
                        print("--------------------------------------")
                    print("Clump", desc, "is halo?", desc_is_halo )

                if prog > 0:    # if progenitor isn't jumper
                    # find progenitor's index in previous snapshot
                    p_snap_ind = desc_snap_ind + 1
                    pind = np.where(mtd.descendants[p_snap_ind]==prog)[0]

                    if pind.shape[0] != 1:
                        print('got pind.shape=', pind.shape)
                        print('something went wrong, exiting')
                        raise IndexError


                    prog_is_halo = prog in cd.haloid[p_snap_ind]
                    if debug:
                        print("Prog", prog, "is halo?", prog_is_halo)

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
    return p.noutput - snap
 

