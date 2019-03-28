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

    for c, clump in enumerate(mtd.descendants[p.z0]):
        if clump > 0:
            # clump = 0 is removed jumper;
            # clump < 0 is merger; we're only taking main branch here

            desc_snap_ind = p.z0
            dind = c
            prog = mtd.progenitors[p.z0][c]

            loop = True

            while loop:

                desc_is_halo = clump in cd.haloid[desc_snap_ind]

                if prog > 0:    # if progenitor isn't jumper
                    # find progenitor's index in previous snapshot
                    p_snap_ind = desc_snap_ind + 1
                    pind = np.where(mtd.descendants[p_snap_ind]==prog)[0]

                    if pind.shape[0] != 1:
                        print('got pind.shape=', pind.shape)
                        print('something went wrong, exiting')
                        raise IndexError


                    prog_is_halo = prog in cd.haloid[p_snap_ind]

                    md = mtd.mass[desc_snap_ind][dind]
                    mp = mtd.mass[p_snap_ind][pind]
                    td = sd.times[desc_snap_ind]
                    tp = sd.times[p_snap_ind]

                    calc_mgrowth = True
                    if prog_is_halo:
                        calc_mgrowth = calc_mgrowth and (mp >= p.mth_main)
                    else:
                        calc_mgrowth = calc_mgrowth and (mp >= p.mth_sub)

                    if desc_is_halo:
                        calc_mgrowth = calc_mgrowth and (md >= p.mth_main)
                    else:
                        calc_mgrowth = calc_mgrowth and (md >= p.mth_sub)

                    if calc_mgrowth:
                        mgrowth = _calc_mass_growth(md, mp, td, tp)
                        
                        if (prog_is_halo and desc_is_halo):
                            r.add_halo_growth(mgrowth)
                        elif ((not prog_is_halo) and (not desc_is_halo)):
                            r.add_subhalo_growth(mgrowth)

                        r.add_any_growth(mgrowth)
             

                
                    # prepare next loop

                    desc_snap_ind = p_snap_ind
                    dind = pind
                    prog = mtd.progenitors[p_snap_ind][pind]


                elif prog == 0:
                    loop = False

                else:
                    loop = False


                    



#==========================================================
def _calc_mass_growth(md, mp, td, tp):
#==========================================================
    """
    Calculate the logarithmic mass growth
    md, td: descendant mass/time
    mp, tp: progenitor mass/time
    """

    mass_growth = 2.0/np.pi*np.arctan((md - mp)*(tp + td)/(mp + md)/(tp-td))

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
 

