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

    pass


    

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
 

