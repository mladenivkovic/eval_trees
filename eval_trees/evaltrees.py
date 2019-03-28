#!/usr/bin/env python3

#======================================================
# Evaluate trees for paper.
# New version. Cleaner version.
# Usage: ./evaltrees.py
# use in directory where output_XXXXX dirs are
#======================================================


import numpy as np
import readwrite as rw
import cosmo
import tree
import classes as cl
import visualize as vs


#====================================
# MANUALLY SET PARAMETERS
#====================================

# mass threshold for halos and subhalos
mthresh_main = 2e11
mthresh_sub  = 2e11









#===================================
if __name__ == '__main__' :
#===================================

    p = cl.params(mthresh_main, mthresh_sub)
    c = cl.constants()
    r = cl.results()

    # Read available output, get global parameters
    rw.get_output_info(p, c)
    sd, mtd = cl.init_lists(p)
    hd = cl.halodata(p)

    rw.read_mergertree_data(p, sd, mtd, c)
    rw.read_halo_data(p, hd)

    cosmo.compute_cosmo_quantities(p, c, sd) # in particular get p.z0 and rho_crit


    #  tree.count_pruned_trees(p, r, mtd, sd)
    tree.clean_jumpers(p, r, mtd, sd)
    tree.get_mass_flucts(p, r, mtd, hd, sd)

    vs.plot_mass_growth(r)
