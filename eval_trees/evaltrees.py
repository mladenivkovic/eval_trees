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
# default value, you can ignore it
#  mthresh_main = 2e11
#  mthresh_sub  = 2e11
mthresh_main = 0.
mthresh_sub  = 0.









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
    
    tree.clean_jumpers(p, r, mtd, sd)
    tree.determine_if_halo(p, mtd, hd)
    #  tree.get_mass_thresholds(p, mtd, hd)
    #  tree.count_pruned_trees(p, r, mtd, sd)

    # get geometry
    tree.get_main_branch_lengths(p, r, mtd, hd, sd)
    tree.get_nr_of_branches(p, r, mtd, hd)

    # get mass evolution
    tree.get_mass_evolution(p, r, mtd, hd, sd)

    # quick plot for checking
    #  vs.plot_mass_growth(r)
    #  vs.plot_geometry(p, r)

    rw.write_results(p, sd, r)

    print('Finished.')
