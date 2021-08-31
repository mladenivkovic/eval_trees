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




#===================================
if __name__ == '__main__' :
#===================================

    p = cl.params()
    c = cl.constants()
    r = cl.results()

    # Read available output, get global parameters
    rw.get_output_info(p, c)
    sd, mtd = cl.init_lists(p)
    hd = cl.halodata(p)
    cd = cl.clumpdata(p)

    rw.read_mergertree_data(p, sd, mtd, c)
    rw.read_halo_data(p, hd)
    if p.do_subhalo_investigation:
        rw.read_clump_data(p, cd)

    cosmo.compute_cosmo_quantities(p, c, sd) # in particular get p.z0 and rho_crit

    tree.determine_if_halo(p, mtd, hd)
    tree.clean_jumpers(p, r, mtd, sd)
    tree.get_mass_thresholds(p, mtd, hd)

    # get geometry
    #  tree.count_pruned_trees(p, r, mtd, sd)
    #  tree.count_statistics_at_z0(p, r, mtd)
    #  tree.count_total_descendants(p, mtd)
    #  tree.get_main_branch_lengths(p, r, mtd, hd, sd)
    #  tree.get_nr_of_branches(p, r, mtd, hd)
    #  tree.get_number_direct_progenitors(p, r, mtd, sd)

    # get mass evolution
    tree.get_mass_evolution(p, r, mtd, hd, sd, cd)

    # quick plot for checking
    #  vs.plot_mass_growth(r)
    #  vs.plot_geometry(p, r)
    #  vs.plot_displacements(r)

    #  rw.write_results(p, sd, r)

    print('Finished.')
