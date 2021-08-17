#!/usr/bin/env python3


#===================================================
# This script plots the results of eval_tree.py
#===================================================

import numpy as np
import matplotlib as mpl
#  mpl.use('Agg')
import matplotlib.pyplot as plt
import pickle

# py files in this directory
from tools_for_plotting import file_selection, print_evaltrees_parameters

#-----------------------------------------
# set which dataset you want to plot
#-----------------------------------------
from tools_for_plotting import plot_selection # set it in tools_for_plotting.py !!


# set which plots to make
plot_mass_growth_halo_subhalo = False
plot_mass_fluctuations_halo_subhalo = False
plot_mass_growth_and_fluctuations_vertical = False
plot_geometry_partbins = False
plot_main_brench_length_halo_subhalo_no_partbins = False
plot_main_brench_length_halo_subhalo_partbins = False
plot_nbranches_halo_subhalo_no_partbins = False
plot_nbranches_all_no_partbins = False
plot_nbranches_halo_subhalo_particle_bins = False
plot_displacements = True
plot_dlogMdlogt = False

print_tree_statistics_table = False




params = {
    'axes.labelsize': 18,
    'font.size': 16,
    'font.family': 'serif',
    'font.serif': 'Computer Modern',
    'legend.fontsize': 16,
    'text.usetex': True,
    'lines.linewidth': 3,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
}

mpl.rcParams.update(params)


hist_bins = 100
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#e377c2', '#bcbd22', '#17becf']
alpha = 0.6

# how many particle bins do we have for tree geometry?
binnames = [r'$< 100$ particles', r'$100-500$ particles', r'$500-1000$ particles', r'$> 1000$ particles']
npartbins = len(binnames)





#===================
def main():
#===================
    
    allfiles, labelnames, suffix, linestyle = file_selection(plot_selection)
    print_evaltrees_parameters(allfiles[0])

    # columns for legends
    if len(allfiles) == 4:
        ncols = 2
    else:
        ncols = 3


    #==========================
    # Set up
    #==========================

    # mass growth with halo-subhalo divide
    if plot_mass_growth_halo_subhalo:
        # mass_growth-halo-subhalo.png
        fig1 = plt.figure(1, figsize=(18, 6))
        ax1 = fig1.add_subplot(131)
        ax2 = fig1.add_subplot(132)
        ax3 = fig1.add_subplot(133)

    # mass fluctuations with halo-subhalo divide
    if plot_mass_fluctuations_halo_subhalo:
        # mass_fluctuations-halo-subhalo.png
        fig2 = plt.figure(2, figsize=(18, 6))
        ax4 = fig2.add_subplot(131)
        ax5 = fig2.add_subplot(132)
        ax6 = fig2.add_subplot(133)

    # tree geometry
    if plot_geometry_partbins:
        # tree_geometry.png
        fig3 = plt.figure(3, figsize=(14,10)) # for 4 plots
        for b in range(2*npartbins):
            fig3.add_subplot(npartbins, 2, b+1)

    # tree geometry; halo-subhalo divide main branch lengths
    if plot_main_brench_length_halo_subhalo_partbins:
        # main-branch-lengths-halo-subhalo-partbins.png
        fig4 = plt.figure(4, figsize=(16,9))
        for b in range(2*npartbins):
            fig4.add_subplot(npartbins, 2, b+1)

    # tree geometry; halo-subhalo divide number of branches
    if plot_nbranches_halo_subhalo_particle_bins:
        # number-of-branches-halo-subhalo-partbins.png
        fig5 = plt.figure(5, figsize=(15,8))
        for b in range(2*npartbins):
            fig5.add_subplot(npartbins, 2, b+1)

    # mass growth and fluctuations, no divide
    if plot_mass_growth_and_fluctuations_vertical:
        # mass-statistics.png
        fig6 = plt.figure(6, figsize=(6, 12))
        ax7 = fig6.add_subplot(2, 1, 1)
        ax8 = fig6.add_subplot(2, 1, 2)

    # main brench length for all particle bins
    if plot_main_brench_length_halo_subhalo_no_partbins:
        # plot main branch length with halo/subhalo divide, without particle bins
        # main-branch-lengths-no-partbins.png
        fig7 = plt.figure(7, figsize=(6, 12))
        ax9 = fig7.add_subplot(2, 1, 1)
        ax10 = fig7.add_subplot(2, 1, 2)

    # number of branches for all particle bins
    if plot_nbranches_halo_subhalo_no_partbins:
        # plot number of branches with halo/subhalo divide, without particle bins
        # number-of-branches-halo-subhalo-no-partbins.png
        fig8 = plt.figure(8, figsize=(6, 12))
        ax11 = fig8.add_subplot(2, 1, 1)
        ax12 = fig8.add_subplot(2, 1, 2)

    # plot number of branches without halo/subhalo divide, without particle bins
    if plot_nbranches_all_no_partbins:
        # nbranches-all-no-partbins.png
        fig9 = plt.figure(9, figsize=(6, 6))
        ax13 = fig9.add_subplot(111)

    # displacement statistic
    if plot_displacements:
        # displacement.png
        fig10 = plt.figure(10, figsize=(6, 6))
        ax14 = fig10.add_subplot(111)

    # d logM/ dlog t
    if plot_dlogMdlogt:
        # dlogMdlogt.png
        fig11 = plt.figure(11, figsize=(7,7))
        ax15 = fig11.add_subplot(111)



    # tweaking stuff to store
    mblmin = 1000
    nbranchmin = 1000
    fluctmin = 10000
    growthmin = 10000
    displacements_max = 0
    displacements_count_max = 0

    #  nbins_nbranches = 100 # how many bins to use for number of branches; = max number of branches expected
    nbins_nbranches = 1200 # how many bins to use for number of branches; = max number of branches expected

    #=============================
    # Read in data
    #=============================

    for f,srcfname in enumerate(allfiles):

        srcfile = open(srcfname, 'rb')
        p, sd, r = pickle.load(srcfile) 
        srcfile.close()


        if plot_mass_growth_halo_subhalo:
            #------------------------------------------------
            # Plot Mass Growth with halo/subhalo separation
            #------------------------------------------------

            # Logarithmic Mass Growth
            #------------------------
 
            for arr, n, ax, norm in [(r.mg, r.mg_free, ax1, r.mg_free), (r.hmg, r.hmg_free, ax2, r.mg_free), (r.shmg, r.shmg_free, ax3, r.mg_free)] :
                cut = arr[:n]
                hist, bin_edges = np.histogram(cut, bins=hist_bins, range=(-1, 1))
                hist = hist/norm # normalize histogram
                growthmin = min(growthmin, hist[hist>0].min())
                bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])

                ax.semilogy(bin_centers, hist,
                    label = labelnames[f],
                    c=colors[f],
                    ls = linestyle[f],
                    alpha = alpha,
                    )


            # Mass Fluctuations
            #------------------
 
            for arr, n, ax, norm in [(r.mf, r.mf_free, ax4, r.mf_free), (r.hmf, r.hmf_free, ax5, r.mf_free), (r.shmf, r.shmf_free, ax6, r.mf_free)] :
                cut = arr[:n]
                hist, bin_edges = np.histogram(cut, bins=hist_bins, range=(-1, 1))
                hist = hist/norm # normalize histogram
                fluctmin = min(fluctmin, hist[hist>0].min())
                bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])

                ax.semilogy(bin_centers, hist,
                    label = labelnames[f],
                    c=colors[f],
                    ls = linestyle[f],
                    alpha = alpha,
                    )


        if plot_mass_growth_and_fluctuations_vertical:
            #------------------------------------------------------------
            # Mass growth and fluctuations without halo/subhalo divisions
            #------------------------------------------------------------

            for arr, n, ax, norm in [(r.mg, r.mg_free, ax7, r.mg_free+1)] :
                cut = arr[:n]
                hist, bin_edges = np.histogram(cut, bins=hist_bins, range=(-1, 1))
                hist += 1
                #  hist = hist/norm # normalize histogram
                growthmin = min(growthmin, hist[hist>0].min())
                bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])

                ax.semilogy(bin_centers, hist,
                    label = labelnames[f],
                    c=colors[f],
                    ls = linestyle[f],
                    alpha = alpha,
                    )


            for arr, n, ax, norm in [(r.mf, r.mf_free, ax8, r.mf_free+1)] :
                cut = arr[:n]
                hist, bin_edges = np.histogram(cut, bins=hist_bins, range=(-1, 1))
                hist += 1
                #  hist = hist/norm # normalize histogram
                fluctmin = min(fluctmin, hist[hist>0].min())
                bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])

                ax.semilogy(bin_centers, hist,
                    label = labelnames[f],
                    c=colors[f],
                    ls = linestyle[f],
                    alpha = alpha,
                    )



        if plot_dlogMdlogt:
            # ---------------------------------------
            # d log M / d log t
            # ---------------------------------------

            cut = r.logM[:r.logM_free]

            edge = 0.01
            bins1 = np.logspace(np.log10(edge), np.log10(15), 200)
            bins1 = np.flip(bins1[:-1]) * -1
            bins2 = np.logspace(np.log10(edge), np.log10(15), 200, endpoint=False)
            bins_mid = np.linspace(-edge, edge, 100, endpoint=False)

            bins = np.concatenate((bins1, bins_mid, bins2))
            hist, bin_edges = np.histogram(cut, bins=bins)
            hist += 1
            bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])

            #  print("alpha m min max", cut.min(), cut.max())

            ax15.semilogy(bin_centers, hist, 
                    label = labelnames[f],
                    c=colors[f],
                    ls = linestyle[f],
                    alpha = alpha,
                    )



        if plot_geometry_partbins:

            # Plot tree geometry: Main branch lengths, number of branches
            #--------------------------------------------------------------

            # figure subplot setup:
            #   1   2
            #   3   4
            #   5   6
            #   7   8
            # -> nr of branches: 2*(b+1)    -1 for list index of fig3.axes
            # -> MBL: 2*b + 1               -1 for list index of fig3.axes


            # plot branch lengths
            #--------------------------

            for b in range(npartbins):
                ax = fig3.axes[b*2]
                norm = r.branchlen_free[b]

                if norm == 0:
                    print("Found no entries for branch lengths in bin", binnames[b], "for", labelnames[f])
                    continue

                cut = r.branchlengths[b][:r.branchlen_free[b]]
                hist, bin_edges = np.histogram(cut, bins = p.nout, range=(1, p.nout))
                hist += 1 # plot N + 1

                #  hist = hist/norm
                mblmin = min(mblmin, hist[hist>0].min())
                #  bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
                bins_left = bin_edges[:-1]

                ax.semilogy(bins_left, hist,
                    label = labelnames[f],
                    c=colors[f],
                    ls = linestyle[f],
                    alpha = alpha,
                    )



            # plot number of branches
            #--------------------------
            for b in range(npartbins):
                ax = fig3.axes[b*2+1]
                norm = r.nbranch_free[b]

                if norm == 0:
                    print("Found no entries for number of branches in bin", binnames[b], "for", labelnames[f])
                    continue

                cut = r.nbranches[b][:r.nbranch_free[b]]
                hist, bin_edges = np.histogram(cut, bins = nbins_nbranches, range=(1, nbins_nbranches))
                hist += 1 # plot N + 1
                #  hist = hist/norm
                nbranchmin = min(nbranchmin, hist[hist>0].min())
                #  bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
                bins_left = bin_edges[:-1]

                ax.loglog(bins_left, hist,
                    label = labelnames[f],
                    c=colors[f],
                    ls = linestyle[f],
                    alpha = alpha,
                    )




        if plot_main_brench_length_halo_subhalo_partbins:
            #-----------------------------------------------------
            # Plot main branch lengths with halo/subhalo divide;
            # separate clumps by particle number bins
            #-----------------------------------------------------

            # plot branch lengths for main haloes first
            #------------------------------------------
            for b in range(npartbins):
                ax = fig4.axes[b*2]
                norm = r.branchlen_free_main[b]
                if norm == 0:
                    print("main haloes branch lengths got norm=0, bin is", binnames[b])
                    continue

                cut = r.branchlengths_main[b][:r.branchlen_free_main[b]]
                hist, bin_edges = np.histogram(cut, bins = p.nout, range=(1, p.nout))
                mblmin = min(mblmin, hist[hist>0].min())
                bins_left = bin_edges[:-1]

                ax.semilogy(bins_left, hist,
                    label = labelnames[f],
                    c=colors[f],
                    ls = linestyle[f],
                    alpha = alpha,
                    )


            # plot branch lengths for sub haloes now
            #---------------------------------------
            for b in range(npartbins):
                ax = fig4.axes[b*2+1]
                norm = r.branchlen_free_sub[b]
                if norm == 0:
                    continue

                cut = r.branchlengths_sub[b][:r.branchlen_free_sub[b]]
                hist, bin_edges = np.histogram(cut, bins = p.nout, range=(1, p.nout))
                mblmin = min(mblmin, hist[hist>0].min())
                bins_left = bin_edges[:-1]

                ax.semilogy(bins_left, hist,
                    label = labelnames[f],
                    c=colors[f],
                    ls = linestyle[f],
                    alpha = alpha,
                    )



        if plot_main_brench_length_halo_subhalo_no_partbins:
            #-----------------------------------------------------
            # Plot main branch lengths with halo/subhalo divide;
            # no particle number bins
            #-----------------------------------------------------

            # Main Haloes first
            #------------------

            hbranchlen_all = np.ones(p.nout) # plot N + 1

            for b in range(npartbins):
                norm = r.branchlen_free_main[b]
                if norm == 0:
                    print("main haloes branch lengths got norm=0, bin is", binnames[b])
                    continue

                cut = r.branchlengths_main[b][:r.branchlen_free_main[b]]
                hist, bin_edges = np.histogram(cut, bins = p.nout, range=(1, p.nout))
                hbranchlen_all += hist
                mblmin = min(mblmin, hist[hist>0].min())
                bins_left = bin_edges[:-1]

            ax9.semilogy(bins_left, hbranchlen_all,
                            label = labelnames[f],
                            c=colors[f],
                            ls = linestyle[f],
                            alpha = alpha,
                        )


            # Subhaloes now
            #------------------

            sbranchlen_all = np.ones(p.nout) # plot N + 1

            for b in range(npartbins):
                norm = r.branchlen_free_sub[b]
                if norm == 0:
                    continue

                cut = r.branchlengths_sub[b][:r.branchlen_free_sub[b]]
                hist, bin_edges = np.histogram(cut, bins = p.nout, range=(1, p.nout))
                sbranchlen_all += hist
                mblmin = min(mblmin, hist[hist>0].min())
                bins_left = bin_edges[:-1]

            ax10.semilogy(bins_left, sbranchlen_all,
                label = labelnames[f],
                c=colors[f],
                ls = linestyle[f],
                alpha = alpha,
                )


        if plot_nbranches_halo_subhalo_particle_bins:
            #-------------------------------------------------
            # Plot number of branches for haloes and subhaloes
            # separately, and use particle number bins
            #-------------------------------------------------

            # plot number of branches for main haloes first
            #----------------------------------------------
            for b in range(npartbins):
                ax = fig5.axes[b*2]
                norm = r.nbranch_free_main[b]
                if norm == 0:
                    continue

                cut = r.nbranches_main[b][:r.nbranch_free_main[b]]
                hist, bin_edges = np.histogram(cut, bins = nbins_nbranches, range=(1, nbins_nbranches))
                nbranchmin = min(nbranchmin, hist[hist>0].min())
                bins_left = bin_edges[:-1]

                ax.loglog(bins_left, hist,
                    label = labelnames[f],
                    c=colors[f],
                    ls = linestyle[f],
                    alpha = alpha,
                    )

            # plot number of branches for sub haloes now
            #-------------------------------------------
            for b in range(npartbins):
                ax = fig5.axes[b*2+1]
                norm = r.nbranch_free_sub[b]
                if norm == 0:
                    continue

                cut = r.nbranches_sub[b][:r.nbranch_free_sub[b]]
                hist, bin_edges = np.histogram(cut, bins = nbins_nbranches, range=(1, nbins_nbranches))
                nbranchmin = min(nbranchmin, hist[hist>0].min())
                bins_left = bin_edges[:-1]

                ax.loglog(bins_left, hist,
                    label = labelnames[f],
                    c=colors[f],
                    ls = linestyle[f],
                    alpha = alpha,
                    )




        if plot_nbranches_halo_subhalo_no_partbins:
            #----------------------------------------------------------
            # Plot number of branches with halo subhalo separation,
            # no particle bins
            #----------------------------------------------------------

            # First do subhaloes
            snbranches_all = np.ones(nbins_nbranches) # plot N + 1
     
            for b in range(npartbins):
                norm = r.nbranch_free_sub[b]
                if norm == 0:
                    continue

            cut = r.nbranches_sub[b][:r.nbranch_free_sub[b]]
            hist, bin_edges = np.histogram(cut, bins = nbins_nbranches, range=(1, nbins_nbranches))
            snbranches_all += hist
            nbranchmin = min(nbranchmin, hist[hist>0].min())
            bins_left = bin_edges[:-1]

            ax12.loglog(bins_left, snbranches_all,
                label = labelnames[f],
                c=colors[f],
                ls = linestyle[f],
                alpha = alpha,
                )


            # Now do main haloes
            hnbranches_all = np.ones(nbins_nbranches) # plot N + 1

            for b in range(npartbins):
                norm = r.nbranch_free_main[b]
                if norm == 0:
                    continue

            cut = r.nbranches_main[b][:r.nbranch_free_main[b]]
            hist, bin_edges = np.histogram(cut, bins = nbins_nbranches, range=(1, nbins_nbranches))
            hnbranches_all += hist
            nbranchmin = min(nbranchmin, hist[hist>0].min())
            bins_left = bin_edges[:-1]

            ax11.loglog(bins_left, hnbranches_all,
                label = labelnames[f],
                c=colors[f],
                ls = linestyle[f],
                alpha = alpha,
                )






        if plot_nbranches_all_no_partbins:
            #------------------------------------------------------------
            # Number of Branches of all clumps, no particle bins
            #------------------------------------------------------------

            #  plot number of branches without halo/subhalo divide, 
            #  without particle bins

            norm = r.dirprog_free
            cut = r.dirprogs[:r.dirprog_free]
            hist, edges = np.histogram(cut, bins=int(cut.max()))
            bins_left = edges[:-1]
            ax13.semilogy(bins_left, hist/norm,
                label = labelnames[f],
                c=colors[f],
                ls = linestyle[f],
                alpha = alpha,
                )




        if plot_displacements:
            # -------------------------------
            # Plot displacement statistics
            # -------------------------------

            cut = r.displacements[:r.displacement_free]
            displacements_max = max(displacements_max, cut.max())
            displacements_max = min(displacements_max, 10)
            hist, bin_edges = np.histogram(cut, bins=hist_bins, range=(1e-6, 10))

            #  find 90 and 99 percentile bins
            displacements_ninety_percent = None
            displacements_ninety_nine_percent = None

            cumulation = 0.
            hist_tot = hist.sum()
            for i in range(hist_bins):
                cumulation += hist[i]
                if cumulation >= 0.9 * hist_tot:
                    if displacements_ninety_percent is None:
                        displacements_ninety_percent = bin_edges[i]
                if cumulation >= 0.99 * hist_tot:
                    if displacements_ninety_nine_percent is None:
                        displacements_ninety_nine_percent = bin_edges[i]
                        break


            hist += 1 # add 1 so we can do logarithmic plots
            displacements_count_max = max(displacements_count_max, hist.max())
            bins_left = bin_edges[:-1]

            ax14.loglog(bins_left, hist,
                label = labelnames[f],
                c=colors[f],
                ls = linestyle[f],
                alpha = alpha,
                )



        if print_tree_statistics_table:
            #----------------------------------------------
            # Print the table containing average values
            # for the tree
            #----------------------------------------------

            print("Tree Geometry Table for", srcfname)
            print("-------------------------------------------------")

            print("Total clumps at z=0:", r.clumps_at_z0)
            print("Median particles in clumps at z=0:", r.median_clump_particlecount_at_z0)
            print("-------------------------------------------------")
            print("Average branch lengths")

            for b in range(npartbins):
                cut = r.branchlengths[b][:r.branchlen_free[b]]
                av = np.average(cut)
                print("{0:25} {1:.1f}".format(binnames[b], av))

            print("-------------------------------------------------")
            print("Average number of branches")

            for b in range(npartbins):
                cut = r.nbranches[b][:r.nbranch_free[b]]
                av = np.average(cut)
                print("{0:25} {1:.1f}".format(binnames[b], av))
            print("-------------------------------------------------")
            print("-------------------------------------------------")









    #==================================
    # Tweak plots and save figures
    #==================================

    plotnames = ['A: All Haloes', 'B: Main Haloes', 'C: Subhaloes']

    #-------------------------
    # Mass growth
    #-------------------------

    if plot_mass_growth_halo_subhalo:
        for i, ax in enumerate(fig1.axes):
            # for inclusive/exclusive
            ax.set_ylim(growthmin, 1e-1)
            ax.set_xlim(-1.05, 1.05)
            ax.set_xlabel(r'$\beta_M = \frac{2}{\pi}\arctan\frac{d \log M}{d \log t}$')
            ax.grid()
            if i == 0:
                ax.legend(loc='lower center', ncol=1,)
                #  ax.legend(loc='lower center', ncol=ncols,)
            ax.set_title(plotnames[i])

        fig1.axes[0].set_ylabel(r"$N/N_{tot}$")

        # in case you want a title
        fig1.suptitle("Logarithmic Mass Growth", fontsize=21)
        fig1.tight_layout(pad=0.5, w_pad=1.5, h_pad=1.0, rect=(0,0,1,0.9))

        # in case you don't want a title
        #  fig1.tight_layout(pad=0.5, w_pad=1.5, h_pad=1.0)

        figname="mass_growth-"+suffix+".png"
        fig1.savefig(figname, dpi=300)



    #-------------------------
    # Mass growth flucts
    #-------------------------

    if plot_mass_fluctuations_halo_subhalo:
        for i, ax in enumerate(fig2.axes):
            # for inclusive/exclusive
            ax.set_ylim(fluctmin, 1e-1)
            ax.set_xlim(-1.05, 1.05)
            ax.set_xlabel(r'$\xi_M = \frac{\beta_M(k, k+1) - \beta_M(k-1, k)}{2}$')
            ax.grid()
            if i == 0:
                ax.legend(loc='lower center', ncol=1,)
                #  ax.legend(loc='lower center', ncol=ncols,)
            ax.set_title(plotnames[i])

        fig2.axes[0].set_ylabel(r"$N/N_{tot}$")

        #  # in case you want a title
        fig2.suptitle("Mass Growth Fluctuations", fontsize=21)
        fig2.tight_layout(pad=0.5, w_pad=1.5, h_pad=1.0, rect=(0,0,1,0.9))

        #  # in case you don't want a title
        #  fig2.tight_layout(pad=0.5, w_pad=1.5, h_pad=1.0)

        figname="mass_fluctuations-halo-subhalo-"+suffix+".png"
        fig2.savefig(figname, dpi=300)















    if plot_geometry_partbins:
        #-----------------------------------------------------
        # Tree Geometry: Main Branch Length with particle bins
        #-----------------------------------------------------

        for b in range(npartbins):
            ax = fig3.axes[b*2]
            ax.grid()
            ax.set_xlim(0.9, p.nout)
            #  ax.set_ylim(mblmin, 0.5)
            ax.set_ylabel(r"$N + 1$")
            if b == npartbins - 1:
                ax.set_xlabel("Main Branch Length", labelpad=18, fontsize=22)

            #  add redshift axis
            ax.set_xticks([i for i in range(0,p.nout, 10)])
            ax.set_xticklabels(["{0:2d}".format(i) for i in range(0,p.nout, 10)])
            if b==0:
                axt = ax.twiny()
                axt.set_xlim(0.9, p.nout)
                axt.set_xticks([i for i in range(0, p.nout, 10)])
                axt.set_xticklabels(["%.2f" % abs(sd.redshift[i+p.z0]) for i in range(0, p.nout, 10)])
                axt.set_xlabel(r'redshift $z$')

            # label on right side of plot
            axtwin = ax.twinx()
            axtwin.set_yticks([])
            axtwin.set_ylabel(binnames[b], fontsize=12)


        #---------------------------------------
        # Number of Branches with particle bins
        #---------------------------------------

        for b in range(npartbins):
            ax = fig3.axes[b*2+1]
            ax.grid()
            ax.set_xlim(0.9, nbins_nbranches)
            #  ax.set_ylim(nbranchmin, 1)
            ax.set_ylabel(r"$N + 1$")
            #  ax.set_ylabel(r"$N/N_{tot}$ ")
            if b == 0:
                ax.legend()
                #  ax.set_title("Number of Branches", pad=48, fontsize=22)
            if b == npartbins - 1:
                ax.set_xlabel("Number of Branches", labelpad=18, fontsize=22)
            # label on right side of plot
            axtwin = ax.twinx()
            axtwin.set_yticks([])
            axtwin.set_ylabel(binnames[b], fontsize=12)



        fig3.tight_layout(pad=0.1, w_pad=2.5, h_pad=0.5)

        figname = 'tree_geometry-'+suffix+'.png'
        fig3.savefig(figname, dpi=300)






    if plot_main_brench_length_halo_subhalo_partbins:
        #-----------------------------------------------------------------
        # Main Branch Lengths for main haloes/subhaloes and particle bins
        #-----------------------------------------------------------------


        fig4.axes[0].set_title("Main Haloes")
        fig4.axes[1].set_title("Subhaloes")

        for b,ax in enumerate(fig4.axes):
            ax.grid()
            ax.set_xlim(0.9, p.nout)
            #  ax.set_ylim(mblmin, 0.5)
            ax.set_ylabel(r"$N + 1$")

            if b >= 2*npartbins - 2:
                ax.legend()

            #  add redshift axis
            ax.set_xticks([i for i in range(0,p.nout, 10)])
            ax.set_xticklabels(["{0:2d}".format(i) for i in range(0,p.nout, 10)])
            if b < 2:
                axt = ax.twiny()
                axt.set_xlim(0.9, p.nout)
                axt.set_xticks([i for i in range(0, p.nout, 10)])
                axt.set_xticklabels(["%.2f" % abs(sd.redshift[i+p.z0]) for i in range(0, p.nout, 10)])
                axt.set_xlabel(r'redshift $z$')

            # label on right side of plot
            axtwin = ax.twinx()
            axtwin.set_yticks([])
            axtwin.set_ylabel(binnames[b//2])

        figname = 'main-branch-lengths-halo-subhalo-partbins-'+suffix+'.png'

        fig4.tight_layout(pad=0.1, w_pad=2.5, h_pad=1.5)
        fig4.savefig(figname, dpi=300)





    if plot_nbranches_halo_subhalo_particle_bins:
        #---------------------------------------------------------------
        # Number of Branches for main haloes/subhaloes and particle bins
        #---------------------------------------------------------------

        fig5.axes[0].set_title("Main Haloes")
        fig5.axes[1].set_title("Subhaloes")

        for b,ax in enumerate(fig5.axes):
            ax.grid()
            ax.set_xlim(0.9, nbins_nbranches)
            #  ax.set_ylim(nbranchmin, 1)
            #  ax.set_ylabel(r"$N/N_{tot}$ ")
            ax.set_ylabel(r"$N + 1$")
            if b < 2:
                ax.legend()
            if b >= 2*npartbins - 2:
                ax.set_xlabel("Number of Branches", fontsize=16, labelpad=18)

            # label on right side of plot
            axtwin = ax.twinx()
            axtwin.set_yticks([])
            axtwin.set_ylabel(binnames[b//2])


        figname = 'number-of-branches-halo-subhalo-partbins-'+suffix+'.png'
        fig5.tight_layout(pad=0.1, w_pad=2.5, h_pad=1.5)
        fig5.savefig(figname, dpi=300)





    if plot_mass_growth_and_fluctuations_vertical:
        #----------------------------------------------------------------
        # Mass growth and fluctuations without halo/subhalo divide
        #----------------------------------------------------------------


        ax7.set_title("Logarithmic Mass Growth")
        ax8.set_title("Mass Growth Fluctuations")

        ax7.set_ylim(1, 5e3)
        ax8.set_ylim(1, 5e3)
        ax7.set_xlim(-1.05, 1.05)
        ax8.set_xlim(-1.05, 1.05)
        ax7.set_xlabel(r'$\beta_M = \frac{2}{\pi}\arctan\frac{d \log M}{d \log t}$',labelpad=18)
        ax8.set_xlabel(r'$\xi_M = \frac{\beta_M(k, k+1) - \beta_M(k-1, k)}{2}$',labelpad=18)
        ax7.set_ylabel(r"$N+1$")
        ax8.set_ylabel(r"$N+1$")
        ax7.legend(loc='lower center')
        ax8.legend(loc='lower center')
        ax7.grid()
        ax8.grid()

        figname = 'mass-statistics-'+suffix+'.png'
        fig6.tight_layout(pad=0.1, w_pad=2.5, h_pad=1.5)
        fig6.savefig(figname, dpi=300)





    if plot_main_brench_length_halo_subhalo_no_partbins:
        #----------------------------------------------
        # Main Branch Lengths, no particle bins separation
        #----------------------------------------------

        for ax in (ax9, ax10):
            ax.set_ylabel(r"$N+1$")
            ax.set_xlim(0.9, p.nout)
            #  ax.set_ylim(mblmin, 0.5)
            ax.grid()
            ax.legend()
            ax.set_xlabel("Main Branch Length", fontsize=16, labelpad=18)

            #  add redshift axis
            ax.set_xticks([i for i in range(0,p.nout, 10)])
            ax.set_xticklabels(["{0:2d}".format(i) for i in range(0,p.nout, 10)])
            axt = ax.twiny()
            axt.set_xlim(0.9, p.nout)
            axt.set_xticks([i for i in range(0, p.nout, 10)])
            axt.set_xticklabels(["%.2f" % abs(sd.redshift[i+p.z0]) for i in range(0, p.nout, 10)])
            axt.set_xlabel(r'redshift $z$')

        ax9.set_title("Main Haloes")
        ax10.set_title("Subhaloes")

        figname = 'main-branch-lenghts-no-partbins-'+suffix+'.png'

        fig7.tight_layout(pad=0.1, w_pad=2.5, h_pad=1.5)
        fig7.savefig(figname, dpi=300)




    if plot_nbranches_halo_subhalo_no_partbins:
        #--------------------------------------------------
        # Number of Branches for all particle bins
        #--------------------------------------------------

        for ax in (ax11, ax12):
            ax.set_ylabel(r"$N + 1$")
            ax.set_xlim(0.9, nbins_nbranches)
            ax.grid()
            ax.legend()
            ax.set_xlabel("Number of Branches", fontsize=16, labelpad=18)

        ax11.set_title("Main Haloes")
        ax12.set_title("Subhaloes")


        figname = 'number-of-branches-halo-subhalo-no-partbins-'+suffix+'.png'

        fig8.tight_layout(pad=0.1, w_pad=2.5, h_pad=1.5)
        fig8.savefig(figname, dpi=300)







    if plot_nbranches_all_no_partbins:
        #------------------------------------------------------
        # Number of Branches of all clumps, no particle bins
        #------------------------------------------------------

        ax13.set_ylabel(r"$N/N_{tot}$")
        #  ax.set_xlim(0.9, nbins_nbranches)
        ax13.grid()
        ax13.legend()
        ax13.set_xlabel("number of direct progenitors", fontsize=16, labelpad=18)

        figname = 'nbranches-all-no-bins-'+suffix+'.png'

        fig9.tight_layout(pad=0.1, w_pad=2.5, h_pad=1.5)
        fig9.savefig(figname, dpi=300)




    if plot_displacements:
        # -----------------------------
        # Displacements
        # -----------------------------

        ax14.set_ylabel(r"$N + 1$")
        ax14.set_xlabel(r"$\Delta_r$")
        ax14.set_xlim((0.1, max(int(displacements_max+1), 10.0)))

        ymax = 10**int(np.log10(displacements_count_max)+1)
        ax14.set_ylim((1., ymax))

        ax14.plot( [displacements_ninety_percent, displacements_ninety_percent], 
                [1, ymax], 
                'b--', 
                label='90th percentile', 
                linewidth=1
                )
        ax14.plot([displacements_ninety_nine_percent, displacements_ninety_nine_percent], [1, ymax], 'r--', label='99th percentile', linewidth=1)
        ax14.legend()
        ax14.grid()
        figname = 'displacements-'+suffix+'.png'

        fig10.tight_layout()
        fig10.savefig(figname, dpi=300)



    if plot_dlogMdlogt:
        # ---------------------------------------
        # d log M / d log t
        # ---------------------------------------

        ax15.set_ylabel(r"$N + 1$")
        ax15.set_xlabel(r"d log $M$ / d log $t$")
        ax15.set_xscale('symlog', linthresh=edge)
        ax15.set_xlim([-10, 10])
        ax15.legend()
        ax15.grid()
        figname = 'dlogMdlogt-'+suffix+'.png'
        #  fig11.tight_layout()
        fig11.savefig(figname, dpi=300)





if __name__ == "__main__":
    main()
