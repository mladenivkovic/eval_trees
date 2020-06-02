#!/usr/bin/python3


#===================================================
# This script plots the results of eval_tree.py
#===================================================

import numpy as np
import matplotlib as mpl
#  mpl.use('Agg')
import matplotlib.pyplot as plt
#  plt.ioff()
#  from matplotlib.font_manager import FontProperties # for legend
import pickle


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








#===================
def main():
#===================
    

    #================================
    # Select files
    #================================

    #--------------------
    # Debug
    #--------------------
    #  allfiles = ['eval_trees-sussing-criteria.pkl']
    #  labelnames = ['label']
    #  suffix = 'debug'
    #  linestyle = ['-']


    #--------------------
    # Actually in use
    #--------------------

    # For no-thresholds

    ntrace = [1, 10, 100, 1000]
    #  ntrace = [1, 10, 50, 100, 200, 500, 1000]
    allfiles = ['eval_trees-no-threshold-ntrace-'+str(i)+'.pkl' for i in ntrace]
    labelnames = [r'$n_{mb}='+str(i)+'$' for i in ntrace]
    #  linestyle = ['-', '-', '--', '--', '-.', '-.', ':']
    linestyle = ['-']*10
    suffix='ntrace'

    #  allfiles = [ 'eval_trees-no-threshold-inclusive-nosaddle.pkl',
    #              'eval_trees-no-threshold-exclusive-nosaddle.pkl',
    #              'eval_trees-no-threshold-inclusive-saddle.pkl',
    #              'eval_trees-no-threshold-exclusive-saddle.pkl' ]
    #  labelnames = [ 'inclusive loosely bound', 'exclusive loosely bound', 'inclusive strictly bound', 'exclusive strictly bound' ]
    #  linestyle = ['-', '--', '-', '--']
    #  suffix='inc-excl'


    #  For sussing thresholds
    #  ntrace = [1, 10, 50, 100, 200, 500, 1000]
    #  ntrace = [100, 1000]
    #  allfiles = ['eval_trees-sussing-criteria-ntrace-'+str(i)+'.pkl' for i in ntrace]
    #  labelnames = [r'$n_{mb}='+str(i)+'$' for i in ntrace]
    #  linestyle = ['-', '-', '--', '--', '-.', '-.', ':']
    #  suffix='ntrace'

    #  allfiles = [ 'eval_trees-sussing-criteria-inclusive-nosaddle.pkl',
    #              'eval_trees-sussing-criteria-exclusive-nosaddle.pkl',
    #              'eval_trees-sussing-criteria-inclusive-saddle.pkl',
    #              'eval_trees-sussing-criteria-exclusive-saddle.pkl' ]
    #  labelnames = [ 'inclusive loosely bound', 'exclusive loosely bound', 'inclusive strictly bound', 'exclusive strictly bound' ]
    #  #  linestyle = ['--', '--', ':', ':']
    #  suffix='inc-excl'




    #  linestyle = ['-', ':', '--', '-.', '-', ':', '--', '-.']
    #  linestyle = ['-']*10
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#e377c2', '#bcbd22', '#17becf']

    # how many particle bins do we have for tree geometry?
    #  binnames = [r'$100-500$ particles', r'$500-1000$ particles', r'$> 1000$ particles']
    binnames = [r'$< 100$ particles', r'$100-500$ particles', r'$500-1000$ particles', r'$> 1000$ particles']
    npartbins = len(binnames)


    # columns for legends
    if len(allfiles) == 4:
        ncols = 2
    else:
        ncols = 3

    alpha = 0.6


    #==========================
    # Set up
    #==========================

    # mass growth with halo-subhalo divide
    # mass_growth.png
    fig1 = plt.figure(1, figsize=(18, 6))
    ax1 = fig1.add_subplot(131)
    ax2 = fig1.add_subplot(132)
    ax3 = fig1.add_subplot(133)

    # mass fluctuations with halo-subhalo divide
    # mass_fluctuations.png
    fig2 = plt.figure(2, figsize=(18, 6))
    ax4 = fig2.add_subplot(131)
    ax5 = fig2.add_subplot(132)
    ax6 = fig2.add_subplot(133)

    # tree geometry
    # tree_geometry.png
    #  fig3 = plt.figure(3, figsize=(14,8)) # for 3 plots
    fig3 = plt.figure(3, figsize=(14,10)) # for 4 plots
    for b in range(2*npartbins):
        fig3.add_subplot(npartbins, 2, b+1)

    # tree geometry; halo-subhalo divide main branch lengths
    # main-branch-lengths-halo-subhalo.png
    fig4 = plt.figure(4, figsize=(16,9))
    for b in range(2*npartbins):
        fig4.add_subplot(npartbins, 2, b+1)

    # tree geometry; halo-subhalo divide number of branches
    # branching-ratio-halo-subhalo.png
    fig5 = plt.figure(5, figsize=(15,8))
    for b in range(2*npartbins):
        fig5.add_subplot(npartbins, 2, b+1)

    # mass growth and fluctuations, no divide
    # mass-statistics.png
    fig6 = plt.figure(6, figsize=(12, 6))
    ax7 = fig6.add_subplot(1, 2, 1)
    ax8 = fig6.add_subplot(1, 2, 2)

    # main brench length for all particle bins
    # main-branch-lengths-all-bins.png
    fig7 = plt.figure(7, figsize=(12, 6))
    ax9 = fig7.add_subplot(1, 2, 1)
    ax10 = fig7.add_subplot(1, 2, 2)

    # number of branches for all particle bins
    # number-of-branches-all-bins.png
    fig8 = plt.figure(8, figsize=(12, 6))
    ax11 = fig8.add_subplot(1, 2, 1)
    ax12 = fig8.add_subplot(1, 2, 2)



    # tweaking stuff to store
    mblmin = 1000
    nbranchmin = 1000
    fluctmin = 10000
    growthmin = 10000

    #  nbins_nbranches = 100 # how many bins to use for number of branches; = max number of branches expected
    nbins_nbranches = 1200 # how many bins to use for number of branches; = max number of branches expected

    #=============================
    # Read in data
    #=============================

    for f,srcfname in enumerate(allfiles):

        srcfile = open(srcfname, 'rb')
        p, sd, r = pickle.load(srcfile) 
        srcfile.close()


        #-------------------------------
        # Mass Evolution
        #-------------------------------

        norm = r.mg_free
        for arr, n, ax in [(r.mg, r.mg_free, ax1), (r.hmg, r.hmg_free, ax2), (r.shmg, r.shmg_free, ax3), (r.mg, r.mg_free, ax7)] :
            cut = arr[:n]
            hist, bin_edges = np.histogram(cut, bins=100, range=(-1, 1))
            hist = hist/norm # normalize histogram
            growthmin = min(growthmin, hist[hist>0].min())
            bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])

            ax.semilogy(bin_centers, hist,
                label = labelnames[f], 
                c=colors[f], 
                ls = linestyle[f],
                alpha = alpha,
                )


        norm = r.mf_free
        for arr, n, ax in [(r.mf, r.mf_free, ax4), (r.hmf, r.hmf_free, ax5), (r.shmf, r.shmf_free, ax6), (r.mf, r.mf_free, ax8)] :
            cut = arr[:n]
            hist, bin_edges = np.histogram(cut, bins=100, range=(-1, 1))
            hist = hist/norm # normalize histogram
            fluctmin = min(fluctmin, hist[hist>0].min())
            bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
 
            ax.semilogy(bin_centers, hist,
                label = labelnames[f], 
                c=colors[f], 
                ls = linestyle[f],
                alpha = alpha,
                )





        #---------------------------------
        # Geometry
        #---------------------------------

        # figure subplot setup:
        #   1   2
        #   3   4
        #   5   6
        #   7   8
        # -> nr of branches: 2*(b+1)    -1 for list index of fig3.axes
        # -> MBL: 2*b + 1               -1 for list index of fig3.axes


        # plot branch lengths
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






        hbranchlen_all = np.ones(p.nout) # plot N + 1
        sbranchlen_all = np.ones(p.nout) # plot N + 1

        # plot branch lengths for main haloes only
        for b in range(npartbins):
            ax = fig4.axes[b*2]     
            norm = r.branchlen_free_main[b]
            if norm == 0:
                print("main haloes branch lengths got norm=0, bin is", binnames[b])
                continue

            cut = r.branchlengths_main[b][:r.branchlen_free_main[b]]
            hist, bin_edges = np.histogram(cut, bins = p.nout, range=(1, p.nout))
            hbranchlen_all += hist
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
        ax9.semilogy(bins_left, hbranchlen_all,
            label = labelnames[f], 
            c=colors[f], 
            ls = linestyle[f],
            alpha = alpha,
            )


        # plot branch lengths for sub haloes only
        for b in range(npartbins):
            ax = fig4.axes[b*2+1]     
            norm = r.branchlen_free_sub[b]
            if norm == 0:
                continue

            cut = r.branchlengths_sub[b][:r.branchlen_free_sub[b]]
            hist, bin_edges = np.histogram(cut, bins = p.nout, range=(1, p.nout))
            sbranchlen_all += hist
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
        ax10.semilogy(bins_left, sbranchlen_all,
            label = labelnames[f], 
            c=colors[f], 
            ls = linestyle[f],
            alpha = alpha,
            )



        hnbranches_all = np.ones(nbins_nbranches) # plot N + 1
        snbranches_all = np.ones(nbins_nbranches) # plot N + 1

        # plot number of branches for main haloes only
        for b in range(npartbins):
            ax = fig5.axes[b*2]
            norm = r.nbranch_free_main[b]
            if norm == 0:
                continue

            cut = r.nbranches_main[b][:r.nbranch_free_main[b]]
            hist, bin_edges = np.histogram(cut, bins = nbins_nbranches, range=(1, nbins_nbranches))
            hnbranches_all += hist
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

        ax11.loglog(bins_left, hnbranches_all,
            label = labelnames[f], 
            c=colors[f], 
            ls = linestyle[f],
            alpha = alpha,
            )


        # plot number of branches for sub haloes only
        for b in range(npartbins):
            ax = fig5.axes[b*2+1]
            norm = r.nbranch_free_sub[b]
            if norm == 0:
                continue

            cut = r.nbranches_sub[b][:r.nbranch_free_sub[b]]
            hist, bin_edges = np.histogram(cut, bins = nbins_nbranches, range=(1, nbins_nbranches))
            snbranches_all += hist
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
        ax12.loglog(bins_left, snbranches_all,
            label = labelnames[f], 
            c=colors[f], 
            ls = linestyle[f],
            alpha = alpha,
            )



    #==================================
    # Tweak plots and save figures
    #==================================


    #-------------------------
    # Mass growth
    #-------------------------

    plotnames = ['A: All Haloes', 'B: Main Haloes', 'C: Subhaloes']

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
    fig1.savefig(figname, dpi=300, form='png')

    

    #-------------------------
    # Mass growth flucts
    #-------------------------

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

    figname="mass_fluctuations-"+suffix+".png"
    fig2.savefig(figname, dpi=300, form='png')
    













    #-------------------------------------
    # Main Branch Length
    #-------------------------------------

    for b in range(npartbins):
        ax = fig3.axes[b*2]
        ax.grid()
        ax.set_xlim(0.9, p.nout)
        #  ax.set_ylim(mblmin, 0.5)
        ax.set_ylabel(r"$N + 1$")
        #  ax.set_ylabel(r"$N/N_{tot}$")
        #  if b == 0:
        #      ax.set_title("Main Branch Length", pad=18, fontsize=22)
        if b == npartbins - 1:
            ax.set_xlabel("Main Branch Length", labelpad=18, fontsize=22)
            #  ax.legend() # TODO: for my threshold, skip legend on MBL

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


    #-------------------------------------
    # Number of Branches
    #-------------------------------------

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
    fig3.savefig(figname, dpi=300, form='png')






    #------------------------------------------------
    # Main Branch Lengths for main haloes/subhaloes
    #------------------------------------------------

    fig4.axes[0].set_title("Main Haloes")
    fig4.axes[1].set_title("Subhaloes")

    for b,ax in enumerate(fig4.axes):
        ax.grid()
        ax.set_xlim(0.9, p.nout)
        #  ax.set_ylim(mblmin, 0.5)
        ax.set_ylabel(r"$N + 1$")

        if b >= 2*npartbins - 2:
            #  ax.set_xlabel("Main Branch Length", fontsize=16, labelpad=18)
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

    figname = 'main-branch-lengths-halo-subhalo-'+suffix+'.png'

    fig4.tight_layout(pad=0.1, w_pad=2.5, h_pad=1.5)
    fig4.savefig(figname, dpi=300, form='png')




    #------------------------------------------------
    # Number of Branches for main haloes/subhaloes
    #------------------------------------------------

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


    figname = 'branching-ratio-halo-subhalo-'+suffix+'.png'
    fig5.tight_layout(pad=0.1, w_pad=2.5, h_pad=1.5)
    fig5.savefig(figname, dpi=300, form='png')





    #----------------------------------------------------------------
    # Mass growth and fluctuations without halo/subhalo divide
    #----------------------------------------------------------------

    ax7.set_title("Logarithmic Mass Growth")
    ax8.set_title("Mass Growth Fluctuations")

    ax7.set_ylim(growthmin, 1e-1)
    ax8.set_ylim(fluctmin, 1e-1)
    ax7.set_xlim(-1.05, 1.05)
    ax8.set_xlim(-1.05, 1.05)
    ax7.set_xlabel(r'$\beta_M = \frac{2}{\pi}\arctan\frac{d \log M}{d \log t}$',labelpad=18)
    ax8.set_xlabel(r'$\xi_M = \frac{\beta_M(k, k+1) - \beta_M(k-1, k)}{2}$',labelpad=18)
    ax7.set_ylabel(r"$N/N_{tot}$")
    ax8.set_ylabel(r"$N/N_{tot}$")
    ax7.legend(loc='lower center')
    ax8.legend(loc='lower center')
    ax7.grid()
    ax8.grid()


    figname = 'mass-statistics-'+suffix+'.png'
    fig6.tight_layout(pad=0.1, w_pad=2.5, h_pad=1.5)
    fig6.savefig(figname, dpi=300, form='png')





    #-----------------------------------------------------------------------
    # Main Branch Lengths for all particle bins
    #-----------------------------------------------------------------------


    for ax in (ax9, ax10):
        ax.set_ylabel(r"$N$")
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



    figname = 'main-branch-lenghts-all-bins-'+suffix+'.png'

    fig7.tight_layout(pad=0.1, w_pad=2.5, h_pad=1.5)
    fig7.savefig(figname, dpi=300, form='png')




    #-----------------------------------------------------------------------
    # Number of Branches for all particle bins
    #-----------------------------------------------------------------------


    for ax in (ax11, ax12):
        ax.set_ylabel(r"$N$")
        ax.set_xlim(0.9, nbins_nbranches)
        ax.grid()
        ax.legend()
        ax.set_xlabel("Number of Branches", fontsize=16, labelpad=18)

    ax11.set_title("Main Haloes")
    ax12.set_title("Subhaloes")


    figname = 'number-of-branches-all-bins-'+suffix+'.png'

    fig8.tight_layout(pad=0.1, w_pad=2.5, h_pad=1.5)
    fig8.savefig(figname, dpi=300, form='png')










if __name__ == "__main__":
    main()
