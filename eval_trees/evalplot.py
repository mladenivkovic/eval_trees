#!/usr/bin/python3


#===================================================
# This script plots the results of eval_tree.py
#===================================================

import numpy as np
import matplotlib as mpl
#  mpl.use('Agg')
import matplotlib.pyplot as plt
#  plt.ioff()
from matplotlib.font_manager import FontProperties # for legend
import pickle

# use LaTeX text
from matplotlib import rc
rc('font', **{'family':['serif'],
    'serif':['Computer Modern Roman'],
    'monospace': ['Computer Modern Typewriter'], 'size':14})
rc('text', usetex=True)

fontP=FontProperties()
fontP.set_size('small') 






#===============================================
def save_figures(figname, fig, rect=None):
#===============================================
    """
    Saves figure fig as figname as .png
    """

    if rect is not None:
        fig.tight_layout(rect=rect)
    else:
        fig.tight_layout()
    fig.savefig(figname, dpi=300, form='png')
    print("Saved", figname)

    return






#===================
def main():
#===================
    

    #================================
    # Select files
    #================================

    #--------------------
    # Debug
    #--------------------
    #  allfiles = ['eval_trees.pkl']
    #  labelnames = ['label']
    #  suffix = 'debug'


    #--------------------
    # Actually in use
    #--------------------
    #  ntrace = [1, 10, 50, 100, 200, 500, 1000]
    #  allfiles = ['eval_trees-ntrace-'+str(i)+'.pkl' for i in ntrace]
    #  labelnames = [r'$n_{mb}='+str(i)+'$' for i in ntrace]
    #  linestyle = ['-', '-', '--', '--', '-.', '-.', ':']
    #  suffix='ntrace'

    allfiles = [ 'eval_trees-inclusive-nosaddle.pkl',  'eval_trees-exclusive-nosaddle.pkl', 'eval_trees-inclusive-saddle.pkl', 'eval_trees-exclusive-saddle.pkl' ]
    labelnames =[ 'inclusive loosely bound', 'exclusive loosely bound', 'inclusive strictly bound', 'exclusive strictly bound' ]
    linestyle = ['--', '--', ':', ':']
    suffix='inc-excl'

    #  linestyle = ['-', ':', '--', '-.', '-', ':', '--', '-.']
    #  linestyle = ['-']*10
    #  linestyle = [':']*8
    #  linestyle = ['--']*8
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#e377c2', '#bcbd22', '#17becf']

    # columns for legends
    if len(allfiles) == 4:
        ncols = 2
    else:
        ncols = 3

    linewidth = 2


    #==========================
    # Set up
    #==========================

    fig1 = plt.figure(1, figsize=(18, 6))
    ax1 = fig1.add_subplot(131)
    ax2 = fig1.add_subplot(132)
    ax3 = fig1.add_subplot(133)

    fig2 = plt.figure(2, figsize=(18, 6))
    ax4 = fig2.add_subplot(131)
    ax5 = fig2.add_subplot(132)
    ax6 = fig2.add_subplot(133)

    fig3 = plt.figure(3, figsize=(18,12))
    npartbins = 3   # too lazy to read in results object here; see visualize.py
    for b in range(2*(npartbins+1)):
        fig3.add_subplot(npartbins+1, 2, b+1)


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
        for arr, n, ax in [(r.mg, r.mg_free, ax1), (r.hmg, r.hmg_free, ax2), (r.shmg, r.shmg_free, ax3)] :
            cut = arr[:n]
            hist, bin_edges = np.histogram(cut, bins=100, range=(-1, 1))
            hist = hist/norm # normalize histogram
            growthmin = min(growthmin, hist[hist>0].min())
            bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])

            ax.semilogy(bin_centers, hist,
                label = labelnames[f], c=colors[f], ls = linestyle[f],
                lw=linewidth)


        norm = r.mf_free
        for arr, n, ax in [(r.mf, r.mf_free, ax4), (r.hmf, r.hmf_free, ax5), (r.shmf, r.shmf_free, ax6)] :
            cut = arr[:n]
            hist, bin_edges = np.histogram(cut, bins=100, range=(-1, 1))
            hist = hist/norm # normalize histogram
            fluctmin = min(fluctmin, hist[hist>0].min())
            bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
 
            ax.semilogy(bin_centers, hist,
                label = labelnames[f], c=colors[f], ls = linestyle[f],
                lw=linewidth)





        #---------------------------------
        # Geometry
        #---------------------------------

        npartbins = len(r.branchlengths)

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

            cut = r.branchlengths[b][:r.branchlen_free[b]]
            hist, bin_edges = np.histogram(cut, bins = p.nout, range=(1, p.nout))
            hist = hist/norm
            mblmin = min(mblmin, hist[hist>0].min())
            #  bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
            bins_left = bin_edges[:-1]

            ax.semilogy(bins_left, hist,
                label = labelnames[f], c=colors[f], ls = linestyle[f],
                lw=linewidth)


        # plot number of branches
        for b in range(npartbins):
            ax = fig3.axes[b*2+1]
            norm = r.nbranch_free[b]

            cut = r.nbranches[b][:r.nbranch_free[b]]
            hist, bin_edges = np.histogram(cut, bins = nbins_nbranches, range=(1, nbins_nbranches))
            hist = hist/norm
            nbranchmin = min(nbranchmin, hist[hist>0].min())
            #  bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
            bins_left = bin_edges[:-1]

            ax.loglog(bins_left, hist,
                label = labelnames[f], c=colors[f], ls = linestyle[f],
                lw=linewidth)






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
        ax.legend(loc='lower center', ncol=ncols,prop=fontP)
        ax.set_title(plotnames[i])

    fig1.suptitle("Logarithmic Mass Growth")

    figname="mass_growth-"+suffix+".png"
    save_figures(figname, fig1, rect=[0, 0.03, 1, 0.95])

    

    #-------------------------
    # Mass growth flucts
    #-------------------------

    for i, ax in enumerate(fig2.axes):
        # for inclusive/exclusive
        ax.set_ylim(fluctmin, 1e-1)
        ax.set_xlim(-1.05, 1.05)
        ax.set_xlabel(r'$\xi_M = \frac{\beta_M(k, k+1) - \beta_M(k-1, k)}{2}$')
        ax.grid()
        ax.legend(loc='lower center', ncol=ncols,prop=fontP)
        ax.set_title(plotnames[i])

    fig2.suptitle("Mass Growth Fluctuations")

    figname="mass_fluctuations-"+suffix+".png"
    save_figures(figname, fig2, rect=[0, 0.03, 1, 0.95])













    binnames = [r'$< 100$ particles', r'$100-500$ particles', r'$500-1000$ particles', r'$> 1000$ particles']

    #-------------------------------------
    # Main Branch Length
    #-------------------------------------

    for b in range(npartbins):
        ax = fig3.axes[b*2]
        ax.grid()
        ax.set_xlim(0.9, p.nout)
        ax.set_ylim(mblmin, 0.5)
        ax.set_ylabel(r"$N/N_{tot}$")
        if b == npartbins - 1:
            ax.set_xlabel("Main Branch Length")
            ax.legend()

        #  add redshift axis
        ax.set_xticks([i for i in range(0,p.nout, 10)])
        ax.set_xticklabels(["{0:2d}".format(i) for i in range(0,p.nout, 10)])
        if b==0:
            axt = ax.twiny()
            axt.set_xlim(-0.5, p.nout)
            axt.set_xticks([i for i in range(0, p.nout, 10)])
            axt.set_xticklabels(["%.2f" % abs(sd.redshift[i+p.z0]) for i in range(0, p.nout, 10)])
            axt.set_xlabel(r'redshift $z$')

        # label on right side of plot
        axtwin = ax.twinx()
        axtwin.set_yticks([])
        axtwin.set_ylabel(binnames[b])


    #-------------------------------------
    # Number of Branches
    #-------------------------------------

    for b in range(npartbins):
        ax = fig3.axes[b*2+1]
        ax.grid()
        ax.set_xlim(0.9, nbins_nbranches)
        ax.set_ylim(nbranchmin, 1)
        ax.set_ylabel(r"$N/N_{tot}$ ")
        if b == 0:
            ax.legend()
        if b == npartbins - 1:
            ax.set_xlabel("Number of Branches")
        # label on right side of plot
        axtwin = ax.twinx()
        axtwin.set_yticks([])
        axtwin.set_ylabel(binnames[b])


    figname = 'tree_geometry-'+suffix+'.png'
    save_figures(figname, fig3)


if __name__ == "__main__":
    main()
