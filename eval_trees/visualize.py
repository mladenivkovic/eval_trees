#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt


myfigsize=(20,12)
nbins = 50




#==============================
def plot_mass_growth(r):
#==============================
    """
    r: results object
    """


    fig = plt.figure(figsize=myfigsize)
    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    ax1.set_ylabel("Mass Growth")

    norm = r.mg_free

    for arr, n, ax in [(r.mg, r.mg_free, ax1), (r.hmg, r.hmg_free, ax2), (r.shmg, r.shmg_free, ax3)] :
        cut = arr[:n]
        hist, bin_edges = np.histogram(cut, bins=100, range=(-1, 1))
        hist = hist/norm # normalize histogram
        bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])


        ax.semilogy(bin_centers, hist)

    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)
    ax4.set_ylabel("Mass Growth Fluctuations")

    norm = r.mf_free
    for arr, n, ax in [(r.mf, r.mf_free, ax4), (r.hmf, r.hmf_free, ax5), (r.shmf, r.shmf_free, ax6)] :
        cut = arr[:n]
        hist, bin_edges = np.histogram(cut, bins=100, range=(-1, 1))
        hist = hist/norm # normalize histogram
        bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
        ax.semilogy(bin_centers, hist)



    plt.show()

    



#=============================
def plot_geometry(p, r):
#=============================
    """
    Plot main branch length and number of branches histograms 

        p: params object
        r: results object
    """

    npartbins = len(r.branchlengths)
    plotindex = 0

    fig = plt.figure()

    # plot branch lengths
    for b in range(npartbins):
        plotindex += 1
        ax = fig.add_subplot(2, npartbins, plotindex)
        norm = r.branchlen_free[b]

        cut = r.branchlengths[b][:r.branchlen_free[b]]
        hist, bin_edges = np.histogram(cut, bins=p.noutput, range=(0, p.noutput))
        hist = hist/norm
        bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])

        ax.semilogy(bin_centers, hist)
        if b==0:
            ax.set_ylabel("Main Branch Lengths")


    # find max branches for plot boundaries
    maxbranches = max([np.asscalar(r.nbranches[b][:r.nbranch_free[b]].max()) for b in range(npartbins)])

    # plot number of branches
    for b in range(npartbins):
        plotindex += 1
        ax = fig.add_subplot(2, npartbins, plotindex)
        norm = r.nbranch_free[b]

        cut = r.nbranches[b][:r.nbranch_free[b]]
        hist, bin_edges = np.histogram(cut, bins=maxbranches, range=(1, maxbranches))
        hist = hist/norm
        bin_left = bin_edges[:-1]

        ax.loglog(bin_left, hist)
        if b==0:
            ax.set_ylabel("Nr of Branches")


    plt.show()
