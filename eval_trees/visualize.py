#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt


myfigsize=(20,6)
nbins = 100




def plot_mass_growth(r):
    """
    r: results object
    """


    fig = plt.figure(figsize=myfigsize)
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    for arr, n, ax in [(r.mg, r.mg_free, ax1), (r.hmg, r.hmg_free, ax2), (r.shmg, r.shmg_free, ax3)] :
        cut = arr[:n]
        hist, bin_edges = np.histogram(cut, bins=100, range=(-1, 1))
        print(hist)
        hist = hist/n # normalize histogram
        bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])


        ax.plot(bin_centers, hist)


    plt.show()

    


