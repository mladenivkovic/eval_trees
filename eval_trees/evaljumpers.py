#!/usr/bin/env python3


#===================================================
# This script plots analytics for jumpers from
# eval_tree.py
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


# select which plots to make
plot_jumping_distance = True
plot_mass_statistics = True
plot_with_different_mass_thresholds = True




part_thresh = [100, 200, 500, 1000]
#  part_thresh = [0, 100, 200, 500, 1000]
mpart = 1.553e+09


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

hist_bins = 200
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#e377c2', '#bcbd22', '#17becf']
alpha = 0.6






#===================
def main():
#===================
    
    allfiles, labelnames, suffix, linestyle = file_selection(plot_selection)

    # columns for legends
    if len(allfiles) == 4:
        ncols = 2
    else:
        ncols = 3


    if plot_jumping_distance:
        # histogram of snapshots jumped over
        fig1 = plt.figure(1, figsize=(6, 6))
        ax1 = fig1.add_subplot(111)

    if plot_mass_statistics:
        # progenitor and descendant mass histograms
        fig2 = plt.figure(2, figsize=(6, 18))
        ax21 = fig2.add_subplot(311)
        ax22 = fig2.add_subplot(312)
        ax23 = fig2.add_subplot(313)

    if plot_with_different_mass_thresholds:
        thresholdfigs = [plt.figure(666+i, figsize=(12, 12)) for i in range(len(allfiles))]
        for fig in thresholdfigs:
            fig.add_subplot(221)
            fig.add_subplot(222)
            fig.add_subplot(223)
            fig.add_subplot(224)


        # min/max for axes
        descmassmin = 1e30
        descmassmax = 0.
        deschistmax = 1.
        progmassmin = 1e30
        progmassmax = 0.
        proghistmax = 1.
        ratiomin = 1e30
        ratiomax = -1.
        ratiohistmax = 1.
 





    #=============================
    # Read in data
    #=============================

    for f,srcfname in enumerate(allfiles):

        srcfile = open(srcfname, 'rb')
        p, sd, r = pickle.load(srcfile) 
        srcfile.close()



        #--------------------------------
        # Extract data in useable form
        #--------------------------------
        jumps = np.array([j.snapshot_desc - j.snapshot_prog for j in r.jumper_results])
        mprog = np.array([j.mass_prog for j in r.jumper_results])
        mdesc = np.array([j.mass_desc for j in r.jumper_results])
        njumpers = len(r.jumper_results)
        print("Read in", njumpers, "jumpers")


        #---------------------------------------
        # Histogram of snapshots jumped over
        #---------------------------------------

        if plot_jumping_distance:
            hist, bin_edges = np.histogram(jumps, bins=p.noutput, range=(0, p.noutput))
            hist += 1 # plot N + 1
            ax1.semilogy(
                        bin_edges[:-1], hist, 
                        label=labelnames[f], 
                        c=colors[f],
                        ls = linestyle[f],
                        alpha = alpha,
                    )



        #--------------------------------------------
        # Progenitor and descendant mass histograms
        #--------------------------------------------

        if plot_mass_statistics:
            bins = np.logspace(np.log10(mprog.min()*0.5), np.log10(mprog.max()*2), hist_bins)
            hist_prog, bin_edges = np.histogram(mprog, bins=bins)
            hist_prog += 1 # plot N+1
            bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
            mask = hist_prog > 1
            ax21.loglog(
                            bin_centres[mask], hist_prog[mask], 
                            label=labelnames[f], 
                            c=colors[f],
                            ls = linestyle[f],
                            alpha = alpha,
                        )

            bins = np.logspace(np.log10(mdesc.min()*0.5), np.log10(mdesc.max()*2), hist_bins)
            hist_desc, bin_edges = np.histogram(mdesc, bins=bins)
            hist_desc += 1 # plot N+1
            bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
            mask = hist_desc > 1
            ax22.loglog(
                            bin_centres[mask], hist_desc[mask], 
                            label=labelnames[f], 
                            c=colors[f],
                            ls = linestyle[f],
                            alpha = alpha,
                        )

            ratio = np.array([mdesc[i] / mprog[i] for i in range(njumpers)])
            bins = np.logspace(np.log10(ratio.min()*0.5), np.log10(ratio.max()*2), hist_bins)
            hist_ratio, bin_edges = np.histogram(ratio, bins=bins)
            hist_ratio += 1 # plot N+1
            mask = hist_ratio > 1
            bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
            ax23.loglog(
                        bin_centres[mask], hist_ratio[mask], 
                        label=labelnames[f], 
                        c=colors[f],
                        ls = linestyle[f],
                        alpha = alpha,
                        )






        # Do thresholds plot for each file individually
        #----------------------------------------------

        if plot_with_different_mass_thresholds:
            axT1, axT2, axT3, axT4 = thresholdfigs[f].axes
            for i, t in enumerate(part_thresh):
                thresh = t * mpart
                mask = np.logical_and(mprog > thresh, mdesc > thresh)

                label=r"threshold = {0:4d} $m_p$".format(t)

                # plot progenitor masses
                #-------------------------
                bins = np.logspace(np.log10(mprog[mask].min()*0.5), np.log10(mprog[mask].max()*2), hist_bins)
                hist_prog, bin_edges = np.histogram(mprog[mask], bins=bins)
                hist_prog += 1 # plot N+1
                bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
                #  mask_plot = hist_prog > 1
                axT1.loglog(
                                bin_centres, hist_prog, 
                                #  bin_centres[mask_plot], hist_prog[mask_plot],
                                label=label, 
                                c=colors[i],
                                ls = linestyle[i],
                                alpha = alpha,
                            )

                progmassmin = min(progmassmin, bin_centres[0])
                progmassmax = max(progmassmax, bin_centres[-1])
                proghistmax = max(proghistmax, hist_prog.max())

                # plot descendant masses
                #-------------------------
                bins = np.logspace(np.log10(mdesc[mask].min()*0.5), np.log10(mdesc[mask].max()*2), hist_bins)
                hist_desc, bin_edges = np.histogram(mdesc[mask], bins=bins)
                hist_desc += 1 # plot N+1
                bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
                #  mask_plot = hist_desc > 1
                axT2.loglog(
                                bin_centres, hist_desc, 
                                #  bin_centres[mask_plot], hist_desc[mask_plot],
                                label=label, 
                                c=colors[i],
                                ls = linestyle[i],
                                alpha = alpha,
                            )
                descmassmin = min(descmassmin, bin_centres[0])
                descmassmax = max(descmassmax, bin_centres[-1])
                deschistmax = max(deschistmax, hist_desc.max())

                # plot jumping distance
                #-------------------------
                hist, bin_edges = np.histogram(jumps[mask], bins=p.noutput, range=(0, p.noutput))
                hist += 1 # plot N + 1
                axT3.semilogy(
                            bin_edges[:-1], hist, 
                            label=label, 
                            c=colors[i],
                            ls = linestyle[i],
                            alpha = alpha,
                        )

                # plot mass ratios
                #-------------------------
                ratio = mdesc[mask] / mprog[mask]
                bins = np.logspace(np.log10(ratio.min()*0.5), np.log10(ratio.max()*2), hist_bins)
                hist_ratio, bin_edges = np.histogram(ratio, bins=bins)
                hist_ratio += 1 # plot N+1
                bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
                #  mask_plot = hist_ratio > 1
                axT4.loglog(
                            bin_centres, hist_ratio, 
                            #  bin_centres[mask_plot], hist_ratio[mask_plot],
                            label=label, 
                            c=colors[i],
                            ls = linestyle[i],
                            alpha = alpha,
                            )
                ratiomin = min(ratiomin, bin_centres[0])
                ratiomax = max(ratiomax, bin_centres[-1])
                ratiohistmax = max(ratiohistmax, hist_ratio.max())







    #-----------------------------------------
    # TWEAK AND SAVE FIGURES
    #-----------------------------------------


    if plot_jumping_distance:
        smalltitlefontsize = 14

        ax1.set_xlabel(r"snapshot$_{desc}$ - snapshot$_{prog}$")
        ax1.set_ylabel(r"$N+1$")
        ax1.set_title(r"Histogram of snapshot number difference of all jumpers", fontsize=smalltitlefontsize)
        ax1.grid()
        ax1.legend()
        fig1.savefig("jumper_distance-"+suffix+".png", dpi=300)
        fig1.tight_layout()


    if plot_mass_statistics:
        for ax in fig2.axes:
            ax.set_ylabel(r"$N+1$")
            ax.grid()
            ax.legend()

        ax21.set_xlabel(r"mass $[M_\odot]$")
        ax21.set_title(r"Masses of Jumper Progenitors")

        ax22.set_xlabel(r"mass $[M_\odot]$")
        ax22.set_title(r"Masses of Jumper Descendants")

        ax23.set_xlabel(r"$m_{desc}$/$m_{prog}$")
        ax23.set_title(r"Jumper Descendant to Jumper Progenitor Mass Ratio", fontsize=smalltitlefontsize)

        fig2.tight_layout()
        fig2.savefig("jumper_mass_statistics-"+suffix+".png", dpi=300)



    if plot_with_different_mass_thresholds:

        for f, srcfile in enumerate(allfiles):

            # Tweak plots

            fig = thresholdfigs[f]
            axT1, axT2, axT3, axT4 = fig.axes

            for ax in fig.axes:
                ax.set_ylabel(r"$N+1$")
                ax.grid()
                ax.legend()

            smalltitlefontsize = 14

            axT1.set_xlabel(r"mass $[M_\odot]$")
            axT1.set_title(r"Masses of Jumper Progenitors")
            axT1.set_xlim(progmassmin/2, progmassmax*2)
            axT1.set_ylim(0.6, proghistmax*2)

            axT2.set_xlabel(r"mass $[M_\odot]$")
            axT2.set_title(r"Masses of Jumper Descendants")
            axT2.set_xlim(descmassmin/2, descmassmax*2)
            axT2.set_ylim(0.6, deschistmax*2)

            axT3.set_xlabel(r"snapshot$_{desc}$ - snapshot$_{prog}$")
            axT3.set_title(r"Histogram of snapshot number difference of all jumpers", fontsize=smalltitlefontsize)

            axT4.set_xlabel(r"$m_{desc}$/$m_{prog}$")
            axT4.set_title(r"Jumper Descendant to Jumper Progenitor Mass Ratio", fontsize=smalltitlefontsize)
            axT4.set_xlim(ratiomin/2, ratiomax*2)
            axT4.set_ylim(0.6, ratiohistmax*2)


            # Generate title and filename
            if plot_selection.do_debug:
                title = "low res simulation for debugging and dev"
                filename = "jumper_results_with_mass_thresholds-debugging.png"
            elif plot_selection.do_sussing_ntrace:
                print("ERROR: TODO: NEED TO GENERATE FILENAME AND FIGTITLE")
                quit(1)
            elif plot_selection.do_sussing_incexcl:
                print("ERROR: TODO: NEED TO GENERATE FILENAME AND FIGTITLE")
                quit(1)
            elif plot_selection.do_my_threshold_ntrace:
                nmb = [int(s) for s in srcfile[:-4].split('-') if s.isdigit()]
                if len(nmb) != 1:
                    raise ValueError("got len(nmb) != 1. nmb=", nmb)
                else:
                    nmb = nmb[0]
                title = r"Stats using different mass thresholds for $n_{mb}=$"+"{0:d}".format(nmb)
                filename = "jumper_results_with_mass_thresholds-ntrace-{0:d}".format(nmb)
                #  filename = "jumper_results_with_mass_thresholds-ntrace-{0:d}-including-zero-threshold".format(nmb)
            elif plot_selection.do_my_threshold_incexcl:
                print("ERROR: TODO: NEED TO GENERATE FILENAME AND FIGTITLE")
                quit(1)

            fig.suptitle(title)
            fig.tight_layout()
            fig.savefig(filename, dpi=300)
            














if __name__ == "__main__":
    main()
