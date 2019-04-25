#!/usr/bin/python3


#===================================================
# This script plots the results of eval_tree.py
#===================================================

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.font_manager import FontProperties # for legend

# use LaTeX text
from matplotlib import rc
rc('font', **{'family':['serif'],
    'serif':['Computer Modern Roman'],
    'monospace': ['Computer Modern Typewriter'], 'size':14})
rc('text', usetex=True)

fontP=FontProperties()
fontP.set_size('small') 


#=========================================================
def get_line_to_array(linestring, maxcounter, toint=False):
#=========================================================
    """
    Translate string line separated by commas to numpy array of integers (if toint=True) or floats (if toint=False)
    simultaneously find highest value in array
    """

    import numpy as np

    linestring = linestring.strip()
    outlist = linestring.split(" ")

    if toint:
        for i, val in enumerate(outlist):
            outlist[i] = int(val.strip())
            maxcounter = max(maxcounter, outlist[i])
        out = np.array(outlist, dtype='int')

    else:
        for i, val in enumerate(outlist):
            outlist[i] = float(val.strip())
            maxcounter = max(maxcounter, outlist[i])

        out = np.array(outlist, dtype='float')

    return out, maxcounter






#==================================
def save_figures(figname, fig):
#==================================
    """
    Saves figure fig as figname as .pdf
    """

    fig.tight_layout()
    fig.savefig(figname, dpi=300, form='pdf')
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
    #  allfiles = ['eval_trees.txt']
    #  labelnames = ['label']


    #--------------------
    # Actually in use
    #--------------------
    #  ntrace = [1, 10, 50, 100, 200, 500, 1000]
    #  allfiles = ['eval_trees-ntrace-'+str(i)+'.txt' for i in ntrace]
    #  labelnames = [r'$n_{mb}='+str(i)+'$' for i in ntrace]
    #  linestyle = ['-', '-', '--', '--', '-.', '-.', ':']

    allfiles = [ 'eval_trees-inclusive-nosaddle.txt',  'eval_trees-exclusive-nosaddle.txt', 'eval_trees-inclusive-saddle.txt', 'eval_trees-exclusive-saddle.txt' ]
    labelnames =[ 'inclusive loosely bound', 'exclusive loosely bound', 'inclusive strictly bound', 'exclusive strictly bound' ]
    linestyle = ['--', '--', ':', ':']

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



    #==========================
    # Set up
    #==========================



    fig1 = plt.figure(1, figsize=(18,6))
    ax1 = fig1.add_subplot(131)
    ax2 = fig1.add_subplot(132, sharey=ax1)
    ax3 = fig1.add_subplot(133, sharey=ax1)
    axgroup1 = [ax1, ax2, ax3]

    fig2 = plt.figure(2, figsize=(18,6))
    ax4 = fig2.add_subplot(131)
    ax5 = fig2.add_subplot(132, sharey=ax4)
    ax6 = fig2.add_subplot(133, sharey=ax4)
    axgroup2 = [ax4, ax5, ax6]

    maxlen = 0
    maxcount1 = 0
    maxnbranch = 0
    maxcount2 = 0
    mincount1 = 1
    mincount2 = 1
    mgmax = 0
    mgmin = 1
    mgcountmin=1
    mgcountmax=0
    mfmax = 0
    mfmin = 0
    mfcountmin = 1
    mfcountmax = 0
    displmax = 0
    displmin = 1
    dcountmax = 0
    dcountmin = 1





    #=============================
    # Read in data
    #=============================

    for f,srcfname in enumerate(allfiles):

        srcfile = open(srcfname, 'r')

        # skip 1 lines first
        for i in range(1):
            srcfile.readline()


        #-------------------------
        # Mass growth
        #-------------------------
        
        for ax in axgroup1:

            mass_growth_counts, mgcountmax = get_line_to_array(srcfile.readline(), mgcountmax, False)
            mgcountmin = min(mass_growth_counts[mass_growth_counts>0].min(), mgcountmin)
            # replace zeroes with value orders of magnitude smaller
            mass_growth_counts[mass_growth_counts==0] = mgcountmin*1e-6

            mass_growth, mgmax = get_line_to_array(srcfile.readline(), mgmax, False)
            mgmin = min(mgmin, mass_growth.min())

            ax.semilogy(mass_growth, mass_growth_counts, label=labelnames[f], lw=4, c=colors[f], ls=linestyle[f], alpha=0.7)


        #----------------------
        # mass fluctuations
        #----------------------

        for ax in axgroup2:

            mass_fluct_counts, mfcountmax = get_line_to_array(srcfile.readline(), mfcountmax, False)
            mfcountmin = min(mass_fluct_counts[mass_fluct_counts>0].min(), mfcountmin)


            # replace zeroes with value orders of magnitude smaller
            mass_fluct_counts[mass_fluct_counts==0] = mfcountmin*1e-6


            mass_flucts, mfmax = get_line_to_array(srcfile.readline(), mfmax, False)
            mfmin = min(mass_flucts.min(), mfmin)

            ax.semilogy(mass_flucts, mass_fluct_counts, label=labelnames[f], lw=4, c=colors[f], ls=linestyle[f], alpha=0.7)


        srcfile.close()









 

    #==================================
    # Tweak plots and save figures
    #==================================



    #-------------------------
    # Mass growth and flucts
    #-------------------------

    for ax in axgroup1:
        ax.set_xlim(-1, 1)
        # for inclusive/exclusive
        ax.set_ylim(1e-4, 1e-1)
        #  ax.set_title("Logarithmic Mass Growth")
        ax.set_xlim(-1.05, 1.05)
        xticks = np.linspace(-1,1, 6)
        ax.set_xticks(xticks)
        ax.set_xticklabels(["%.2f" % i for i in xticks])
        ax.set_xlabel(r'$\beta_M = \frac{2}{\pi}\arctan\frac{d \log M}{d \log t}$')
        ylabelmin = int(np.log10(mgcountmin))
        ylabelmax = int(np.log10(mgcountmax)+0.5)
        ax.set_ylim(mgcountmin/2, mgcountmax*2)
        ax.set_yticks([10.0**i for i in range(ylabelmin, ylabelmax, 1)])
        ax.set_yticklabels([r"$10^{"+str(i)+"}$" for i in range(ylabelmin, ylabelmax, 1)])
        ax.set_ylabel(r'$N/N_{A, tot}$')
        ax.grid()
        ax.legend(loc='lower center', ncol=ncols,prop=fontP)
    ax.legend(loc='upper center', ncol=ncols,prop=fontP)

    


    for ax in axgroup2:
        ax.set_xlim(-1, 1)
        ax.set_ylim(1e-4, 1e-1)


        ax.set_xlim(-1.05, 1.05)
        xticks = np.linspace(-1,1, 6)
        ax.set_xticks(xticks)
        ax.set_xticklabels(["%.2f" % i for i in xticks])
        ax.set_xlabel(r'$\xi_M = \frac{\beta_M(k, k+1) - \beta_M(k-1, k)}{2}$')
        ylabelmin = int(np.log10(mfcountmin))
        ylabelmax = int(np.log10(mfcountmax)+0.5)
        ax.set_ylim(mfcountmin/2, mfcountmax*2)
        ax.set_yticks([10.0**i for i in range(ylabelmin, ylabelmax, 1)])
        ax.set_yticklabels([r"$10^{"+str(i)+"}$" for i in range(ylabelmin, ylabelmax, 1)])
        ax.set_ylabel(r'$N/N_{A, tot}$')
        ax.grid()
        ax.legend(loc='lower center', ncol=ncols,prop=fontP)
    ax.legend(loc='upper center', ncol=ncols,prop=fontP)

    ax1.set_title('A: All Haloes')
    ax2.set_title('B: Main Haloes')
    ax3.set_title('C: Subhaloes')
    ax4.set_title('A: All Haloes')
    ax5.set_title('B: Main Haloes')
    ax6.set_title('C: Subhaloes')


    # draw 90% lines
    #  ymin, ymax = ax9.get_ylim()
    #  xl, xr = get_99_percent_position(mass_growth_counts, mass_growth)
    #  ax9.plot([xl,xl], (ymin, ymax), 'r')
    #  ax9.plot([xr,xr], (ymin, ymax), 'r')
    #
    #
    #  ymin, ymax = ax10.get_ylim()
    #  xl, xr = get_99_percent_position(mass_fluct_counts, mass_flucts)
    #  ax10.plot([xl,xl], (ymin, ymax), 'r')
    #  ax10.plot([xr,xr], (ymin, ymax), 'r')

    #  plt.sca(ax9)
    #  plt.xscale('symlog')


    figname="mass_growth.pdf"
    save_figures(figname, fig1)

    figname="mass_fluctuations.pdf"
    save_figures(figname, fig2)




if __name__ == "__main__":
    main()
