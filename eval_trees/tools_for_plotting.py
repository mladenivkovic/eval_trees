#!/usr/bin/env python3

#----------------------------------------------------------
# Some tools for evalplot.py and evaljumpers.py that they
# share.
#----------------------------------------------------------


import pickle


class plot_selection():

    # meant for "fieldtest" run debugging
    do_debug = False

    # using sussing criteria
    do_sussing_ntrace = True
    do_sussing_incexcl = False

    # use my arbitrary mass threshold
    do_my_threshold_ntrace = False
    do_my_threshold_incexcl = False

    # use minimum particle threshold
    npart_threshold = 200
    do_npart_debug = False
    do_npart_threshold_ntrace = False
    do_npart_threshold_incexcl = False


    def __init__(self):
        return




#================================
def file_selection(ps):
#================================
    """
    Select file and plot parameters

    Parameters:

        ps: plot_selection object


    Returns:
        allfiles: list of all file names to be read in
        labelnames: list of a legend label for each file
        suffix: suffix for pngs for selected case
        linestyle: list of line styles for each file
    """

    #--------------------
    # Debug
    #--------------------
    if ps.do_debug:
        allfiles = ['eval_trees-no-threshold.pkl']
        #  allfiles = ['eval_trees-sussing-criteria.pkl']
        labelnames = ['label']
        suffix = 'debug'
        linestyle = ['-']


    #--------------------
    # Actually in use
    #--------------------

    # For no-thresholds

    elif ps.do_my_threshold_ntrace:
        #  ntrace = [1, 10, 100, 1000]
        ntrace = [1, 100, 200, 1000]
        #  ntrace = [1, 10, 50, 100, 200, 500, 1000]
        allfiles = ['eval_trees-no-threshold-ntrace-'+str(i)+'.pkl' for i in ntrace]
        labelnames = [r'$n_{mb}='+str(i)+'$' for i in ntrace]
        #  linestyle = ['-', '-', '--', '--', '-.', '-.', ':']
        linestyle = ['-']*10
        suffix='ntrace'

    elif ps.do_my_threshold_incexcl:
        allfiles = [ 'eval_trees-no-threshold-inclusive-nosaddle.pkl',
                    'eval_trees-no-threshold-exclusive-nosaddle.pkl',
                    'eval_trees-no-threshold-inclusive-saddle.pkl',
                    'eval_trees-no-threshold-exclusive-saddle.pkl' ]
        labelnames = [ 'inclusive loosely bound', 'exclusive loosely bound', 'inclusive strictly bound', 'exclusive strictly bound' ]
        linestyle = ['-', '--', '-', '--']
        suffix='inc-excl'


    #  For sussing thresholds
    elif ps.do_sussing_ntrace:
        #  ntrace = [1, 10, 50, 100, 200, 500, 1000]
        ntrace = [100, 1000]
        allfiles = ['eval_trees-sussing-criteria-ntrace-'+str(i)+'.pkl' for i in ntrace]
        labelnames = [r'$n_{mb}='+str(i)+'$' for i in ntrace]
        linestyle = ['-', '-', '--', '--', '-.', '-.', ':']
        suffix='ntrace'


    elif ps.do_sussing_incexcl:
        allfiles = [ 'eval_trees-sussing-criteria-inclusive-nosaddle.pkl',
                    'eval_trees-sussing-criteria-exclusive-nosaddle.pkl',
                    'eval_trees-sussing-criteria-inclusive-saddle.pkl',
                    'eval_trees-sussing-criteria-exclusive-saddle.pkl' ]
        labelnames = [ 'inclusive loosely bound', 'exclusive loosely bound', 'inclusive strictly bound', 'exclusive strictly bound' ]
        #  linestyle = ['--', '--', ':', ':']
        linestyle = ['-']*10
        suffix='inc-excl'

    elif ps.do_npart_debug:
        allfiles = [ 'eval_trees-npart-threshold-'+str(ps.npart_threshold)+'.pkl']
        labelnames = [ 'debug npart threshold' ]
        linestyle = ['-']*10
        suffix='debug-npart'

    elif ps.do_npart_threshold_ntrace:
        ntrace = [1, 10, 100, 1000]
        #  ntrace = [1, 100, 200, 1000]
        #  ntrace = [1, 10, 50, 100, 200, 500, 1000]
        allfiles = ['eval_trees-npart-threshold-'+str(ps.npart_threshold)+'-ntrace-'+str(i)+'.pkl' for i in ntrace]
        labelnames = [r'$n_{mb}='+str(i)+'$' for i in ntrace]
        #  linestyle = ['-', '-', '--', '--', '-.', '-.', ':']
        linestyle = ['-']*10
        suffix='ntrace-npart'

    elif ps.do_npart_threshold_incexcl:
        allfiles = [ 'eval_trees-npart-threshold-'+str(ps.npart_threshold)+'-inclusive-nosaddle.pkl',
                    'eval_trees-npart-threshold-'+str(ps.npart_threshold)+'-exclusive-nosaddle.pkl',
                    'eval_trees-npart-threshold-'+str(ps.npart_threshold)+'-inclusive-saddle.pkl',
                    'eval_trees-npart-threshold-'+str(ps.npart_threshold)+'-exclusive-saddle.pkl' ]
        labelnames = [ 'inclusive loosely bound', 'exclusive loosely bound', 'inclusive strictly bound', 'exclusive strictly bound' ]
        linestyle = ['-', '--', '-', '--']
        suffix='inc-excl-npart'



    else:
        print("You need to select for which dataset I'm plotting")
        quit()


    return allfiles, labelnames, suffix, linestyle




#===============================================
def print_evaltrees_parameters(filename):
#===============================================
    """
    Read in an evaltrees.py output file and print the
    interesting parameters that were used.
    """


    srcfile = open(filename, 'rb')
    p, sd, r = pickle.load(srcfile) 
    srcfile.close()


    print("--------------------------------------")
    print("evaltrees.py parameters:")
    print("    mth_main", p.mth_main)
    print("    mth_sub", p.mth_sub)
    print("    sussing", p.sussing)
    try:
        print("    use_npart_threshold", p.use_npart_threshold)
        print("    npart_thresh_sub", p.npart_thresh_sub)
        print("    npart_thresh_main", p.npart_thresh_main)
    except AttributeError:
        print("params object has missing attributes. re-run evaltree.py")
        quit()
    print("--------------------------------------")


