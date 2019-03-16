#!/usr/bin/env python3


#------------------------------
# I/O related functions
#------------------------------

import numpy as np
import warnings
import gc
import os


np.warnings.filterwarnings('ignore')


#==========================
def get_output_info(p, c):
#==========================
    """
    Read in the output info based on the files in the current
    working directory.
    Reads in last directory, ncpu, noutputs. 

        p: params object
        c: constants object
    """

    from os import getcwd
    from os import listdir

    p.workdir = getcwd()
    filelist = listdir(p.workdir)
    
    outputlist = []
    for filename in filelist:
        if "output_" in filename:
            outputlist.append(filename)

    
    if len(outputlist)<1:
        print( "I didn't find any output_XXXXX directories in current working directory.")
        print( "Are you in the correct workdir?")
        quit()

    outputlist.sort()
    p.lastdir = outputlist[-1]
    p.lastdirnr=int(p.lastdir[-5:])
    p.noutput = len(outputlist)

    # read ncpu from infofile in last output directory
    infofile = p.lastdir+'/'+'info_'+p.lastdir[-5:]+'.txt'
    f = open(infofile, 'r')
    ncpuline = f.readline()
    line = ncpuline.split()
    p.ncpu = int(line[-1])


    # read H0
    for i in range(9):
        f.readline()
    H0line = f.readline()
    line = H0line.split()
    c.H0 = float(line[-1])

    #  read omegas
    om = [c.omega_m, c.omega_l, c.omega_k, c.omega_b]
    for o in om:
        fullline = f.readline()
        line = fullline.split()
        o = float(line[-1])

    f.close()

    return








#============================================
def read_mergertree_data(p, sd, mtd, c):
#============================================
    """
    reads in mergertree data as written by the mergertree patch.
    NOTE: requires mergertree.txtYYYYY files.

        p:      params object
        sd:     snapshotdata object
        mtd:    halodata object
        c:      constants object
    """ 


    print( "Reading in mergertree data.")




    #-------------------
    # Preparation
    #-------------------


    # generate list of outputdirnumbers
    startnr=p.lastdirnr
    p.outputnrs = np.array(range(startnr, startnr-p.noutput, -1))

    # define new datatype for mergertree output
    mtree = np.dtype([  ('clump', 'i4'), 
                        ('prog', 'i4'), 
                        ('prog_outnr', 'i4'),
                        ('mass', 'f8'), 
                        ('npart', 'f8'),
                        ('x', 'f8'), 
                        ('y', 'f8'), 
                        ('z', 'f8'),
                        ('vx', 'f8'), 
                        ('vy', 'f8'), 
                        ('vz', 'f8')
                        ])




    #---------------------------
    # Loop over directories
    #---------------------------

    outcount = 0
    dir_template = 'output_'

    for output in range(p.noutput):
        # loop through every output: Progenitor data only starts at output_00002,
        # but you'll need sd.time/sd.redshift data from output_00001!

        # Start with last directory (e.g. output_00060),
        # work your way to first directory (e.g. output_00001)
        dirnr =  str(startnr - output).zfill(5)
        srcdir = dir_template + dirnr


        # Stop early if you reach a directory that has no mergertree.txt* files
        # (Can happen if there are no halos in the simulation yet)


        #------------------------------
        # Read in mergertree data 
        #------------------------------

        fnames = [srcdir + '/' + "mergertree_"+dirnr+'.txt' + str(cpu+1).zfill(5) for cpu in range(p.ncpu)]

        datalist = [np.zeros((1,3)) for i in range(p.ncpu)]
        i = 0
        nofile = 0
        for f in fnames:
            if os.path.exists(f):
                datalist[i] = np.atleast_1d(np.genfromtxt(f, dtype=mtree, skip_header=1))
                i += 1
            else:
                nofile += 1

        if nofile == p.ncpu:
            print( "Didn't find any mergertree data in", srcdir)


        #---------------------------------
        # Sort out data
        #---------------------------------
        if i > 0:
            fulldata = np.concatenate(datalist[:i], axis=0)

            mtd.descendants[output] = fulldata[:]['clump']
            mtd.progenitors[output] = fulldata[:]['prog']
            mtd.progenitor_outputnrs[output] = fulldata[:]['prog_outnr'] 
            mtd.mass[output] = fulldata[:]['mass']
            mtd.npart[output] = fulldata[:]['npart']
            mtd.x[output] = fulldata[:]['x']
            mtd.y[output] = fulldata[:]['y']
            mtd.z[output] = fulldata[:]['z']
            mtd.vx[output] = fulldata[:]['vx']
            mtd.vy[output] = fulldata[:]['vy']
            mtd.vz[output] = fulldata[:]['vz']

            outcount += 1





        
        try:
            #------------------------------------------------------
            # get time, redshift, and units even for output_00001
            #------------------------------------------------------
            fileloc = srcdir+'/info_'+dirnr+'.txt'
            infofile = open(fileloc)
            for i in range(9):
                infofile.readline() # skip first 9 lines

            # get expansion factor
            aline = infofile.readline()
            astring, equal, aval = aline.partition("=")
            afloat = float(aval)
            sd.aexp[output] = afloat

            for i in range(5):
                infofile.readline() # skip 5 lines

            # get unit_l
            unitline = infofile.readline()
            unitstring, equal, unitval = unitline.partition("=")
            unitfloat = float(unitval)
            sd.unit_l[output] = unitfloat
            
            # get unit_dens
            unitline = infofile.readline()
            unitstring, equal, unitval = unitline.partition("=")
            unitfloat = float(unitval)
            sd.unit_dens[output] = unitfloat

            # get unit_t
            unitline = infofile.readline()
            unitstring, equal, unitval = unitline.partition("=")
            unitfloat = float(unitval)
            sd.unit_t[output] = unitfloat

            infofile.close()
    
        except IOError: # If file doesn't exist
            print( "Didn't find any info data in ", srcdir)
            break


    # keep only entries that contain data
    if outcount > 1:
        mtd.descendants = mtd.descendants[:outcount]
        mtd.progenitors = mtd.progenitors[:outcount]
        mtd.progenitor_outputnrs = mtd.progenitor_outputnrs[:outcount]
        mtd.mass = mtd.mass[:outcount]
        mtd.npart = mtd.npart[:outcount]
        mtd.x = mtd.x[:outcount]
        mtd.y = mtd.y[:outcount]
        mtd.z = mtd.z[:outcount]
        mtd.vx = mtd.vx[:outcount]
        mtd.vy = mtd.vy[:outcount]
        mtd.vz = mtd.vz[:outcount]
        
        p.nout = outcount



    #--------------------------------------
    # Transform units to physical units
    #--------------------------------------

    sd.unit_m = sd.unit_dens*sd.unit_l**3/c.M_Sol
    sd.unit_l /= c.Mpc
    sd.unit_t /= c.Gyr

    # transform units to physical units
    for i in range(len(mtd.descendants)):
        mtd.x[i] *= sd.unit_l[i] # only transform later when needed; Need to check for periodicity first!
        mtd.y[i] *= sd.unit_l[i]
        mtd.z[i] *= sd.unit_l[i]
        mtd.mass[i] *= sd.unit_m[i]
        mtd.vx[i] *= sd.unit_l[i]/sd.unit_t[i]
        mtd.vy[i] *= sd.unit_l[i]/sd.unit_t[i]
        mtd.vz[i] *= sd.unit_l[i]/sd.unit_t[i]


    # collect garbage
    gc.collect()



    return 






#=====================================
def read_halo_data(p, hd):
#=====================================
    """
    reads in clumpfinder data from files.
        p :   params object
        hd :  halodata object
    """


    print("Reading in halo data.")


    # Loop over files
    dir_template = 'output_'
    startnr=p.lastdirnr
    for output in range(p.noutput):

        dirnr =  str(startnr - output).zfill(5)
        srcdir = dir_template + dirnr

        fnames = [srcdir + '/' + "halo_"+dirnr+'.txt' + str(cpu+1).zfill(5) for cpu in range(p.ncpu)]

        try:
            datalist = [np.zeros((1)) for i in range(p.ncpu)]
            i = 0
            for f in fnames:
                datalist[i] = np.atleast_1d(np.loadtxt(f, dtype='int', skiprows=1, usecols=[0]))
                i += 1

            if i > 0:
                fulldata = np.concatenate(datalist[:i], axis=0)
                hd.haloid[output] = fulldata[:]

        except IOError: # If file doesn't exist
            print( "Didn't find any halo data in", srcdir)

    

    return


