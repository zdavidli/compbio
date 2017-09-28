#!/usr/bin/env python
"""
Program RDF_CALC

Calculates radial distribution function from simulation trajectory
file.

Writes out g(r) as a function of pairwise distance for plotting.

Written by Patrick Fleming, Ph.D., 2003                                                      
The Johns Hopkins University                                                                 
Department of Biophysics                                                                     
201 Jenkins Hall                                                                             
3400 North Charles Street                                                                    
Baltimore, MD 21218                                                                          
                                                                                             
"""                                                                                          
                                                                                             
USAGE = """python rdf_calc.py traj.pdb"""                           
                                                                                             
import sys

def get_config(ff, box, npart):
    # read input trajectory file in PDB format  
    # this file will have multiple configurations with a "TER" at the end
    # of each one
    # Example:

    #REMARK                  box Length   Num. Particles    
    #REMARK Restart Data      18.9715455          144        
    #MODEL
    #ATOM      0  O   UNK     0      16.488  18.480  18.805  1.00  0.00
    #ATOM      1  O   UNK     1       0.102  18.632   3.668  1.00  0.00
    #.
    #.
    #.
    #ATOM    142  O   UNK   142       9.824  15.934  13.027  1.00  0.00
    #ATOM    143  O   UNK   143      10.316  15.643  16.447  1.00  0.00
    #ENDMDL
    #MODEL
    #ATOM      0  O   UNK     0      16.488  18.480  18.805  1.00  0.00
    #ATOM      1  O   UNK     1       0.102  18.632   3.668  1.00  0.00
    #etc.

    # the data read from trajectory file are as follows:

    # box           : box length configuration
    # npart         : number of particles 
    # X[0],Y[0],Z[0]: position particle 0
    # ...           : continues for all particles
    # X[n],Y[n],Z[n]: position particle n

    # Initialize lists for coordinates
    X = []
    Y = []
    Z = []

    nextline = ff.readline
    while 1:
        line = nextline()
        if not line:
            # first item is a flag, zero = end of file
            # the rest are dummies
            return 0,0.0,0,0.0,0.0,0.0
        line = line.strip()
        line = line.split()
        if line[0].upper() == 'REMARK':
            if line[1] == 'Restart':
                # Box length was in reduced units - change for argon
                box =  float(line[3])
                npart = int(line[4])
                print(box, npart)
            else:
                continue
        elif line[0].upper() == 'ATOM':
            X.append(float(line[5]))
            Y.append(float(line[6]))
            Z.append(float(line[7]))
        elif line[0].upper() == 'ENDMDL':
            # first item is a flag, 1 = data from current configuration
            return 1,box,npart,X,Y,Z

def rdf_calc(box,npart,X,Y,Z,ngr,delgi,g):
    # accumulates number of particles at distance (r*delgi)
    # for radial distribution function

    # box = length of box side
    # npart = number of particles
    # X,Y,Z = particle coordinate lists
    # ngr = number of configurations calculated (for normalization)
    # delgi = radial shell width (bin increment)
    # g = list of running sums for each shell (bin)

    from math import sqrt

    hbox = box/2.0
    ngr = ngr + 1
    for i in range(0,npart):
        xi = X[i]
        yi = Y[i]
        zi = Z[i]

        # look only at i+1 other particles
        for j in range(i+1,npart):
            dx = xi - X[j]
            dy = yi - Y[j]
            dz = zi - Z[j]

            # periodic boundary conditions
            if dx > hbox:
                dx = dx - box
            elif dx < -hbox:
                dx = dx + box
            if dy > hbox:
                dy = dy - box
            elif dy < -hbox:
                dy = dy + box
            if dz > hbox:
                dz = dz - box
            elif dz < -hbox:
                dz = dz + box
        
            # use square of distance rather than calculating square root
            r2 = dx*dx + dy*dy + dz*dz 

            # limit accumulation to only atoms within half box distance
            # farther away should be same as bulk
            if r2 <= hbox*hbox:
                r = sqrt(r2)
                ig = int(r*delgi)

                # finally update distance bin
                g[ig] = g[ig] + 1


    return ngr,g

def rdf_write(box,npart,ngr,g,rdffile):
    # writes radial distribution data to file

    from math import pi

    hbox = box/2.0
    nhgr = len(g)
    delg = hbox/nhgr
    rhon = npart/(box**3)
    for i in range(1,nhgr):
        # r = distance at center of radial shell
        r = (i + 0.5)*delg

        # dr = volume of shell
        dr = 4.0*pi*(delg**3)*((i+1)**3 - i**3)/3.0

        # g[i] = rdf normalized for number of configurations included
        # and volume of shell, and expected rdf given particle density
        g[i]=2.0*g[i]/(ngr*dr*rhon*npart)
 
        # write it out
        rdffile.write('%10.3f %10.3f\n' % (r,g[i]))

def run(file):
    # open rdf output file
    rdffile = open('rdf.dat', 'w')
    
    # open trajectory file for reading coordinates
    ff = open(file, 'r')

    # initialize lists for coordinates
    X = []
    Y = []
    Z = []

    # need a counter to know when on first file
    n = 0

    # initialize until read from file
    box = 0.0
    npart = 0

    # start main loop
    while 1:
        n = n + 1
        # read configuration from trajectory file
        # ok is a flag: 0=finished, 1=have data
        ok,box,npart,X,Y,Z = get_config(ff,box,npart)

        # if first time do some housekeeping
        if n == 1:
            # save box and npart for later calculation
            box_sav = box
            npart_sav = npart
            # initialize some parameters for rdf calculation
            ngr = 0
            hbox = box/2.0
            nhgr = 250
            delg = hbox/nhgr
            delgi = 1/delg
            # initialize the rdf bins with zeros
            g=[]
            for i in range(0,nhgr):
               g.append(0.0)

        # if no more data exit
        if ok == 0:
            break

        # accumulate data from each configuration in trajectory file
        # for g calculation
        ngr,g = rdf_calc(box,npart,X,Y,Z,ngr,delgi,g)
        print(' Read configuration ', n)

    # write rdf file
    print(' Writing rdf file ',rdffile.name)
    rdf_write(box_sav,npart_sav,ngr,g,rdffile)

    # close files
    ff.close()
    rdffile.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(USAGE)
        sys.exit()

    file = sys.argv[1]
    run(file)
