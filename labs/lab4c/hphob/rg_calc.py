#!/usr/bin/env python
"""
Program RG_CALC

Calculates radius of gyration from PDB file or LINUS simulation trajectory file.

Writes out to file rg.dat column of Rg values for each structure in PDB file.

Written by Patrick Fleming, Ph.D., 2012                                                      
The Johns Hopkins University                                                                 
Department of Biophysics                                                                     
103 Jenkins Hall                                                                             
3400 North Charles Street                                                                    
Baltimore, MD 21218                                                                          
                                                                                             
"""                                                                                          
                                                                                             
USAGE = """python3 rg_calc.py traj.pdb"""                           
                                                                                             
import sys
import math

def get_config(ff):
    # read input trajectory file in PDB format  
    # this file can have multiple configurations with an "END" or "ENDMDL"
    # at the end of each one
    # Example:

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
            return 0,0.0,0.0,0.0
        line = line.strip()
        line = line.split()
        if line[0].upper() == 'REMARK':
            continue
        if line[0].upper() == 'COMPND':
            continue
        if line[0].upper() == 'SEQRES':
            continue
        elif line[0].upper() == 'ATOM':
            X.append(float(line[5]))
            Y.append(float(line[6]))
            Z.append(float(line[7]))
        elif (line[0].upper() == 'ENDMDL') or (line[0].upper() == 'END'):
            # first item is a flag, 1 = data from current configuration
            return 1,X,Y,Z

def radius_of_gyration(X,Y,Z):
    #
    # Rg = sqrt(sum((r - <r>)**2)/N)
    #
    # Initialize
    rg = x = y = z = 0.0
    numatoms = len(X)
    # Calculate center of geometry
    for i in range(numatoms):
        x = x + X[i]
        y = y + Y[i]
        z = z + Z[i]
    n = 1.0/numatoms
    x = x*n
    y = y*n
    z = z*n
    for i in range(numatoms):
        rg = rg + (X[i] - x)**2 + (Y[i] - y)**2 + (Z[i] - z)**2
    rg = math.sqrt(rg * n)
    return rg

def run(file):
    # open rg output file
    rgfile = open('rg.dat', 'w')
    
    # open PDB file for reading coordinates
    ff = open(file, 'r')

    # need a counter for number of frames
    n = 0
 
    # start main loop
    while 1:
        n = n + 1
        # read configuration from PDB file
        # ok is a flag: 0=finished, 1=have data
        ok,X,Y,Z = get_config(ff)

        # if no more data exit
        if ok == 0:
            break

        rg = radius_of_gyration(X,Y,Z)
        rgfile.write('%8.3f \n' % rg)

    # write rg file
    print(' Writing rg data file ',rgfile.name)

    # close files
    ff.close()
    rgfile.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(USAGE)
        sys.exit()

    file = sys.argv[1]
    run(file)
