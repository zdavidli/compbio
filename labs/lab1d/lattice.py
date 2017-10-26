#!/usr/bin/env python
"""
Program LATTICE

Creates a box of particles on a lattice.

This code is mostly translated from the similar FORTRAN program 
by Daan Frenkel and Berend Smit as described in their book,
Understanding Molecular Simulations: From Algorithms to Applications.

Output of configuration coordinates is in init.pdb

Written by Patrick Fleming, Ph.D., 2003
The Johns Hopkins University 
Department of Biophysics
201 Jenkins Hall 
3400 North Charles Street 
Baltimore, MD 21218 
 
Tel: 410-516-7654
Fax: 410-516-4118
"""
USAGE = """python lattice.py """

# Three subroutines are specific for this program
# The first subroutine, param, creates a dictionary of parameter values
# called Data

def param():

    # Put number of particles, volume density and van der Waals radius
    # in a dictionary called Data

    Data = {}
    Data['NPART'] = 144
    Data['RHOV'] = 0.434
    Data['RAD'] = 1.7

    # here's how to extract them later in the code
    #
    # Npart = int(Data['NPART'])
    # RhoV = float(Data['RHOV'])
    # Rad = float(Data['RAD'])

    return Data

# The second subroutine, lattice, actually creates the lattice of particles
# It requires the dictionary Data to be passed to it

def lattice(Data):
    # place 'Npart' particles on a simple cubic
    # lattice with volume density 'RhoV'

    # First extract parameter values from dictionary, Data

    # Number of particles
    Npart = int(Data['NPART'])

    # Volume density
    RhoV = float(Data['RHOV'])

    # van der Waals radius
    Rad = float(Data['RAD'])

    # Calculate radius cubed
    Rad3 = Rad*Rad*Rad

    # Calculate the volume of a particle
    from math import pi
    Volpart = 4.0*pi*Rad3/3.0

    # Calculate the box size needed to give requested density
    # Box is the length of the cubic box side
    Box = ((Npart*Volpart)/RhoV)**(1.0/3.0)

    # The cube root of the number of particles gives us the number along a side
    # Round up to nearest integer
    n=int(Npart**(1.0/3.0)) + 1

    # The distance between particles, delt, is the (length of the box/number per side)
    delt = round(float(Box))/float(n)

    # Initialize a counter
    count = 0

    # Create empty lists for coordinates
    X=[]
    Y=[]
    Z=[]

    # Initialize so after first iteration coordinates are 0.0, 0.0, 0.0
    dx = -delt

    # Iterate over delt to calculate coordinates and append them to lists
    for i in range(0,n):
        dx = dx + delt
        dy = -delt
        for j in range(0,n):
            dy = dy + delt
            dz = -delt
            for k in range(0,n):
                dz = dz + delt
                if count < Npart:
                    count = count + 1
                    X.append(dx)
                    Y.append(dy)
                    Z.append(dz)

    # Retrun coordinates to calling module
    return X,Y,Z

# The third subroutine writes out the particle coordinates in PDB format.
# It requires a file name and the x,y,z coordinates

def pdbout(file,X,Y,Z):
    # write coordinates of single configuration in PDB format 
    
    # The line label, atom name, residue name, occupancy, and B-factor are 
    # hard wired. The atom name used is 'O' for oxygen but others could be used
    record = 'ATOM  '
    atmnam = 'O'
    resnam = 'UNK'
    occ = 1.0
    bfac = 0.0

    # Iterate over the length of the coordinate lists to write out each line
    # Use length of list X for iteration number (could use Y or Z also)
    # Next line describes PDB format
    format = '%6s%5d  %-3s%4s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
    for i in range(len(X)):
        file.write(format % (record,i+1,atmnam,resnam,i+1,X[i], Y[i], Z[i],occ,bfac))

    # Write the word END as the last line
    file.write('END\n')

# Now define the main program which calls the subroutines

def run():
    # Call param to put parameter values in dictionary "Data"
    Data= param()

    # Call lattice to obtain x,y,z coordinates of particles on lattice
    X,Y,Z = lattice(Data)

    # Open a file and write the coordinates to the file
    ifile = open('init.pdb', 'w')
    pdbout(ifile,X,Y,Z)
    ifile.close()

# The next line ensures that the run subroutine will execute only if
# the script runs as a main script. It allows the program to
# be put in a module itself so that it may be imported into another module
# and the subroutines used there.

if __name__ == "__main__":
    # Finally run the program
    run()

