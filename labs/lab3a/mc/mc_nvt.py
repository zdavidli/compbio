#!/usr/bin/env python
"""
Program MC_NVT

Performs Monte Carlo simulation of Lennard-Jones particles
in the NVT ensemble.

This code is mostly translated from the similar FORTRAN program 
by Daan Frenkel and Berend Smit as described in their book,
Understanding Molecular Simulations: From Algorithms to Applications.

The program here is simplified - no tail corrections

Reduced units are used internally

It requires a box of particles (see lattice.py)

It requires the file modules.py

It requires a parameter file containing a list of the following
parameters:

     cycles      : number of Monte Carlo cycles during equilibration
     Nsamp      : number of Monte Carlo cycles between two sampling periods
     Nstruc     : number of Monte Carlo cycles between sampling structures
     Dr         : maximum displacement in units of sigma (actual will be random(0,1)*Dr)
     Ndispl     : number of attemps to displace a particle per MC cycle
     Npart      : total number of particles
     Temp       : temperature
     Rho        : Volume density 
     Eps        : epsilon, max energy of particle attraction 
     Sig        : distance between particles at zero attraction 

for argon example:

cycles 1000
nsamp 2
dr  0.1
ndispl 100
npart  144
temp 298
rho  0.7
eps 0.238 #kcal/mol
sig 3.4

Output of final configuration is in final.pdb      
Output of configuration ensemble in PDB format is in ensem.pdb
Output of restart configuration is in restart.pdb
Output of Energy/particle, pressure and density is in epd.dat

Written by Patrick Fleming, Ph.D., 2003
The Johns Hopkins University 
Department of Biophysics
201 Jenkins Hall 
3400 North Charles Street 
Baltimore, MD 21218 
 
Tel: 410-516-7654
Fax: 410-516-4118
"""
USAGE = """python mc_nvt.py params [PDB file]"""

# get needed modules
import sys, os
from math import fmod, pi
from random import random
from modules import pdbout,trajout,ener,eneri,totergmc, \
 readparam,adjust,mcmove,samplemc,MC_restartout

# main program
def run(file1,file2):
    # define some parameters
    Data= readparam(file1)
    Dr = float(Data['DR'])
    Npart = int(Data['NPART'])
    Rho = float(Data['RHO'])
    sigma = float(Data['SIG']) #Angstrom
    eps = 1000.0 * float(Data['EPS']) #cal/mol
    Nstruc = int(Data['NSTRUC'])

    # read configuration from file
    print((' Reading configuration from file %s' % file2))
    if file2 == None:
        raise ValueError(' Configuration file not specified')

    # All atom coordinates go into lists
    X = []
    Y = []
    Z = []

    ff = open(file2, 'r')
    nextline = ff.readline  
    while 1:
        line = nextline()
        if not line:
            break
        line = line.strip()
        line = line.split()
        if line[0].upper() == 'REMARK':
            if line[1] == 'Max.':
                # If this is a restart file, read the best displacement ratio
                Dr = float(line[3])
            else:
                continue
        elif line[0].upper() == 'ATOM':
            # convert to reduced distance units
            x = float(line[5])/sigma
            y = float(line[6])/sigma
            z = float(line[7])/sigma
            
            X.append(x)
            Y.append(y)
            Z.append(z)
    ff.close()

    # Print input data
    Ncycles = int(Data['CYCLES'])
    print((' Number of equilibration cycles             : %10d' % Ncycles))
    Nsamp = int(Data['NSAMP'])
    print((' Sample frequency                           : %10d' % Nsamp)) 
    Ndispl = int(Data['NDISPL'])
    print((' Number attmpt. to move a molecule/cycle    : %10d' % Ndispl))
    print((' Maximum displacement                       : %10.3f' % Dr))
    Npart = int(Data['NPART'])
    print((' Number of particles                        : %10d' % Npart))
    SI_temp = float(Data['TEMP'])
    print((' Temperature                                : %10.2f' % SI_temp))
    Rho = float(Data['RHO'])
    print((' Density                                    : %10.3f' % Rho))

    # calculate total energy
    En,Vir = totergmc(Data,X,Y,Z)
    # En is in reduced units
    SI_en = (eps * En) / 1000.0 #kcal/mol
    print((' Total energy initial configuration: %8.3f kcal/mol' % SI_en))

    # open ensemble file
    tfile = open('ensemMC.pdb', 'w')
    # write first configuration to trajectory
    icycl = 0
    traj_flag = 1
    trajout(Data,tfile,icycl,En,X,Y,Z,traj_flag)
    traj_flag = 0

    # open energy, pressure, density log file
    epdfile = open('epd.dat', 'w')
    epdfile.write('    Cycle     Energy    Pressure     Density\n')
    epdfile.write('             kcal/mol     psi        volume\n')

    # start MC-cycle
    # the following parameters will be needed below
    Ndispl = int(Data['NDISPL']) 
    Ncycles = int(Data['CYCLES'])
    Nsamp = int(Data['NSAMP'])

    # initialize some counters
    attempp = naccp = 0
    Attemp = Nacc = 0

    # intialize the subroutine that adjusts the maximum displacement
    Dr, naccp, attempp = adjust(Data,Attemp,attempp,Nacc,naccp,Dr)

    # start either equilibration or production cycles
    for icycl in range(1,Ncycles+1):
        for imove in range(0,Ndispl):
            # attempt to displace a particle
            # here is the Monte Carlo move attempt
            X,Y,Z,En,Vir,Attemp,Nacc = mcmove(Data,X,Y,Z,En,Vir,Attemp,Nacc,Dr)

        # sample averages and write ensemble
        # write stats every Nsamp cycles
        if fmod(icycl,Nsamp) == 0:
            samplemc(Data,icycl,En,Vir,epdfile)

        # write ensemble file every Nstruc cycles
        # this will display actual vdw radii
        if fmod(icycl,Nstruc) == 0:
            trajout(Data,tfile,icycl,En,X,Y,Z,traj_flag)

        # store restart configuration every Ncycles/5 cycles
        if fmod(icycl,(Ncycles/5)) == 0:
            sfile = open('restart.pdb', 'w')
            MC_restartout(Data,sfile,icycl,En,X,Y,Z,Dr)
            sfile.close()
            Dr, naccp, attempp = adjust(Data,Attemp,attempp,Nacc,naccp,Dr)
            print(('======>> Done %8d out of %8d' % (icycl,Ncycles))) 
            SI_enpot = (eps * En) / 1000.0 #kcal/mol
            print((' Potential energy %12.3f kcal/mol' % SI_enpot))

    # print intermediate move stats
    if Ncycles != 0:
        print((' Number of att. to displ. a part.  : %10d' % attempp))
        percent = 100.0*float(naccp)/float(attempp)
        print((' Success: %10d =  %5.2f Percent\n' % (naccp, percent)))

    # finished, write restart file and pdb format file of final configuration
    sfile = open('restart.pdb', 'w')
    MC_restartout(Data,sfile,icycl,En,X,Y,Z,Dr)    
    sfile.close()

    # calculate total energy
    En,Vir = totergmc(Data,X,Y,Z)
    SI_en = (eps * En) / 1000.0 #kcal/mol
    print((' Total energy final configuration: %12.5f kcal/mol' % SI_en))  

    ofile = open('finalMC.pdb', 'w')     
    # this will display actual vdw radii
    pdbout(ofile,X,Y,Z,Data)
    ofile.close()
    # close ensemble and energy/pressure/density files
    tfile.close()
    epdfile.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(USAGE)  
        sys.exit()

    # first file is parameter file
    # second file is restart file if used
    file1 = sys.argv[1]
    if len(sys.argv) < 3:
        file2 = None
    elif len(sys.argv) == 3:
        file2 = sys.argv[2]
    run(file1, file2)

