#!/usr/bin/env python
"""
Program MD_NVT

Performs Molecular Dynamics simulation of Lennard-Jones particles
in the NVT ensemble.

This code is mostly translated from the similar FORTRAN program 
by Daan Frenkel and Berend Smit as described in their book,
Understanding Molecular Simulations: From Algorithms to Applications.

The program here is simplified - no tail corrections, no shift 

It requires either a PDB file or a restart.cv file
 from a previous simulation.

Reduced units are used internally..

It requires a parameter file containing a list of the following
parameters:

#    Delt       : time step MD simulation (in units of ~ps)
#    Tmax       : total simulation time (in units of ~ps)
#    Nsamp      : number of steps between two sampling periods
#    Nstruc     : number of steps between ensemble structures
#    Npart      : total number of particles
#    Temp       : temperature (K)
#    Rho        : density (volume density)
#    Eps        : epsilon, max energy of particle attraction (kcal/mol)
#    Sig        : distance between particles at maximum attraction (Angstrom)
#    Mass       : mass of atom (kg/mol)
#    Equil      : flag for equilibration (1) or production (0)

for Argon as example:

delt 0.001
tmax  10.0
nsamp 20
nstruc 50
npart  144
temp 150
rho  0.6
sig 3.4
eps 0.238
mass 0.03994
equil 1

Output of final configuration is in finalMD.pdb      
Output of configuration trajectory in PDB format is in trajMD.pdb
Output of restart configuration is in restart.cv
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
USAGE = """python md_nve.py md_params [restart configuration file]"""

# get needed modules
import sys, os
from math import fmod, pi
from random import random
from modules import pdbout,trajout,ener,eneri,totergmd, \
 readparam,samplemd,setvel,force,solve,velocs,storemd

# main program
def run(file1,file2):
    # define some parameters
    Data= readparam(file1)
    Npart = int(Data['NPART'])
    Rho = float(Data['RHO'])
    sigma = float(Data['SIG']) #Angstrom
    k_eps = float(Data['EPS']) #kcal/mol
    eps = 1000.0 * k_eps  #cal/mol
    Nstruc = int(Data['NSTRUC']) 
    kb = 1.9872 # cal/mol/K
    equil = int(Data['EQUIL'])

       
    # boxf       : Box length old configuration
    # npartf     : number of particles (over rules Npart in parameter file)
    # time       : ps of previously accumulated simulation
    #              indication of how equilibrated the system it)
    # X[0],Y[0],Z[0],VX[0].VY[0],VZ[0]: position & velocities,particle 0 
    # ...                             : continues for all particles
    # X[n],Y[n],Z[n],VX[n].VY[n],VZ[n]: position & velocities,particle n 
        
    if file2 == None:
        raise ValueError('Configuration file not specified, try "restart.cv"')
    elif file2 == 'restart.cv':
            
        print((' Reading MD positions and velocities from file %s' % file2))
    
        # All atom coordinates and velocities go into lists
        X = []
        Y = []
        Z = []
        VX = []
        VY = []
        VZ = []
    
    
        ff = open(file2, 'r')
        nextline = ff.readline
        line =  nextline()
        line = line.strip()
        line = line.split()
    
        if line[0].upper() == 'REMARK':
            if line[1] == 'MD':
                # This is a restart file,
                # read next three lines
                for i in range(3):
                    line = nextline()
                    if not line:
                        break
                    line = line.strip()
                    if line[:1] == '#':
                        continue
                    line = line.split()
                    if i == 0:
                        # Size of box
                        boxf = float(line[0])
                    elif i == 1: 
                        # Number of atoms
                        npartf = int(line[0])
                    elif i == 2: 
                        # Accumulated simulation time
                        Time = int(line[0])
            
                rho = npartf/boxf**3
            
                sumv2 = 0.0
                sumvx = 0.0
                sumvy = 0.0
                sumvz = 0.0
             
                # read rest of lines
                while 1:                                                        
                    line = nextline()                                           
                    if not line:                                                
                        break                                                   
                    line = line.strip()                                   
                    if line[:1] == '#':                                         
                        continue                                                
                    line = line.split()
                    x = float(line[0])                                          
                    y = float(line[1])                                          
                    z = float(line[2])                                          
                    p = float(line[3])                                          
                    q = float(line[4])                                          
                    r = float(line[5])                                          
                                                                                    
                    X.append(x)                                      
                    Y.append(y)                                      
                    Z.append(z)                                      
            
                    VX.append(p)
                    VY.append(q)
                    VZ.append(r)
            
                    sumv2 = sumv2 + p*p + q*q + r*r
                ff.close()                                                      
                Tempact = sumv2/(3.0*Npart)
                SI_tempact = Tempact * (eps/kb)
    
        elif file2 == 'init.pdb':  # This is a PDB file
            print('Reading PDB file to start MD')
            while 1:
                line = nextline()
                if not line:
                    break
                line = line.strip()
                line = line.split()
                if line[0].upper() == 'REMARK':
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

            # set random velocities
            VX,VY,VZ,Tempact = setvel(Data)
            SI_tempact = Tempact * (eps/kb)

            # Initialize simulation time
            Time = 0.0

                                                                        
    # Print input data
    Delt = float(Data['DELT'])
    print((' Time Step                                  : %10.3f' % Delt))
    Tmax = float(Data['TMAX'])
    print((' Total Simulation Time (ps)                 : %10.2f' % Tmax)) 
    print((' Starting with %8.2f (ps) accumulated simulation time' % Time)) 
    Nsamp = int(Data['NSAMP'])
    print((' Statistical Sample Every                     %10d Steps' % Nsamp))
    Nstruc = int(Data['NSTRUC'])
    print((' Trajectory Output Every                      %10d Steps' % Nstruc)) 
    Npart = int(Data['NPART'])
    print((' Number of particles                        : %10d' % Npart))
    Temp = float(Data['TEMP'])
    print((' Requested Temperature                      : %10.2f' % Temp))
    Rho = float(Data['RHO'])
    print((' Density                                    : %10.3f' % Rho))

    # calculate total energy
    Entot, Enkin, Vir = totergmd(Data,X,Y,Z,VX,VY,VZ)
    Enpot = Entot - Enkin
    SI_entot = (k_eps * Entot) #kcal/mol
    SI_enpot = (k_eps * Enpot) #kcal/mol
    print(' ')
    print((' Total energy initial configuration: %12.3f ' % SI_entot))
    print((' Potential energy initial configuration: %12.3f ' % SI_enpot))
    print(' Initial Temperature = %8.2f\n' % SI_tempact)

    # open trajectory file
    tfile = open('trajMD.pdb', 'w')
    # set flag so first trajectory configuration has box size
    traj_flag = 1

    # open pressure, temperature log file
    epdfile = open('epd.dat', 'w')
    epdfile.write('    Step      Energy    Pressure     Density\n')
    epdfile.write('             kcal/mol     psi        volume\n')


    # start MD
    # the following parameters will be needed below
    stepf = 0
    Nstep = int(Tmax/Delt)
    Nstep10 = Nstep/10
    if Nstep == 0:
        Nstep10 = 0

    # initialize some counters
    Step = stepf

    # create list for forces
    Fx = []
    Fy = []
    Fz = []

    # start doing MD steps
    Start_time = 0
    while Start_time < Tmax:
        # calculate forces
        Fx,Fy,Fz,Enpot,Vir = force(Data,X,Y,Z)
        SI_enpot = (k_eps * Enpot) #kcal/mol

        # solve Newton's equations, update positions and velocities
        X,Y,Z,VX,VY,VZ,Enkin = solve(Data,X,Y,Z,VX,VY,VZ,Fx,Fy,Fz,Delt)

        # update energy and time
        Entot = Enpot + Enkin
        Start_time = Start_time + Delt
        Time = Time + Delt
        Step = Step + 1

        # scale velocities every 20 steps if this is equilibration
        if equil:
            if fmod(Step,20) == 0:
                VX,VY,VZ = velocs(Data,VX,VY,VZ)

        # write stats every 25 steps
        if fmod(Step,25) == 0:
            samplemd(Data,Step,Entot,Vir,Enkin,epdfile)

        # write trajectory file 
        if fmod(Step,Nstruc) == 0:
            trajout(Data,tfile,Step,Enpot,X,Y,Z,traj_flag)
            traj_flag = 0

        # store restart configuration every 10 percent of simulation
        if fmod(Step,Nstep10) == 0:
            sfile = open('restart.cv', 'w')
            storemd(Data,X,Y,Z,VX,VY,VZ,sfile,Time)
            sfile.close()
            Entot, Enkin, Vir = totergmd(Data,X,Y,Z,VX,VY,VZ)
            SI_entot = (k_eps * Entot) #kcal/mol
            print((' Total energy at time %10.2f: %12.3f kcal/mol' % (Time,SI_entot)))
            Enpot = Entot - Enkin
            SI_enpot = (k_eps * Enpot) #kcal/mol
            print((' Potential energy: %12.3f kcal/mol' % SI_enpot))
            Tempact = (2.0 * Enkin)/(3.0*Npart)
            SI_tempact = Tempact * (eps/kb)
            print((' Temperature: %12.3f K\n' % SI_tempact))

    # done
    # close trajectory and pressure files
    tfile.close()
    epdfile.close()

    # calculate total energy
    Entot, Enkin, Vir = totergmd(Data,X,Y,Z,VX,VY,VZ)
    SI_entot = (eps * Entot) / 1000.0 #kcal/mol
    print((' Total energy final configuration: %12.3f kcal/mol' % SI_entot))
    Enpot = Entot - Enkin
    SI_enpot = (eps * Enpot) / 1000.0 #kcal/mol
    print((' Potential energy final configuration: %12.3f ' % SI_enpot))

    # write coordinates of final configuration
    ofile = open('finalMD.pdb', 'w')     
    pdbout(ofile,X,Y,Z,Data)
    ofile.close()


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

