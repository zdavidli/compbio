"""
Modules for simulation programs

Written by Patrick Fleming, Ph.D., 2003

The Johns Hopkins University 
Department of Biophysics
201 Jenkins Hall 
3400 North Charles Street 
Baltimore, MD 21218 
 
"""

def readparam(file):
    # read the parameter values - put in dictionary, Data

    from math import pi
    import string
    if type(file) == type(''):
        ff = open(file, 'r')
    else:
        ff = file
    nextline = ff.readline
#   Strip = string.strip
#   Split = string.split
#   Upper = string.upper
#   Int = string.atoi

    Data = {}
    while 1:
        line = nextline()
        if not line:
            break
        line = line.strip()
        line = line.upper()
        if line[:1] == '#':
            continue
        line = line.split()
        if len(line) < 1:
            continue
        param = line[0]
        Data[param]=line[1]

    return Data

def pdbout(file,X,Y,Z,Data):
    # write coordinates of single configuration in PDB format 
    # call all atoms oxygen with unknown residue name
    sigma = float(Data['SIG'])
    record = 'ATOM  '
    atmnam = 'O'
    resnam = 'UNK'
    occ = 1.0
    bfac = 0.0
    for i in range(len(X)):
        x = X[i] * sigma
        y = Y[i] * sigma
        z = Z[i] * sigma
        file.write('%6s%5d  %-3s%4s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n' %
        (record,i,atmnam,resnam,i,x,y,z,occ,bfac))
    file.write('END\n')

def trajout(Data,file,icycl,En,X,Y,Z,traj_flag):
    # write coordinates of multiple configurations in PDB format
    # call all atoms oxygen with unknown residue name
#   k_eps = float(Data['EPS']) #kcal/mol
#   SI_en = k_eps * En #kcal/mol
#   title1 = 'REMARK Cycle %d        Energy %g\n' % (icycl,SI_en)
#   file.write(title1)

    sigma = float(Data['SIG'])
    Npart = int(Data['NPART'])
    Rho = float(Data['RHO'])
    Box = (Npart/Rho)**(1.0/3.0)
    Box = Box * sigma

    if traj_flag:
        file.write('REMARK                  Box Length   Num. Particles\n')
        file.write('REMARK Restart Data   %13.7f     %8d\n' % (Box,Npart))
    file.write('MODEL\n')

    record = 'ATOM  '
    atmnam = 'O'
    resnam = 'UNK'
    occ = 1.0
    bfac = 0.0
    for i in range(len(X)):
        # use argon sigma for display
        x = X[i] * sigma
        y = Y[i] * sigma
        z = Z[i] * sigma
        file.write('%6s%5d  %-3s%4s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n' %
        (record,i,atmnam,resnam,i,x,y,z,occ,bfac))
    file.write('ENDMDL\n')     

def MC_restartout(Data,file,icycl,En,X,Y,Z,Dr):
    # write coordinates of configuration in PDB format
    # call all atoms oxygen with unknown residue name

    title1 = 'REMARK Cycle %d        Energy %g\n' % (icycl,En)
    file.write(title1) 
    Npart = int(Data['NPART'])
    Rho = float(Data['RHO'])

    # hard shell  radius
    sigma = float(Data['SIG'])
    Rad = sigma/2.0

    # Calculate radius cubed
    Rad3 = Rad*Rad*Rad

    # Calculate the volume of a particle
    from math import pi
    Volpart = 4.0*pi*Rad3/3.0

    # Calculate the box size needed to give requested density
    # Box is the length of the cubic box side
    SI_Box = ((Npart*Volpart)/Rho)**(1.0/3.0)


    file.write('REMARK          Box Length       Num. Particles\n')
    file.write('REMARK Data   %13.7f     %8d\n' % (SI_Box,Npart))
    file.write('REMARK Max. Displ. %10.4f\n' % (Dr))
    record = 'ATOM  '
    atmnam = 'O'
    resnam = 'UNK'
    occ = 1.0 
    bfac = 0.0
    for i in range(len(X)):
        # keep reduced units for restart file
        x = X[i] * sigma
        y = Y[i] * sigma
        z = Z[i] * sigma
        file.write('%6s%5d  %-3s%4s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n' %
        (record,i,atmnam,resnam,i,x,y,z,occ,bfac))
    file.write('TER\n')

def ener(Data,R2):
    # calculate energy between two particles
    # units are reduced
    # use half box length as the cut off distance
    Npart = int(Data['NPART'])
    Rho = float(Data['RHO'])
    Box = (Npart/Rho)**(1.0/3.0)
    Hbox = Box/2.0
    Rc = Hbox
    Rc2 = Rc*Rc

    r2i = 1.0/R2
    r6i = r2i*r2i*r2i
    En = 4*(r6i*r6i-r6i)
    Vir =  48*(r6i*r6i-0.5*r6i)
    return En, Vir

def eneri(Data,X,Y,Z,xi,yi,zi,i,Jb):
    #calculates the energy of particle I with particles j=Jb,npart
 
    # xi (input) x coordinate particle I
    # yi (input) y coordinate particle I
    # zi (input) z coordinate particle I
    # i  (input) particle number
    # Jb (input) = 0 calculates energy particle I with all other particle
    #            = Jb calculates energy particle I with all particles j > Jb
    # En  (output) energy particle i
    # Vir (output) virial particle i

    En = 0.0
    Vir = 0.0
    Npart = int(Data['NPART'])
    Rho = float(Data['RHO'])
    Box = (Npart/Rho)**(1.0/3.0)
    Hbox = Box/2.0
    Rc = Hbox
    Rc2 = Rc*Rc
    for j in range(Jb,Npart):
        if j != i:
            dx = xi - X[j]
            dy = yi - Y[j]
            dz = zi - Z[j]

            # Periodic boundary
            if dx > Hbox:
                dx = dx - Box
            else:
                if dx < -Hbox:
                    dx = dx + Box
            if dy > Hbox:
                dy = dy - Box
            else:
                if dy < -Hbox:
                    dy = dy + Box
            if dz > Hbox:
                dz = dz - Box
            else:
                if dz < -Hbox:
                    dz = dz + Box
            R2 = dx*dx + dy*dy + dz*dz
            if R2 <= Rc2:
                enij, virij = ener(Data,R2)
                En = En + enij
                Vir = Vir + virij
    return En,Vir


def totergmc(Data,X,Y,Z):
    # sums energy and virial
    # virial not needed, old code needs to be cleaned up

    Ener = 0.0
    Vir = 0.0
    Npart = int(Data['NPART'])
    for i in range(0,Npart-1):
        xi = X[i]
        yi = Y[i]
        zi = Z[i]
        Jb = i +1
        eni,viri = eneri(Data,X,Y,Z,xi,yi,zi,i,Jb)
        Ener = Ener + eni
        Vir = Vir + viri
    return Ener, Vir


def adjust(Data,Attemp,attempp,Nacc,naccp,Dr):
    # adjusts maximum displacement such that 50% of the
    # moves will be accepted in MC simulation
 
    # Attemp (input) number of attemps that have been performed to displace a particle
    # attempp (input) number of previous attemps that have been performed 
    # Nacc   (input) number of successful attemps to displace a particle
    # naccp  (input) number of previous successful attemps 
    # Dr     (output) new maximum displacement
 
    format1 = ' Max. displ. set to : %6.3f (old : %6.3f)'
    format2 = ' Frac. acc.: %4.2f attempts: %7d success: %7d'
 
    Npart = int(Data['NPART'])
    Rho = float(Data['RHO'])
    Box = (Npart/Rho)**(1.0/3.0)

    Hbox = Box/2.0
    frac = 0.0
    if Attemp == 0 or attempp >= Attemp:
        naccp = Nacc
        attempp = Attemp
    else:
        frac = float(Nacc-naccp)/float(Attemp-attempp)
        dro = Dr
        Dr = Dr*abs(frac/0.5)
        # limit the range
        if (Dr/dro) > 1.5:
            Dr = dro*1.5
        if (Dr/dro) < 0.5:
            Dr = dro*0.5
        if Dr > (Hbox/2.0):
            Dr = Hbox/2.0
        print((format1 % (Dr,dro)))
        print((format2 % (frac, Attemp - attempp, Nacc - naccp)))
        # store nacc and attemp for next use
        naccp = Nacc
        attempp = Attemp
    return Dr, naccp, attempp

def mcmove(Data,X,Y,Z,En,Vir,Attemp,Nacc,Dr):
    # attempts to displace a randomly selected particle
    # applies the Metropolis criterion for acceptance
 
    # En     (input/output) : total energy
    # Vir    (input/output) : total virial
    # Attemp (input/output) number of attempts that have been
    #                       performed to displace a particle
    # Nacc   (input/output) number of successful attempts
    #                       to displace a particle

    from math import exp
    from random import random

    Npart = int(Data['NPART'])
    Rho = float(Data['RHO'])
    Box = (Npart/Rho)**(1.0/3.0)
    SI_temp = float(Data['TEMP'])
    kb = 1.9872 # cal/mol/K
    eps = 1000.0 * float(Data['EPS']) #cal/mol
    Temp = SI_temp/(eps/kb)
    Beta = 1.0/Temp

    Attemp = Attemp + 1
    Jb = 0

    # select a particle at random
    o = int(Npart*random())
  
    # calculate energy old configuration
    xo = X[o]
    yo = Y[o]
    zo = Z[o]
    Eno,Viro = eneri(Data,X,Y,Z,xo,yo,zo,o,Jb) 

    # give particle a random displacement
    xn = X[o] + (random()-0.5)*Dr
    yn = Y[o] + (random()-0.5)*Dr
    zn = Z[o] + (random()-0.5)*Dr

    # calculate energy new configuration
    Enn,Virn = eneri(Data,X,Y,Z,xn,yn,zn,o,Jb)

    # Metropolis acceptance test
    if random() < exp(-Beta*(Enn-Eno)):
        # accepted
        Nacc = Nacc + 1
        En = En + (Enn-Eno)
        Vir = Vir + (Virn-Viro)
        # put particle in simulation box in case random displacement 
        # moved it out of box (periodic boundary conditions)
        # (box origin is 0,0,0 and particle coords are positive)
        if xn < 0.0:
            xn = xn + Box
        if xn > Box:
            xn = xn - Box
        if yn < 0.0:
            yn = yn + Box
        if yn > Box:
            yn = yn - Box
        if zn < 0.0:
            zn = zn + Box
        if zn > Box:
            zn = zn - Box
        X[o] = xn
        Y[o] = yn
        Z[o] = zn

    return X,Y,Z,En,Vir,Attemp,Nacc

def samplemc(Data,icycl,En,Vir,prtfile):
    # write quantities (pressure, energy and density) to file

    from math import pi
    Npart = int(Data['NPART'])
    Rho = float(Data['RHO'])
    Box = (Npart/Rho)**(1.0/3.0)
    sigma = float(Data['SIG'])
    sig3 = sigma * sigma * sigma
    SI_temp = float(Data['TEMP'])
    kb = 1.9872 # cal/mol/K
    eps = 1000.0 * float(Data['EPS']) #cal/mol
    Temp = SI_temp/(eps/kb)
    Beta = 1.0/Temp
    k_eps = float(Data['EPS']) #kcal/mol
    SI_en = (k_eps * En)  #kcal/mol
    eps = 1000.0 * k_eps


    if Npart != 0:
        vol = Box**3
        press = (Npart/vol)/Beta + Vir/(3.0*vol)
        SI_press = press * (sig3/eps) #psi
        rho = Npart/vol
       
    else:
        enp = 0.0
        press = 0.0
        rho = 0.0
    prtfile.write('%9d %10.3f %10.3f %10.3f\n' % (icycl,SI_en,SI_press,rho))

# Modules below for MD

def totergmd(Data,X,Y,Z,VX,VY,VZ):
    # sums energy and virial

    Enpot = 0.0
    Virtot = 0.0
    Enkin = 0.0
    Npart = int(Data['NPART'])
    for i in range(0,Npart-1):
        xi = X[i]
        yi = Y[i]
        zi = Z[i]
        Jb = i +1
        eni,viri = eneri(Data,X,Y,Z,xi,yi,zi,i,Jb)
        Enpot = Enpot + eni
        Virtot = Virtot + viri
        Enkin = Enkin + (VX[i]**2 + VY[i]**2 + VZ[i]**2)

    # correct for double counting each pair
    Enkin = 0.5*Enkin
    Entot = Enpot + Enkin
    return Entot, Enkin, Virtot

def samplemd(Data,Step,Entot,Vir,Enkin,prtfile):
    # write quantities (pressure, energy and density) to file

    from math import pi
    Npart = int(Data['NPART'])
    Rho = float(Data['RHO'])
    Box = (Npart/Rho)**(1.0/3.0)
    sigma = float(Data['SIG'])
    sig3 = sigma * sigma * sigma
    SI_temp = float(Data['TEMP'])
    kb = 1.9872 # cal/mol/K
    eps = 1000.0 * float(Data['EPS']) #cal/mol
    Temp = SI_temp/(eps/kb)
    Beta = 1.0/Temp
    k_eps = float(Data['EPS']) #kcal/mol
    eps = 1000.0 * k_eps

    if Npart != 0:
        Enpot = (Entot-Enkin)
        SI_enpot = (k_eps * Enpot)  #kcal/mol
        vol = Box**3
        press = (Npart/vol)/Beta + Vir/(3.0*vol)
        SI_press = press * (sig3/eps) #psi

    else:
        SI_enpot = 0.0
        SI_press = 0.0
        Rho = 0.0
    prtfile.write('%9d %10.3f %10.3f %10.3f\n' % (Step,SI_enpot,SI_press,Rho))


def force(Data,X,Y,Z):
    # calculate force on atom given pairwise distance to
    # all other atoms
    Enpot = 0.0
    Vir = 0.0
    Npart = int(Data['NPART'])
    Rho = float(Data['RHO'])
    Box = (Npart/Rho)**(1.0/3.0)

    Hbox = Box/2.0
    Rc = Hbox
    Rc2 = Rc*Rc

    # create local lists for force and initialize with zeros
    Fx = []
    Fy = []
    Fz = []
    for i in range(0,Npart):
        Fx.append(0.0)
        Fy.append(0.0)
        Fz.append(0.0)

    for i in range(0,Npart-1):
        xi = X[i]
        yi = Y[i]
        zi = Z[i]

        # look only at i+1 other particles
        for j in range(i+1,Npart):
            dx = xi - X[j]
            dy = yi - Y[j]
            dz = zi - Z[j]

            # periodic boundary conditions
            if dx > Hbox:
                dx = dx - Box
            elif dx < -Hbox:
                dx = dx + Box
            if dy > Hbox:
                dy = dy - Box
            elif dy < -Hbox:
                dy = dy + Box
            if dz > Hbox:
                dz = dz - Box
            elif dz < -Hbox:
                dz = dz + Box
            r2 = dx*dx + dy*dy + dz*dz
            if r2 <= Rc2:
                enij, virij = ener(Data,r2)
                Enpot = Enpot + enij
                Vir = Vir + virij
                fr = virij/r2
                Fx[i] = Fx[i] + fr*dx
                Fy[i] = Fy[i] + fr*dy
                Fz[i] = Fz[i] + fr*dz
                Fx[j] = Fx[j] - fr*dx
                Fy[j] = Fy[j] - fr*dy
                Fz[j] = Fz[j] - fr*dz

    return(Fx,Fy,Fz,Enpot,Vir)

def setvel(Data):
    # assign random velocities to each atom

    from math import sqrt
    from random import random
    Npart = int(Data['NPART'])
    SI_temp = float(Data['TEMP'])
    kb = 1.9872 # cal/mol/K
    eps = 1000.0 * float(Data['EPS']) #cal/mol
    Temp = SI_temp/(eps/kb)

    vx0 = 0.0
    vy0 = 0.0
    vz0 = 0.0
    v2 = 0.0
    VX = []
    VY = []
    VZ = []
    for i in range(0,Npart):
        VX.append(random() - 0.5)
        VY.append(random() - 0.5)
        VZ.append(random() - 0.5)
        vx0 = vx0 + VX[i]
        vy0 = vy0 + VY[i]
        vz0 = vz0 + VZ[i]
        v2 = v2 + (VX[i])**2 + (VY[i])**2 + (VZ[i])**2

    # set center of mass movement to zero
    vx0 = vx0/Npart
    vy0 = vy0/Npart
    vz0 = vz0/Npart
    f = sqrt(3.0*Npart*Temp/v2)
    v2 = 0.0
    for i in range(0,Npart):
        VX[i] = (VX[i] - vx0)*f
        VY[i] = (VY[i] - vy0)*f
        VZ[i] = (VZ[i] - vz0)*f
        v2 = v2 + (VX[i])**2 + (VY[i])**2 + (VZ[i])**2

    Tempact = v2/(3.0*Npart)

    return VX,VY,VZ,Tempact

def solve(Data,X,Y,Z,VX,VY,VZ,Fx,Fy,Fz,Delt):
    # solve equation of motion for kinetic energy
    Npart = int(Data['NPART'])
    Rho = float(Data['RHO'])
    Box = (Npart/Rho)**(1.0/3.0)

    v2 = 0.0

    for i in range(0,Npart):
        vxt = VX[i]
        vyt = VY[i]
        vzt = VZ[i]
        VX[i] = VX[i] + Delt*Fx[i]
        VY[i] = VY[i] + Delt*Fy[i]
        VZ[i] = VZ[i] + Delt*Fz[i]
        v2 = v2 + ((VX[i]+vxt)**2)/4.0 + \
                         ((VY[i]+vyt)**2)/4.0 + ((VZ[i]+vzt)**2)/4.0

        # put particle in simulation box in case displacement
        # moved it out of box
        # (box origin is 0,0,0 and particle coords are positive)
        X[i] = X[i] + Delt*VX[i]
        if X[i] < 0.0:
            X[i] = X[i] + Box
        if X[i] > Box:
            X[i] = X[i] - Box
        Y[i] = Y[i] + Delt*VY[i]
        if Y[i] < 0.0:
            Y[i] = Y[i] + Box
        if Y[i] > Box:
            Y[i] = Y[i] - Box
        Z[i] = Z[i] + Delt*VZ[i]
        if Z[i] < 0.0:
            Z[i] = Z[i] + Box
        if Z[i] > Box:
            Z[i] = Z[i] - Box

    # this Enkin is in cal/mol
    Enkin = v2/2.0

    return X,Y,Z,VX,VY,VZ,Enkin

def velocs(Data,VX,VY,VZ):
    # scale velocities to maintain temperature
    from math import sqrt
    Npart = int(Data['NPART'])
    SI_temp = float(Data['TEMP'])
    kb = 1.9872 # cal/mol/K
    eps = 1000.0 * float(Data['EPS']) #cal/mol
    Temp = SI_temp/(eps/kb)
    v2 = 0.0
    for i in range(0,Npart):
        vx10 = 10.0 * VX[i]
        vy10 = 10.0 * VY[i]
        vz10 = 10.0 * VZ[i]
        v2 = v2 + (VX[i]**2 + VY[i]**2 + VZ[i]**2)
    Tempact = v2/(3.0*Npart)
    fac = sqrt(Temp/Tempact)
    for i in range(0,Npart):
        VX[i] = VX[i]*fac
        VY[i] = VY[i]*fac
        VZ[i] = VZ[i]*fac

    return VX,VY,VZ

def storemd(Data,X,Y,Z,VX,VY,VZ,sfile,Step):

    # writes configuration to disk
    # output is in reduced units

    Npart = int(Data['NPART'])
    Rho = float(Data['RHO'])
    Box = (Npart/Rho)**(1.0/3.0)

    sfile.write('REMARK  MD Restart Configuration\n')
    sfile.write('%13.9f\n' % Box)
    sfile.write('%8d\n' % Npart)
    sfile.write('%10d\n' % Step)
    for i in range(0,Npart):
        sfile.write('%13.9f %13.9f %13.9f %13.9f %13.9f %13.9f\n' \
                    % (X[i],Y[i],Z[i],VX[i],VY[i],VZ[i]))

