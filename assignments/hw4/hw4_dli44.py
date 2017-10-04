# parameters
# Npart : total number of particles
# Temp  : temperature
# Rho   : volume density
# Eps   : epsilon, max energy of particle attraction (kcal/mol)
# Sig   : distanece between particles at zeo attraction
Data={'NPART': 144, 'TEMP': 298, 'RHO': 0.6, 'EPS': 0.238, 'SIG': 3.4}
Npart = int(Data['NPART'])
Temp = int(Data['TEMP'])
Rho = float(Data['RHO'])
sigma = float(Data['SIG']) #Angstrom
eps = 1000.0 * float(Data['EPS']) #cal/mol
# any local variables that were used often and unchanged between subroutines were instead declared as globals here
# elements of periodic boundary conditions
Box = (Npart/Rho)**(1.0/3.0)
Hbox = Box/2.0
Rc = Hbox
Rc2 = Rc*Rc
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

def main():
    # initialize empty lists for coordinates
    X = []
    Y = []
    Z = []
    # open and read pdb file, write to coordinate lists
    ff = open('Hmwk4.pdb', 'r')
    nextline = ff.readline  
    while 1:
        line = nextline()
        if not line:
            break
        line = line.strip()
        line = line.split()
        if line[0].upper() == 'ATOM':
            # convert to reduced distance units
            x = float(line[5])/sigma
            y = float(line[6])/sigma
            z = float(line[7])/sigma
            
            X.append(x)
            Y.append(y)
            Z.append(z)
    ff.close()

    # calculate total energy
    En,Vir = totergmc(Data,X,Y,Z)
    # En is in reduced units
    SI_en = (eps * En) / 1000.0 #kcal/mol
    print(('Total potential energy of the configuration: %8.3f kcal/mol' % SI_en))
main()