"""Script to process the structures written out during a LINUS simulation
and calculate the torsion angles using utils/torsion.py.

This version should print to the screen the mean and angular deviation of
phi and psi angles for each residue in a LINUS *.out file.

"""

USAGE = """
python print_phi_psi_mean_SD.py [trajectory].pdb > [filename].angstats

  [trajectory].pdb = name of file to which sampled structures were written
            during the course of the simulation (trajectory file).
"""
            
# Import the LINUS modules 
from pylinus_1_0 import linusSim
from pylinus_1_0.utils import torsion

# Import some useful modules
import sys, os, math

# Assign useful variables for conversion
RADIANS_TO_DEGREES = 180.0/math.pi
DEGREES_TO_RADIANS = math.pi/180.0

def get_conformation(out_file):
    # This module processes LINUS trajectory files
    # containing multiple PDB structures.
    #
    # It assumes that each new structure starts with the
    # keyword, MODEL and ends with keyword, ENDMDL

    # Define a usefule variable for reading the file
    nextline = out_file.readline

    # Read through the trajectory to find each structure
    while 1:
        line = nextline()

        # If at start of a structure break to below for reading rest
        # of conformation
        if line[:5] == 'MODEL':
            break

        # Exit while loop if at end of file
        if not line:
            break

    # Exit module if at end of file and return not true
    if not line:
        return 0

    # Open temporary file for writing this conformation as PDB
    tmpfile = open('tmp.pdb', 'w')

    # Loop through each atom in current structure
    while 1:
        line = nextline()
        if line[:6] == 'ENDMDL':
            # If at end of structure return with true
            tmpfile.write('END\n')
            tmpfile.close()
            return 1
        elif line[:4] == 'ATOM':
            # Else write this line to the tmp.pdb file.
            tmpfile.write(line)

    # Notice that there is no return statement at the end of this
    # module. The possible return states (0 or 1) are in the conditional
    # statements above.
    
def get_num_res(out_file):
    # This module processes LINUS trajectory files
    # containing multiple PDB structures.
    #
    # It's sole purpose is to get the number of residues
    #
    # It reads only the first structure in the trajectory
    #
    # It assumes that each new structure starts with the
    # keyword, MODEL and ends with keyword, ENDMDL

    # Define a usefule variable for reading the file
    nextline = out_file.readline

    # Initialize variable
    num_res = 0

    # Read through the trajectory to find each structure
    while 1:
        line = nextline()

        # If at start of a structure break to below for reading rest
        # of conformation
        if line[:5] == 'MODEL':
            break

        # If at end of file exit with message
        # Should never get here if more than one structure in trajectory
        if not line:
            print "No conformations in out file!"
            sys.exit()

    # Loop through each atom in current structure
    while 1:
        line = nextline()
        if line[:6] == 'ENDMDL':
            # If at end of structure return with number of residues
            return num_res
        elif line[:4] == 'ATOM' and line[13:15] == 'CA':
            num_res += 1

def tor(out_file):
    # Read first structure in trajectory to get number of residues
    file = open(out_file, 'rb')
    rn = get_num_res(file)
    file.close()

    # Initialize 2D arrays to store phi and psi for each residue for
    # each conformation
    #
    # Make large enougth for 10000 structures 
    resrange = range(rn)
    phi = []
    for i in resrange:
        phi.append([0.0]*10000)        
    psi = []
    for i in resrange:
        psi.append([0.0]*10000)        

    # Initialize counter for conformation number
    conf = 0

    # Open trajectory again
    file = open(out_file, 'rb')
    while 1:
        # Get a structure and write it as tmp.pdb
        ok = get_conformation(file)

        # If at end of trajectory break out of while loop.
        if not ok: break

        # Run the LINUS module torsion on the structure in tmp.pdb
        # This will make a file called angles.dat 
        torsion.run('tmp.pdb')

        # Open the file with the angle information for this structure
        torfile = open('angles.dat', 'r')

        nextline = torfile.readline
  
        # Initialize counter for residue number
        res = 0
        while 1:
            line = nextline()
            if not line:
                break

            # Don't read the line with labels
            if (line[10:12] != 'SS'):

                # Assign phi and psi to appropriate 2D array index
                # But set first phi and last psi to zero
                phi[res][conf]=(float(line[14:23]))
                if (phi[res][conf] > 999.0):
                    phi[res][conf] = 0.0
                psi[res][conf]=(float(line[23:32]))
                if (psi[res][conf] > 999.0):
                    psi[res][conf] = 0.0

                # Increment residue number
                res += 1

        # Increment conformation number
        conf = conf + 1
 
        # Close and remove the tmp.pdb file
        torfile.close()
        os.system('rm tmp.pdb')

    # Close the trajectory file
    file.close()

    # Print the labels data columns
    print '  Res       N    Meanphi AngDevphi Meanpsi AngDevpsi'

    # Calculate the circular statistics
    for i in range((rn)):
        # Here sumCOS = sum of cos(angle)
        # Here sumSIN = sum of sin(angle)
        sumCOSphi = 0.0
        sumSINphi = 0.0
        sumCOSpsi = 0.0
        sumSINpsi = 0.0

        for j in range(conf):
            sumCOSphi += (math.cos(DEGREES_TO_RADIANS * phi[i][j]))
            sumSINphi += (math.sin(DEGREES_TO_RADIANS * phi[i][j]))
            sumCOSpsi += (math.cos(DEGREES_TO_RADIANS * psi[i][j]))
            sumSINpsi += (math.sin(DEGREES_TO_RADIANS * psi[i][j]))

        # Average of COSphi and SINphi
        COSphi = sumCOSphi / float(conf)
        SINphi = sumSINphi / float(conf)
   
        # Calculate angular deviation
        Rphi = math.sqrt((COSphi*COSphi) + (SINphi*SINphi))
        AngDevphi = (180.0/math.pi) * math.sqrt(2.0 * (1.0 - Rphi))

        # Average of COSpsi and SINpsi
        COSpsi = sumCOSpsi / float(conf)
        SINpsi = sumSINpsi / float(conf)

        # Calculate angular deviation
        Rpsi = math.sqrt((COSpsi*COSpsi) + (SINpsi*SINpsi))
        AngDevpsi = (180.0/math.pi) * math.sqrt(2.0 * (1.0 - Rpsi))

        # Calculate mean of phi
        if sumCOSphi > 0.0:
            meanphi = RADIANS_TO_DEGREES * math.atan(sumSINphi/sumCOSphi)
        elif sumCOSphi < 0.0:
            meanphi = RADIANS_TO_DEGREES * (math.pi + math.atan(sumSINphi/sumCOSphi))
        else:
            meanphi = 999.9

        # Change to familiar phi range
        if meanphi > 180.0:
            meanphi = meanphi - 360.0

        # Calculate mean of psi
        if sumCOSpsi > 0.0:
            meanpsi = RADIANS_TO_DEGREES * math.atan(sumSINpsi/sumCOSpsi)
        elif sumCOSpsi < 0.0:
            meanpsi = RADIANS_TO_DEGREES * (math.pi + math.atan(sumSINpsi/sumCOSpsi))
        else:
            meanpsi = 999.9

        # Change to familiar psi range
        if meanpsi > 180.0:
            meanpsi = meanpsi - 360.0

        # Print formatted results to screen 
        print ('%5d %7d  %8.3f %8.3f %8.3f %8.3f' % (
               (i+1),conf,
               meanphi,AngDevphi,meanpsi,AngDevpsi))

# Run main function with trajectory file as argument input
tor(sys.argv[1])

 
