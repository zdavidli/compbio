"""
LINUS command file to equilibrate structure.

Requires a PDB format file created by RIBOSOME
as an input file.

Output will include:
 [basename].min - PDB format file containing minimal energy structure
 [basename].out - PDB format file containing multiple structures saved
                every n cycles as determined by "numsave" below.
 [basename_xxx_xxx].SWT - File containing secondary strucuture weights
 [basename_xxx_xxx].ps - Postscript plot of data in SWT file
"""

USAGE = """
python linus.py [basename].pdb
"""

# import two required modules
from pylinus_1_0 import linusSim, ran

def gen_seeds():
    # generate seed values for random number generator
    from random import uniform
    s1 = long(uniform(1.0, ran.M1))
    s2 = long(uniform(1.0, ran.M1))
    s3 = long(uniform(1.0, ran.M1))
    s4 = long(uniform(1.0, ran.M2))
    s5 = long(uniform(1.0, ran.M2))
    s6 = long(uniform(1.0, ran.M2))
    return s1, s2, s3, s4, s5, s6

def run(file,seed=None):
    # the main list of commands to give to the program
 
    # define random generator and which seeds to use
    if seed is None:
        gen = apply(ran.generator, gen_seeds())
    else:
        gen = apply(ran.generator, seed)

    # assign the simulation object to variable sim
    sim = linusSim.Linus(file, gen)

    # define a variable equal to the number of residues
    NUMRES = sim.protein.num_residues

    # Do not change anything above this line

    #Beginning of parameter section (these parameters may be edited)

    # save trajectory every numsave cycles and run for numcycles
    sim.set_simulation_parameters(numsave=100,
                                  numcycle=25000)

  # use hydrogen bonds in scoring potential (1) or not (0)
  # score hydrogen bonds only between residue(i) and residue(j)
  # if j-i > hbond_winmin and j-1 < hbond_winmax
    sim.set_hbond_parameters(use_hbond=1,
                             hbond_winmin=5,
                             hbond_winmax=NUMRES)

  # use hydrophobic interactions in scoring potential (1) or not (0)
    sim.set_contact_parameters(use_contact=0,
                               contact_winmax=NUMRES)

    #End of parameter section
    sim.run()

if __name__ == "__main__":
    import sys, os
    filename = sys.argv[1]
    file = os.path.basename(filename)

    basedir = os.getcwd()
    curdir = basedir + os.sep + 'sim_hbond'
    os.mkdir(curdir)
    os.link(filename, curdir+os.sep+file)
    os.chdir(curdir)
    run(file)

