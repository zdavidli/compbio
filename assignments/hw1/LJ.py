#!/usr/bin/env python3
"""
Program LJ

Calculates components of Lennard-Jones potential and
displays results to screen.

Column 1 = Distance between atoms
Column 2 = Attractive Potential
Column 3 = Repulsive Potential
Column 4 = Total Potential

Original code by Jonathan Farrell, TA 2004

Modified by Patrick Fleming

The Johns Hopkins University
Department of Biophysics
201 Jenkins Hall
3400 North Charles Street
Baltimore, MD 21218

Email: pat.fleming@jhu.edu

"""
# Define the module
def LJ():

     # Set epsilon (units = kcal/mol)
     eps = 0.238
     # Set sigma (units = Angstroms)
     sigma = 3.4

     for r in range(200,800):
          # For a sigma of 3.4 we want to examine interatomic
          # distances from about 2 to 8.
          # A range can only be in integers (line above),
          # so here we make it the float range wanted
          rij = float(r/100.0)

          # Define the pre-exponential factor
          factor = 4.0*eps
          # Calculate the terms
          attractive = factor*(-1*(sigma/rij)**6)
          repulsive = factor*((sigma/rij)**12)
          tot_LJ = attractive + repulsive
#         tot_LJ = 100.0 * tot_LJ

          # Cut off output values at a maximum of 10 for
          # plotting purposes
          if tot_LJ > 10.0:
                tot_LJ = 10.0
          if attractive < -10.0:
                attractive = -10.0
          if repulsive > 10.0:
                repulsive = 10.0

          # Return the values
          print(rij,attractive,repulsive,tot_LJ)

# Call the module for execution
LJ()
