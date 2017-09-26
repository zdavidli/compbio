"""
David Li
dli44@jhu.edu
Computational Biology HW2

Generating random Gaussian data
"""

#Import random module for random.gauss() function
import random

#set up variables
mu = 8
sigma = 3
count=10000

#List comprehension generates a new random number $count times.
values = [(random.gauss(mu, sigma)) for i in range(count)]

#list with same length as values of 1s and 0s
#Index i is 0 if it's less than or equal to mu+sigma, 1 if greater
onesigma = [int(num > (mu+sigma)) for num in values]
greaterone = sum(onesigma)

#list with same length as values of 1s and 0s
#Index i is 0 if it's less than or equal to mu+2*sigma, 1 if greater
twosigma = [int(num > (mu+2*sigma)) for num in values]
greatertwo = sum(twosigma)

#Print percentage greater than mu+sigma and mu+2sigma with formatting
print("{:06.2f}% of the generated numbers are greater than mu+sigma".format(greaterone/count * 100)) #multiplying by 100 to get percentages
print("{:06.2f}% of the generated numbers are greater than mu+2*sigma".format(greatertwo/count * 100))
