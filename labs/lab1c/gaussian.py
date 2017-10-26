import random

def Gaussian():
    mu = 3
    sigma = 1
    for count in range(1000):
        print(random.gauss(mu, sigma))

Gaussian()
