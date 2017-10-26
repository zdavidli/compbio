import math

def torsion():
    for theta in range(0, 360):
        x = math.radians(theta)
        y = 5.0 + 0.5*math.cos(x) + 0.8*math.cos(3.0*x)
        print(theta, y)

torsion()
