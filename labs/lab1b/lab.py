#!/usr/bin/env/ python3

#Part 1
print("Hello, World!", type("Hello, World!"))
print(17, type(17))
print(3.14, type(3.14))

#Part 2
message = 'Hello, World!'
n = 17
pi = 3.14
print(message)
print(n)
print(pi)

#Part 3
x = 17+2
print(x)

y = pi/2
print(y)

fruit = 'banana'
bakedgood = ' nut bread'
print(fruit + bakedgood)

#Part 4
print(n, type(n))
n = float(n)
print(n, type(n))

print(pi, type(pi))
pi = int(pi)
print(pi, type(pi))
pi = str(pi)
print(pi, type(pi))

print(str(pi) + str(pi))


#Part 5

import math
print(math.sqrt(25))
pi = math.pi
print(pi)
angle = ((2.0*pi)/360.0) * 180.0
print(angle)
print(math.cos(angle))


#Part 6
def newline():
    print("------")

newline()

print(bakedgood, pi)


#Part 7
def printTwice(arg):
    two = arg + arg
    print(two)

printTwice('spam')

#print(two)


#Part 8

x = 1
y = 2
if x < y:
    print('x < y')
elif x == y:
    print('x = y')
else:
    print('x > y')

#Part 9

#Group 1
def area(radius):
    return math.pi * radius**2

print("Area of circle= ", area(1.4141))

#Group 2
def dist2D(x1, y1, x2, y2):
    dx = x2 - x1
    dy = y2 - y1
    dist2 = dx**2 + dy**2
    dist = math.sqrt(dist2)
    return dist

print("Distance between points = ", dist2D(1,1,2,2))

def dist3D(x1, y1, z1, x2, y2, z2):
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    dist2 = dx**2 + dy**2 + dz**2
    dist = math.sqrt(dist2)
    return dist


#Part 10

def area2(x1, y1, x2, y2):
    radius = dist2D(x1,y1,x2,y2)
    result = area(radius)
    return result

print(area2(1,1,2,2))


#Part 11

for char in fruit:
    print(char)

def countdown(n):
    while n > 0:
        print(n)
        n = n-1

countdown(10)



#Part 12
list1 = [1,2,3,4]
print(list1)

print(list1[0])
print(list1[3])

print(list1[0:-2])


print(" ")
for i in range(len(list1)):
    print(list1[i])

list3 = ['one', 'two', 'three']
print(list3)
print(list3[1])

if 'two' in list3:
    print(1)
else:
    print(0)

list4 = ['four', 'five']
all = list3 + list4
print(all)

print(all[0:4])

del all[1:4]

print(all)

#Part 13

print(list1)
list1.append(5)
print(list1)

print(list1.index(3))

#Part 14

f = open('test.dat', 'w')
f.write('# Mydata \n')
f.write('1.0    1.0\n2.0    4.0')
f.close()


f = open('test.dat', 'r')
print(f.readline())
print(f.readline())
f.close()
f1 = open('test.dat', 'r')
f2 = open('copy.dat', 'w')
while 1:
    line = f1.readline()
    if line == '':
        break
    if line[0] == '#':
        continue
    f2.write(line)
f1.close()
f2.close()




