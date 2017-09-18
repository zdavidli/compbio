import numpy as np
from matplotlib import pyplot as plt
with open("./GitHub/compbio/assignments/hw1/LJoutput.dat") as f:
    data = f.read()

data = data.split('\n')
data = [row.replace(" ", "")[1:-1].split(',') for row in data][:-1]

data = np.array([[float(el) for el in row] for row in data])
col1, col2, col3, col4 = [data[:,i] for i in range(4)]


fig1 = plt.figure(figsize=(20,10))
ax1 = fig1.add_subplot(111)
ax1.set_title("Atom Separation vs LJ Potential")    
ax1.set_xlabel('Atom Separation')
ax1.set_ylabel('Potential Energy')
ax1.plot(col1, col2, label='Attractive Potential')
ax1.plot(col1, col3, label='Repulsive Potential')
ax1.plot(col1, col4, label='Overall Energy')
leg = ax1.legend()
    
plt.show()

fig1.savefig('./GitHub/compbio/assignments/hw1/graph1.png')