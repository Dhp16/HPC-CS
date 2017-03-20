# Imperial College London
# HPC Assignment Task 3
# Dominic Pickford 
# 01272723
# 19/03/2017

# Script to plot the solution required for Task 3 

import matplotlib.pyplot as plt

# The text file includes a single column of data including first of all every time step followed by the displacement of the mid point at each time step
# The first step is to import all the data and subsequently devide this accross two seperate lists then plot one against the other

ALL_DATA = []											
with open("Task3_HT.txt") as file:
    for line in file:
        line = line.strip()
        ALL_DATA.append(line)

length = len(ALL_DATA)

Time = []
for i in range(0, length/2):
	Time.append(ALL_DATA[i])

Displacement = []
for i in range(length/2, length):
	Displacement.append(ALL_DATA[i])

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(Time, Displacement)
ax1.set_ylabel('Displacement [m]')
ax1.set_xlabel('time [s]')
ax1.yaxis.grid(color='gray', linestyle='dashed')
ax1.xaxis.grid(color='gray', linestyle='dashed')
fig1.patch.set_facecolor('white')

plt.show()


# End