# Imperial College London
# HPC Assignment Task 3
# Dominic Pickford 
# 01272723
# 19/03/2017

import matplotlib.pyplot as plt

# The text file includes a single column of data including first of all every time step followed by the displacement of the mid point at each time step
# The first step is to import all the data and subsequently devide this accross two seperate lists then plot one against the other




########################################## Force Increased from [0 T] ####################################
ALL_DATA = []
with open("Task2_FT.txt") as file:
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

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
ax1.plot(Time, Displacement, 'b')

ax1.set_ylabel('Displacement [m]')
ax1.set_xlabel('Time [s]')
ax1.yaxis.grid(color='gray', linestyle='dashed')
ax1.xaxis.grid(color='gray', linestyle='dashed')
fig1.patch.set_facecolor('white')

########################################## Force Increased from [0 t) ####################################
ALL_DATA2 = []
with open("Task2_HT.txt") as file:
    for line in file:
        line = line.strip() 
        ALL_DATA2.append(line)

length2 = len(ALL_DATA2)

Time2 = []
for i in range(0, length2/2):
	Time2.append(ALL_DATA2[i])

Displacement2 = []
for i in range(length2/2, length2):
	Displacement2.append(ALL_DATA2[i])

fig2 = plt.figure(2)
ax2 = fig1.add_subplot(111)
ax2.plot(Time2, Displacement2)

ax2.set_ylabel('Displacement [m]')
ax2.set_xlabel('Time [s]')
ax2.yaxis.grid(color='gray', linestyle='dashed')
ax2.xaxis.grid(color='gray', linestyle='dashed')
fig2.patch.set_facecolor('white')
############################################################################################################

plt.show()
