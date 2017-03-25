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
with open("Task2_oscillations.txt") as file:
    for line in file:
        line = line.strip()
        ALL_DATA.append(line)

length = len(ALL_DATA)-1

Loading_Time = []
for i in range(length/2, length):
	Loading_Time.append(ALL_DATA[i])


Amplitude = []
for i in range(0, length/2):
	Amplitude.append(ALL_DATA[i])

print "Loading Time"
for i in Loading_Time:
	print i

print "Amplitude"
for i in Amplitude:
	print i 

print("ALL DATA length: ", length)
print("Amplitude length: ", len(Amplitude))
print("Loading Time length: ", len(Loading_Time))

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.loglog(Loading_Time, Amplitude)
ax1.set_ylabel('Osciallation Amplitude [m]')
ax1.set_xlabel('Loading Time [s]')
ax1.yaxis.grid(color='gray', linestyle='dashed')
ax1.xaxis.grid(color='gray', linestyle='dashed')
fig1.patch.set_facecolor('white')
plt.show()


# End