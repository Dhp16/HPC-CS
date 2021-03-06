# Imperial College London
# HPC Assignment Task 3
# Dominic Pickford 
# 01272723
# 19/03/2017

# Script to plot the solution required for Task 3 
 
import matplotlib.pyplot as plt

# The text file includes a single column of data including first of all the displacement at every point followed by X coordinate associated 
# The first step is to import all the data and subsequently devide this accross two seperate lists then plot one against the other

ALL_DATA = []											
with open("Task4_data.txt") as file:
    for line in file:
        line = line.strip()
        ALL_DATA.append(line)

length = len(ALL_DATA)
print "ALL DATA"
for i in ALL_DATA:
	print i

Displacement = []
for i in range(0, length/2):
	Displacement.append(ALL_DATA[i])

X = []
for i in range(length/2, length):
	X.append(ALL_DATA[i])

print("ALL DATA length: ", length)
print("X length: ", len(X))
print("Displacement length: ", len(Displacement))

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(X, Displacement)
ax1.set_ylabel('Displacement [m]')
ax1.set_xlabel('X [m]')
ax1.yaxis.grid(color='gray', linestyle='dashed')
ax1.xaxis.grid(color='gray', linestyle='dashed')
fig1.patch.set_facecolor('white')
plt.show()

# End