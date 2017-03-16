# script to read and plot the solution from task 1

import matplotlib.pyplot as plt



#file_object = open("T1_sol.txt", "r")
#file_object.read()
#print "done"

comp_sol = []
with open("T1_csol.txt") as file:
    for line in file:
        line = line.strip() #or someother preprocessing
        comp_sol.append(line)
Xc = []
with open("X_csol.txt") as file:
	for line in file:
		x = line.strip()
		Xc.append(x)


ana_sol = []
with open("T1_asol.txt") as file:
	for line in file:
		line = line.strip()
		ana_sol.append(line)

Xa = []
with open("X.txt") as file:
	for line in file:
		x = line.strip()
		Xa.append(x)

print "Xc: ", len(Xc)
print "ana_sol ", len(ana_sol)
print "comp_sol ", len(comp_sol)

n_elements = len(ana_sol)

#print "Computational Solution", "   Analytical Solution"
#for i in range(0, n_elements):
#	print comp_sol[i], "  ", ana_sol[i]

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(Xc, comp_sol, label="Computational Solution")
ax1.hold(True)
ax1.plot(Xa, ana_sol, label = "Analytical Solution")

ax1.set_ylabel('Displacement [m]')
ax1.set_xlabel('x [m]')
ax1.yaxis.grid(color='gray', linestyle='dashed')
ax1.xaxis.grid(color='gray', linestyle='dashed')
ax1.legend(loc=1)
fig1.patch.set_facecolor('white')


ax1.hold(False)
plt.show()

# comma delimited