# script to read the solution from task 1


#file_object = open("T1_sol.txt", "r")
#file_object.read()
#print "done"


lines = []
with open("T1_sol.txt") as file:
    for line in file:
        line = line.strip() #or someother preprocessing
        comp_sol.append(line)

n_elements = len(comp_sol)/3;

with open("T1_asol.txt") as file:
	for line in file:
		line = line.strip()
		ana_sol.append(line)

n_elements = len(ana_sol)

print n_elements




