CC=g++
PARA_CC=mpicxx
CXXFLAGS= -std=c++11 
LDLIBS=-llapack -lblas -lboost_program_options
T5LDLIBS=-lscalapack-openmpi -lblacs-openmpi -lblacsCinit-openmpi -llapack -lblas -lboost_program_options
TARGET=compile
OBJ=main.cpp tools.cpp task1.cpp task2.cpp task3.cpp task4.cpp task5.cpp

default: $(TARGET)
all: $(TARGET)

$(TARGET): $(OBJ)
	$(PARA_CC) $(OBJ) $(CXXFLAGS) $(T5LDLIBS) -o $(TARGET)

task1: $(TARGET)
	mpiexec -np 1 ./$(TARGET) --L 10.0 --Nx 24 --A 0.012 --I 0.0000144 --E 210000000000.0 --T 0.1 --Nt 10000 --rho 7850.0 --equation s --scheme n --parallel false
  
task2: $(TARGET)
	mpiexec -np 1 ./$(TARGET) --L 10.0 --Nx 24 --A 0.012 --I 0.0000144 --E 210000000000.0 --T 1.0 --Nt 500000 --rho 7850.0 --equation d --scheme e --parallel false

task3: $(TARGET)
	mpiexec -np 1 ./$(TARGET) --L 10.0 --Nx 24 --A 0.012 --I 0.0000144 --E 210000000000.0 --T 1.0 --Nt 10000 --rho 7850.0 --equation d --scheme i --parallel false

task4: $(TARGET)
	mpiexec -np 2 ./$(TARGET) --L 10.0 --Nx 24 --A 0.012 --I 0.0000144 --E 210000000000.0 --T 1.0 --Nt 500000 --rho 7850.0 --equation d --scheme e --parallel true

task5: $(TARGET)
	mpiexec -np 4 ./$(TARGET) --L 10.0 --Nx 24 --A 0.012 --I 0.0000144 --E 210000000000.0 --T 1.0 --Nt 10000 --rho 7850.0 --equation d --scheme i --parallel true

# Paramter notes:
# 	equation: 's' = static   'd' = dynamic
# 	scheme: 'e' = explicit  'i' = implicit 'n' = neither
#   parallel: false || true