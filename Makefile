CC=g++
PARA_CC=mpicxx
CXXFLAGS= -std=c++11 
LDLIBS=-llapack -lblas -lboost_program_options
T5LDLIBS=-lscalapack-openmpi -lblacs-openmpi -lblacsCinit-openmpi -llapack -lblas -lboost_program_options
TARGET=compile
OBJ=main.cpp

default: $(TARGET)
all: $(TARGET)

$(TARGET): $(OBJ)
	$(PARA_CC) $(OBJ) $(CXXFLAGS) $(T5LDLIBS) -o $(TARGET)
	
task1: $(TARGET)
	mpiexec -np 1 ./$(TARGET) --L 10.0 --Nx 24 --A 0.012 --I 0.0000144 --E 210000000000.0 --T 0.1 --Nt 1 --rho 0.0 --equ static --scheme neither
  
task2: $(TARGET)
	mpiexec -np 1 ./$(TARGET) --L 10.0 --Nx 24 --A 0.012 --I 0.0000144 --E 210000000000.0 --T 1.0 --Nt 10000 --rho 7850.0 --equ static --scheme explicit

task3: $(TARGET)
	mpiexec -np 1 ./$(TARGET) --L 10.0 --Nx 24 --A 0.012 --I 0.0000144 --E 210000000000.0 --T 1.0 --Nt 10000 --rho 7850.0 --equ dynamic --scheme implicit

task4: $(TARGET)
	mpiexec -np 1 ./$(TARGET) --L 10.0 --Nx 24 --A 0.012 --I 0.0000144 --E 210000000000.0 --T 1.0 --Nt 10000 --rho 7850.0 --equ dynamic --scheme explicit --parallel true

task5: $(TARGET)
	mpiexec -np 1 ./$(TARGET) --L 10.0 --Nx 24 --A 0.012 --I 0.0000144 --E 210000000000.0 --T 1.0 --Nt 10000 --rho 7850.0 --equ dynamic --scheme implicit --parallel true