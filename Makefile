CC=g++
CXXFLAGS= -std=c++11 
LDLIBS= -llapack -lblas -lboost_program_options
TARGET=main
task1=main
task2=main

default: $(TARGET) 
all: $(TARGET)
task1: $(task1)
task2: $(task2)

$(TARGET): $(TARGET).cpp
	$(CC) $(TARGET).cpp $(CXXFLAGS) $(LDLIBS) -o $(TARGET) 
	./$(TARGET) --L 10.0 --Nx 4 --A 0.012 --I 0.0000144 --E 210000000000 --T 1.0 --Nt 10000 --rho 7850 

$(task1): $(task1).cpp
	$(CC) $(task1).cpp -std=c++11 -llapack -lblas -lboost_program_options -o $(task1)
	./$(task1) --L 10.0 --Nx 4 --A 0.012 --I 0.0000144 --E 210000000000

$(task2): $(task2).cpp
	$(CC) $(task2).cpp -std=c++11 -llapack -lblas -lboost_program_options -o $(task2)
	./$(task2) --L 10.0 --Nx 4 --A 0.012 --I 0.0000144 --E 210000000000 --T 1.0 --Nt 2000 --rho 7850

.PHONY: clean

clean: 
	rm -f *.o

#.PHONY: clean run

# example
#run: $(TARGET)
#	./$(TARGET) --L 5.0 --Nx 3

#clean:
#	rm -f *.o
#	$(C) -o t1 main.o -std=c++11 -llapack -lblas
