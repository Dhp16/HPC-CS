CC=g++
PARA_CC=mpicxx
CXXFLAGS= -std=c++11 
LDLIBS= -llapack -lblas -lboost_program_options
TARGET=main
task4=main
task1=main
task2=main
task3=main


default: $(TARGET) 
all: $(TARGET)
task4: $(task4)
task1: $(task1)
task2: $(task2)
task3: $(task3)

$(TARGET): $(TARGET).cpp
	$(PARA_CC) $(TARGET).cpp $(CXXFLAGS) $(LDLIBS) -o $(TARGET) 
	./$(TARGET) --L 10.0 --Nx 6 --A 0.012 --I 0.0000144 --E 210000000000 --T 1.0 --Nt 10000 --rho 7850 

$(task1): $(task1).cpp
	$(PARA_CC) $(task1).cpp -$(CXXFLAGS) $(LDLIBS) -o $(task1)
	./$(task1) --L 12.0 --Nx 6 --A 0.012 --I 0.0000144 --E 210000000000 --T 1.0 --Nt 10000 --rho 7850
 
$(task2): $(task2).cpp 
	$(CC) $(task2).cpp -$(CXXFLAGS) $(LDLIBS) -o $(task2)
	./$(task2) --L 10.0 --Nx 6 --A 0.012 --I 0.0000144 --E 210000000000 --T 1.0 --Nt 10000 --rho 7850

$(task3): $(task3).cpp
	$(PARA_CC) $(task3).cpp -$(CXXFLAGS) $(LDLIBS) -o $(task3)
	./$(task3) --L 10.0 --Nx 6 --A 0.012 --I 0.0000144 --E 210000000000 --T 1.0 --Nt 10000 --rho 7850

$(task4): $(task4).cpp
	$(PARA_CC) $(task4).cpp -$(CXXFLAGS) $(LDLIBS) -o $(task4)
	mpiexec -np 1 ./$(task4) --L 10.0 --Nx 6 --A 0.012 --I 0.0000144 --E 210000000000 --T 1.0 --Nt 10000 --rho 7850

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
