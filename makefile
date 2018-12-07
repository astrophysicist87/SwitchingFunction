CC = g++
#DEBUG = 
DEBUG = -g
#OMPFLAGS = 
OMPFLAGS = -fopenmp
#OPTFLAGS =
OPTFLAGS = -O3 #optimize
CFLAGS = -c -Wall $(OMPFLAGS) $(OPTFLAGS) -std=c++11 $(DEBUG)
LFLAGS =  $(OMPFLAGS) -lgsl -lgslcblas -Wall $(DEBUG) 
EXECUTABLE = program  #name of executable produced

all:build run

build:
	$(CC) $(CFLAGS) constants.cpp -o constants.o
	$(CC) $(CFLAGS) crossover_line.cpp -o crossover_line.o
	$(CC) $(CFLAGS) main.cpp -o main.o
	$(CC) $(CFLAGS) particle.cpp -o particle.o
	$(CC) $(CFLAGS) particlelist.cpp -o particlelist.o
	$(CC) $(CFLAGS) thermo_pt.cpp -o thermo_pt.o
	$(CC) $(CFLAGS) thermo_exI.cpp -o thermo_exI.o
	$(CC) $(CFLAGS) thermo_exII.cpp -o thermo_exII.o
	$(CC) $(CFLAGS) thermo_pqcd.cpp -o thermo_pqcd.o
	$(CC) $(CFLAGS) thermo_crossover.cpp -o thermo_crossover.o
	$(CC) $(CFLAGS) util.cpp -o util.o
	$(CC) constants.o thermo_pt.o thermo_exI.o crossover_line.o \
	thermo_exII.o thermo_pqcd.o thermo_crossover.o particlelist.o particle.o \
        util.o main.o $(LFLAGS) -o $(EXECUTABLE)


run:
	./$(EXECUTABLE)

clean:
	rm -rf *.o *.out $(EXECUTABLE)




