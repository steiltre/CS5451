CC = gcc
CCFLAGS = -lm -Wall -march=native -O2

all: rs_openmp rs_mpi

debug: CCFLAGS += -g
debug: rs_openmp rs_mpi

rs_openmp: rs_openmp.c
	$(CC) -o rs_openmp rs_openmp.c -fopenmp $(CCFLAGS)

rs_mpi: rs_mpi.c
	mpicc -o rs_mpi rs_mpi.c $(CCFLAGS)
