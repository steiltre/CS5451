CC = icc
CCFLAGS = -O2 -lm -Wall -march=native

all: rs_openmp rs_mpi

debug: CCFLAGS += -g
debug: rs_openmp rs_mpi

rs_openmp: rs_openmp.c rs_shared.h
	$(CC) -o rs_openmp rs_openmp.c rs_shared.h -fopenmp $(CCFLAGS)

rs_mpi: rs_mpi.c rs_shared.h
	mpicc -o rs_mpi rs_mpi.c rs_shared.h
