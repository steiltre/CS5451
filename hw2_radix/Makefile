CC = icc
CCFLAGS = -lm -Wall -march=native
RELEASEFLAGS = -O2

all: rs_openmp rs_mpi

debug: CCFLAGS += -g
debug: rs_openmp rs_mpi

rs_openmp: rs_openmp.c rs_shared.h
	$(CC) -o rs_openmp rs_openmp.c rs_shared.h -fopenmp $(CCFLAGS) $(RELEASEFLAGS)

rs_mpi: rs_mpi.c rs_shared.h
	mpicc -o rs_mpi rs_mpi.c rs_shared.h $(CCFLAGS)
