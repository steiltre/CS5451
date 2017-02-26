CC = gcc
CCFLAGS = -O2 -lm -Wall -march=native

all: rs_openmp

debug: CCFLAGS += -g
debug: rs_openmp

rs_openmp: rs_openmp.c rs_shared.h
	$(CC) -o rs_openmp rs_openmp.c rs_shared.h -fopenmp $(CCFLAGS)
