all: km_pthreads km_openmp

km_pthreads: km_pthreads.c timing.h
	gcc -o km_pthreads km_pthreads.c timing.h -Wall -O2 -march=native -pthread -lm

km_openmp: km_openmp.c timing.h
	gcc -o km_openmp km_openmp.c timing.h -O2 -Wall -march=native -fopenmp -lm
