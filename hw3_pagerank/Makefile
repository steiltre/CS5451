CC = mpicc

pagerank: pagerank.c pr_graph.c pr_utils.c pr_radix_sort.c pr_accum.c
	$(CC) -std=c99 -march=native -O3 -Wall -Wno-unused-result *.c -o pagerank -lm

debug: pagerank.c pr_graph.c pr_utils.c pr_radix_sort.c pr_accum.c
	$(CC) -g -std=c99 -march=native -O0 -Wall -Wno-unused-result *.c -o pagerank -lm

clean:
	rm -f pagerank
