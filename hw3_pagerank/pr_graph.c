

/* ensure we have `getline()` */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "pr_graph.h"
#include "pr_utils.h"



pr_graph * pr_graph_load(
    char const * const ifname)
{

  int pid;
  MPI_Status status;

  pr_graph *graph = malloc(sizeof(*graph));

  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  if (pid == 0)
  {
    pr_int *pvtxs, *pedges; /* arrays for storing number of vertices and edges for each process */
    int npes, ideal_vtxs;

    char * line = malloc(1024 * 1024);
    size_t len = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    pvtxs = malloc(npes * sizeof(*pvtxs));
    pedges = malloc(npes * sizeof(*pedges));

    FILE * fin = fopen(ifname, "r");
    if(!fin) {
      fprintf(stderr, "ERROR: could not open '%s' for reading.\n", ifname);
      return NULL;
    }

    pr_graph * send_graph = malloc(sizeof(*send_graph));

    /* read nvtxs and nedges */
    fscanf(fin, "%lu", &(graph->tvtxs));
    fscanf(fin, "%lu", &(graph->tedges));
    fscanf(fin, "\n"); /* make sure we process the newline, too. */

    /* ideal number of vertices per process */
    ideal_vtxs = (pr_int) ceil( ((double) graph->tvtxs) / npes );

    for (int i=0; i<npes; i++)
    {
      /* number of vertices assigned to process i */
      pvtxs[i] = GetChunkSize(i, ideal_vtxs, graph->tvtxs);
      pedges[i] = 0;

      for (pr_int v=0; v < pvtxs[i]; ++v) {
        ssize_t read = getline(&line, &len, fin);

        /* For each edge in line */
        char *ptr = strtok(line, " ");
        while(ptr != NULL) {
          char *end = NULL;
          pr_int const e_id = strtoull(ptr, &end, 10);
          if (ptr == end) {
            break;
          }
          pedges[i]++;
          ptr = strtok(NULL, " ");
        }
      }
    }

    rewind(fin); /* Start reading from beginning of file */
    ssize_t read = getline(&line, &len, fin); /* Read first line containing V and E */

    MPI_Bcast( &(graph->tvtxs), 1, pr_mpi_int, 0, MPI_COMM_WORLD );
    MPI_Bcast( &(graph->tedges), 1, pr_mpi_int, 0, MPI_COMM_WORLD );

    MPI_Scatter( pvtxs, 1, pr_mpi_int, &(graph->nvtxs), 1, pr_mpi_int, 0, MPI_COMM_WORLD );
    MPI_Scatter( pedges, 1, pr_mpi_int, &(graph->nedges), 1, pr_mpi_int, 0, MPI_COMM_WORLD );

    graph->xadj = malloc((graph->nvtxs + 1) * sizeof(*graph->xadj));
    graph->nbrs = malloc(graph->nedges * sizeof(*graph->nbrs));

    /* How many edges we have read. */
    pr_int edge_ptr = 0;


    /* Read in graph one vertex at a time. */
    for(pr_int v=0; v < graph->nvtxs; ++v) {
      ssize_t read = getline(&line, &len, fin);
      if(read == -1) {
        fprintf(stderr, "ERROR: premature EOF at line %lu\n", v+1);
        pr_graph_free(graph);
        return NULL;
      }

      /* Store the beginning of the adjacency list. */
      graph->xadj[v] = edge_ptr;

      /* Check for sinks -- these make pagerank more difficult. */
      if(read == 1) {
        fprintf(stderr, "WARNING: vertex '%lu' is a sink vertex.\n", v+1);
        continue;
      }

      /* Foreach edge in line. */
      char * ptr = strtok(line, " ");
      while(ptr != NULL) {
        char * end = NULL;
        pr_int const e_id = strtoull(ptr, &end, 10);
        /* end of line */
        if(ptr == end) {
          break;
        }
        assert(e_id > 0 && e_id <= graph->tvtxs);

        graph->nbrs[edge_ptr++] = e_id - 1; /* 1 indexed */
        ptr = strtok(NULL, " ");
      }
    }
    assert(edge_ptr == graph->nedges);
    graph->xadj[graph->nvtxs] = graph->nedges;

    for( int i=1; i<npes; i++ ) {
      send_graph->nvtxs = pvtxs[i];
      send_graph->nedges = pedges[i];
      send_graph->xadj = malloc((send_graph->nvtxs + 1) * sizeof(*send_graph->xadj));
      send_graph->nbrs = malloc(send_graph->nedges * sizeof(*send_graph->nbrs));

      /* How many edges we have read. */
      pr_int edge_ptr = 0;

      /* Read in graph one vertex at a time. */
      for(pr_int v=0; v < send_graph->nvtxs; ++v) {
        ssize_t read = getline(&line, &len, fin);
        if(read == -1) {
          fprintf(stderr, "ERROR: premature EOF at line %lu\n", v+1);
          pr_graph_free(send_graph);
          return NULL;
        }

        /* Store the beginning of the adjacency list. */
        send_graph->xadj[v] = edge_ptr;

        /* Check for sinks -- these make pagerank more difficult. */
        if(read == 1) {
          fprintf(stderr, "WARNING: vertex '%lu' is a sink vertex.\n", v+1);
          continue;
        }

        /* Foreach edge in line. */
        char * ptr = strtok(line, " ");
        while(ptr != NULL) {
          char * end = NULL;
          pr_int const e_id = strtoull(ptr, &end, 10);
          /* end of line */
          if(ptr == end) {
            break;
          }
          assert(e_id > 0 && e_id <= graph->tvtxs);

          send_graph->nbrs[edge_ptr++] = e_id - 1; /* 1 indexed */
          ptr = strtok(NULL, " ");
        }
      }
      assert(edge_ptr == send_graph->nedges);
      send_graph->xadj[send_graph->nvtxs] = send_graph->nedges;

      /* Send graph to appropriate process */
      MPI_Send( send_graph->xadj, send_graph->nvtxs + 1, pr_mpi_int, i, 0, MPI_COMM_WORLD );
      MPI_Send( send_graph->nbrs, send_graph->nedges, pr_mpi_int, i, 0, MPI_COMM_WORLD );

      free(send_graph->xadj);
      free(send_graph->nbrs);
    }

    free(line);
    free(pvtxs);
    free(pedges);
    free(send_graph);
  }

  else {
    MPI_Bcast( &(graph->tvtxs), 1, pr_mpi_int, 0, MPI_COMM_WORLD );
    MPI_Bcast( &(graph->tedges), 1, pr_mpi_int, 0, MPI_COMM_WORLD );

    MPI_Scatter( NULL, 0, pr_mpi_int, &(graph->nvtxs), 1, pr_mpi_int, 0, MPI_COMM_WORLD );
    MPI_Scatter( NULL, 0, pr_mpi_int, &(graph->nedges), 1, pr_mpi_int, 0, MPI_COMM_WORLD );

    graph->xadj = malloc( (graph->nvtxs+1) * sizeof(pr_int) );
    graph->nbrs = malloc( graph->nedges * sizeof(pr_int) );

    MPI_Recv( graph->xadj, graph->nvtxs + 1, pr_mpi_int, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
    MPI_Recv( graph->nbrs, graph->nedges, pr_mpi_int, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
  }

  return graph;
}


void pr_graph_free(
    pr_graph * const graph)
{
  free(graph->xadj);
  free(graph->nbrs);
  free(graph);
}


