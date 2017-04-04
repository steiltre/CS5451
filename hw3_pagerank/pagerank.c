
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "pr_graph.h"
#include "pr_utils.h"
#include "pr_accum.h"


/**
* @brief Compute the PageRank (PR) of a graph.
*
* @param graph The graph.
* @param damping Damping factor (or, 1-restart). 0.85 is typical.
* @param max_iterations The maximium number of iterations to perform.
*
* @return A vector of PR values.
*/
double * pagerank(
    pr_graph const * const graph,
    double const damping,
    int const max_iterations);

/**
 * @brief Creates arrays for All-to-All personalized communication
 *
 * @param accum Process accumulator
 * @param graph The graph.
 * @param sendcounts Number of elements to send to each process
 * @param sdispls Displacement of values to be sent to each process
 * @param recvcounts Number of elements to receive from each process
 * @param rdispls Displacement of values to receive from each process
 */
void CreateCommArrays(
    pr_accum const * const accum,
    pr_graph const * const graph,
    int ** sendcounts,
    int ** sdispls,
    int ** recvcounts,
    int ** rdispls);

/**
 * @brief Write pagerank vector to file
 *
 * @param filename Name of output file
 * @param graph The graph
 * @param PR Pagerank vector
 */
void WriteOutput(
    char const *filename,
    pr_graph * graph,
    double * PR);


int main(
    int argc,
    char * * argv)
{
  MPI_Init(&argc, &argv);

  if(argc == 1) {
    fprintf(stderr, "usage: %s <graph> [output file]\n", *argv);
    return EXIT_FAILURE;
  }

  char * ifname = argv[1];
  char * ofname = NULL;
  if(argc > 2) {
    ofname = argv[2];
  }

  pr_graph * graph = pr_graph_load(ifname);
  if(!graph) {
    return EXIT_FAILURE;
  }

  double * PR = pagerank(graph, 0.85, 100);

  /* write pagerank values */
  if(ofname) {
    WriteOutput( ofname, graph, PR );
  }

  pr_graph_free(graph);
  free(PR);

  MPI_Finalize();

  return EXIT_SUCCESS;
}



double * pagerank(
    pr_graph const * const graph,
    double const damping,
    int const max_iterations)
{
  /* grab graph structures to save typing */
  pr_int const nvtxs = graph->nvtxs;
  pr_int const tvtxs = graph->tvtxs;
  pr_int const * const restrict xadj = graph->xadj;
  pr_int const * const restrict nbrs = graph->nbrs;

  int npes, pid;

  MPI_Comm_size( MPI_COMM_WORLD, &npes );
  MPI_Comm_rank( MPI_COMM_WORLD, &pid );

  int ideal_vtxs = (int) ceil( ((double) tvtxs) / npes );


  /* Create accumulator */
  pr_accum * accum = pr_accum_build(graph, npes);

  int * sendcounts, * sdispls, * recvcounts, * rdispls;
  CreateCommArrays(accum, graph, &sendcounts, &sdispls, &recvcounts, &rdispls);

  pr_int * recv_inds = malloc( rdispls[npes] * sizeof( *recv_inds ) );
  double * recv_vals = malloc( rdispls[npes] * sizeof( *recv_vals ) );

  /* Communication pattern is static, so indices can be sent in advance */
  MPI_Alltoallv( accum->send_ind, sendcounts, sdispls, pr_mpi_int, recv_inds, recvcounts, rdispls, pr_mpi_int, MPI_COMM_WORLD );

  /* Initialize pageranks to be a probability distribution. */
  double * PR = malloc(nvtxs * sizeof(*PR));
  for(pr_int v=0; v < nvtxs; ++v) {
    PR[v] = 1. / (double) tvtxs;
  }

  /* Probability of restart */
  double const restart = (1 - damping) / (double) tvtxs;

  /* Convergence tolerance. */
  double const tol = 1e-12;

  double start = MPI_Wtime();

  for(int i=0; i < max_iterations; ++i) {

    pr_accum_zero_vals(accum);
    for (pr_int v=0; v < nvtxs; ++v) {
      double const num_links = (double)(xadj[v+1] - xadj[v]);
      double const pushing_val = PR[v] / num_links;

      for (pr_int e = xadj[v]; e < xadj[v+1]; ++e) {
        accum->vals[ accum->local_nbrs[e] ] += pushing_val;
      }
    }

    /* Communicate values */
    MPI_Alltoallv( accum->vals, sendcounts, sdispls, MPI_DOUBLE, recv_vals, recvcounts, rdispls, pr_mpi_int, MPI_COMM_WORLD );

    /* Initialize new PR values */
    double norm_changed = 0.;
    double global_norm_changed = 0;

    int * next_ind = malloc( npes * sizeof(*next_ind) );

    for (int j=0; j<npes; j++) {
      next_ind[j] = rdispls[j];
    }

    for (pr_int v=0; v<nvtxs; ++v) {
      double old = PR[v];
      PR[v] = restart;

      for (int j=0; j<npes; j++) {
        if (recv_inds[ next_ind[j] ] == v + pid*ideal_vtxs && next_ind[j] < rdispls[j+1])
        {
          PR[v] += damping * recv_vals[ next_ind[j]++ ];
        }
      }
      norm_changed += (PR[v] - old) * (PR[v] - old);
    }

    free(next_ind);

    /* Communicate global Frobenius norm */
    MPI_Allreduce( &norm_changed, &global_norm_changed, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    global_norm_changed = sqrt(global_norm_changed);

    if( (i > 1 && global_norm_changed < tol) || i == max_iterations - 1 ) {
      if (pid == 0) {
        double end = MPI_Wtime();
        printf("Number of iterations: %i average time: %0.03fs\n", i+1, (end-start)/(i+1));
      }
      break;
    }
  }

  pr_accum_free(accum);
  free(recv_inds);
  free(recv_vals);
  free(sendcounts);
  free(sdispls);
  free(recvcounts);
  free(rdispls);

  return PR;
}


void CreateCommArrays(
    pr_accum const * const accum,
    pr_graph const * const graph,
    int ** sendcounts,
    int ** sdispls,
    int ** recvcounts,
    int ** rdispls)
{
  int npes, pid;

  MPI_Comm_size( MPI_COMM_WORLD, &npes);
  MPI_Comm_rank( MPI_COMM_WORLD, &pid);

  //int ideal_vtxs = (int) ceil( ((double) graph->tvtxs) / npes );

  *sendcounts = malloc( npes * sizeof(int) );
  *sdispls = malloc( npes * sizeof(int) );
  *recvcounts = malloc( npes * sizeof(int) );
  *rdispls = malloc( (npes+1) * sizeof(int) );

  *(sdispls[0]+0) = 0; /* First process's values are stored at beginning of accumulator */

  int p_curr = 0;
  int count = 0;
  pr_int ideal_vtxs = (int) ceil( ((double) graph->tvtxs) / npes );
  for (pr_int v = 0; v < accum->nvals; v++)
  {
    while (accum->send_ind[v] >= (p_curr+1) * ideal_vtxs) {
      *(sendcounts[0]+p_curr) = count;
      p_curr++;
      *(sdispls[0]+p_curr) = *(sdispls[0]+p_curr-1) + count;
      count = 0;
    }
    count++;
  }
  for (int i = p_curr; i<npes; i++) {  // Fill in rest of sendcounts
    *(sendcounts[0]+i) = count;
    *(sdispls[0]+i+1) = *(sdispls[0]+i) + count;
    count = 0;
  }

  MPI_Alltoall( *sendcounts, 1, MPI_INT, *recvcounts, 1, MPI_INT, MPI_COMM_WORLD );

  *rdispls[0] = 0;
  for (int i=1; i<npes+1; i++) {
    *(rdispls[0]+i) = *(rdispls[0]+i-1) + *(recvcounts[0]+i-1);
  }

}

void WriteOutput(
    char const * filename,
    pr_graph * graph,
    double * PR)
{
  int pid, npes;
  MPI_Status status;

  MPI_Comm_size( MPI_COMM_WORLD, &npes );
  MPI_Comm_rank( MPI_COMM_WORLD, &pid );

  if (pid == 0) {
   FILE * fout = fopen(filename, "w");
    if(!fout) {
      fprintf(stderr, "ERROR: could not open '%s' for writing.\n", filename);
    }

    int vtxs;
    int ideal_vtxs = (int) ceil( ((double) graph->tvtxs) / npes );

    double * PR_local = malloc( ideal_vtxs * sizeof(*PR_local) );

    for(pr_int v=0; v < graph->nvtxs; ++v) {
      fprintf(fout, "%0.3e\n", PR[v]);
    }

    for (int i=1; i<npes; ++i) {
      vtxs = GetChunkSize(i, ideal_vtxs, graph->tvtxs);
      MPI_Recv( PR_local, vtxs, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

      for (int v=0; v<vtxs; ++v) {
        fprintf(fout, "%0.3e\n", PR_local[v]);
      }
    }

    fclose(fout);

    free(PR_local);
  }
  else {
    MPI_Send( PR, graph->nvtxs, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
  }
}
