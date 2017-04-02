
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

  int npes, pid;

  double start = MPI_Wtime();

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
  double * PR_odd = malloc(nvtxs * sizeof(*PR_odd));
  double * PR_even = malloc(nvtxs * sizeof(*PR_even));
  double * PR_old, * PR_new;
  for(pr_int v=0; v < nvtxs; ++v) {
    PR_odd[v] = 1. / (double) tvtxs;
  }

  PR_old = PR_odd;

  /* Probability of restart */
  double const restart = (1 - damping) / (double) tvtxs;


  /* Convergence tolerance. */
  double const tol = 1e-12;

  //double * PR_accum = malloc(nvtxs * sizeof(*PR));
  for(int i=0; i < max_iterations; ++i) {

    int * accum_ind = malloc( npes * sizeof(*accum_ind) );
    for (int i = 0; i<npes; i++) {
      accum_ind[i] = accum->bdry[i];
    }

    if ( i%2 == 0 ) {
      PR_new = PR_even;
      PR_old = PR_odd;
    }
    else {
      PR_new = PR_odd;
      PR_old = PR_even;
    }

    /* Each vertex pushes PR contribution to all outgoing links */
    for(pr_int v=0; v < nvtxs; ++v) {
      double const num_links = (double)(xadj[v+1] - xadj[v]);
      double const pushing_val = PR_old[v] / num_links;

      for(pr_int e=xadj[v]; e < xadj[v+1]; ++e) {
        //pr_accum_add_val(accum, pushing_val, nbrs[e]);

        accum->vals[ accum_ind[ accum->send_proc_ind[e] ] ] = pushing_val;
        accum_ind[ accum->send_proc_ind[e] ]++;
      }
    }

    /* Communicate values */
    MPI_Alltoallv( accum->vals, sendcounts, sdispls, MPI_DOUBLE, recv_vals, recvcounts, rdispls, pr_mpi_int, MPI_COMM_WORLD );

    /* Initialize new PR values */
    double norm_changed = 0.;
    double global_norm_changed = 0;
    for(pr_int v=0; v < nvtxs; ++v) {
      PR_new[v] = restart; // + (damping * PR_accum[v]);
    }

    /* Add contributions to new pagerank */
    for (pr_int j=0; j < rdispls[npes]; ++j) {
      PR_new[ recv_inds[j] - pid*ideal_vtxs ] += damping * recv_vals[j];
    }

    /* Calculate local contribution to Frobenius norm */
    for (pr_int v=0; v < nvtxs; ++v) {
      norm_changed += (PR_new[v] - PR_old[v]) * (PR_new[v] - PR_old[v]);
    }

    /* Communicate global Frobenius norm */
    MPI_Allreduce( &norm_changed, &global_norm_changed, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    global_norm_changed = sqrt(global_norm_changed);

    free(accum_ind);

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

  if (max_iterations == 0)
    return PR_old;
  else
    return PR_new;
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


  *sendcounts = malloc( npes * sizeof(int) );
  *sdispls = malloc( npes * sizeof(int) );
  *recvcounts = malloc( npes * sizeof(int) );
  *rdispls = malloc( (npes+1) * sizeof(int) );

  *(sdispls[0]+0) = 0; /* First process's values are stored at beginning of accumulator */

  for (int i = 0; i < npes; i++) {
    *(sendcounts[0]+i) = accum->bdry[i+1] - accum->bdry[i];
    *(sdispls[0]+i) = accum->bdry[i];
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
