
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

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
    pr_int ** sendcounts,
    pr_int ** sdispls,
    pr_int ** recvcounts,
    pr_int ** rdispls);


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
    FILE * fout = fopen(ofname, "w");
    if(!fout) {
      fprintf(stderr, "ERROR: could not open '%s' for writing.\n", ofname);
      return EXIT_FAILURE;
    }
    for(pr_int v=0; v < graph->nvtxs; ++v) {
      fprintf(fout, "%0.3e\n", PR[v]);
    }
    fclose(fout);
  }

  free(PR);

  return EXIT_SUCCESS;
}



double * pagerank(
    pr_graph const * const graph,
    double const damping,
    int const max_iterations)
{
  /* grab graph structures to save typing */
  pr_int const nvtxs = graph->nvtxs;
  pr_int const * const restrict xadj = graph->xadj;
  pr_int const * const restrict nbrs = graph->nbrs;

  /* Create accumulator */
  //pr_accum * accum = malloc(sizeof(pr_accum));
  pr_accum * accum = pr_accum_build(graph);
  pr_accum_zero_vals(accum);
  //accum->nvals = 0;
  //pr_accum_free(accum);

  pr_int * sendcounts, * sdispls, * recvcounts, * rdispls;
  CreateCommArrays(accum, graph, &sendcounts, &sdispls, &recvcounts, &rdispls);

  /* Initialize pageranks to be a probability distribution. */
  double * PR = malloc(nvtxs * sizeof(*PR));
  for(pr_int v=0; v < nvtxs; ++v) {
    PR[v] = 1. / (double) nvtxs;
  }

  /* Probability of restart */
  double const restart = (1 - damping) / (double) nvtxs;


  /* Convergence tolerance. */
  double const tol = 1e-12;

  //double * PR_accum = malloc(nvtxs * sizeof(*PR));
  //for(int i=0; i < max_iterations; ++i) {

  //  for(pr_int v=0; v < nvtxs; ++v) {
  //    PR_accum[v] = 0.;
  //  }

  //  /* Each vertex pushes PR contribution to all outgoing links */
  //  for(pr_int v=0; v < nvtxs; ++v) {
  //    double const num_links = (double)(xadj[v+1] - xadj[v]);
  //    double const pushing_val = PR[v] / num_links;

  //    for(pr_int e=xadj[v]; e < xadj[v+1]; ++e) {
  //      PR_accum[nbrs[e]] += pushing_val;
  //    }
  //  }

  //  /* Finalize new PR values */
  //  double norm_changed = 0.;
  //  for(pr_int v=0; v < nvtxs; ++v) {
  //    double const old = PR[v];
  //    PR[v] = restart + (damping * PR_accum[v]);

  //    norm_changed = (PR[v] - old) * (PR[v] - old);
  //  }
  //  norm_changed = sqrt(norm_changed);

  //  if(i > 1 && norm_changed < tol) {
  //    break;
  //  }
  //}

  //free(PR_accum);
  pr_accum_free(accum);
  return PR;
}


void CreateCommArrays(
    pr_accum const * const accum,
    pr_graph const * const graph,
    pr_int ** sendcounts,
    pr_int ** sdispls,
    pr_int ** recvcounts,
    pr_int ** rdispls)
{
  int npes;

  MPI_Comm_size( MPI_COMM_WORLD, &npes);

  pr_int ideal_vtxs = graph->tvtxs / npes;

  *sendcounts = malloc( npes * sizeof(pr_int) );
  *sdispls = malloc( npes * sizeof(pr_int) );
  *recvcounts = malloc( npes * sizeof(pr_int) );
  *rdispls = malloc( npes * sizeof(pr_int) );

  *(sdispls[0]+0) = 0; /* First process's values are stored at beginning of accumulator */

  pr_int p_curr = 0;
  pr_int count = 0;
  for (pr_int v = 0; v < accum->nvals; v++)
  {
    while (accum->send_ind[v] > (p_curr+1) * ideal_vtxs) {
      *(sendcounts[0]+p_curr) = count;
      p_curr++;
      *(sdispls[0]+p_curr) = *(sdispls[0]+p_curr-1) + count;
      count = 0;
    }
    count++;
  }
  *(sendcounts[0]+p_curr) = count;

  MPI_Alltoall( *sendcounts, 1, pr_mpi_int, *recvcounts, 1, pr_mpi_int, MPI_COMM_WORLD );

  *rdispls[0] = 0;
  for (pr_int i=1; i<npes; i++) {
    *(rdispls[0]+i) = *(rdispls[0]+i-1) + *(recvcounts[0]+i-1);
  }

}
