#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pr_accum.h"
#include "pr_utils.h"
#include "pr_radix_sort.h"

#define SORT

void pr_accum_add_vtx(
    pr_accum * accum,
    pr_int const vtx)
{
  accum->send_ind[accum->nvals] = vtx;
  ++accum->nvals;
}

void pr_accum_add_val(
    pr_accum * accum,
    double const val,
    pr_int const vtx)
{
  pr_int ind = binary_search(accum->send_ind, accum->nvals, vtx);

  if (ind == -1) {
    fprintf(stderr, "ERROR: could not locate '%lu' in send_ind array.\n", vtx);
  }
  else {
    accum->vals[ind] += val;
  }
}

void pr_accum_zero_vals(
    pr_accum * accum)
{
  for (pr_int i=0; i<accum->nvals; i++)
  {
    accum->vals[i] = 0;
  }
}

pr_accum *  pr_accum_build(
    pr_graph const * const graph,
    const int npes)
{
  pr_accum * accum = malloc( sizeof(pr_accum) );
  accum->send_ind = malloc( graph->nedges * sizeof(*accum->send_ind) );
  accum->vals = malloc( graph->nedges * sizeof(*accum->vals) );
  accum->bdry = malloc( (npes+1) * sizeof(*accum->bdry) );
  accum->send_proc_ind = malloc( graph->nedges * sizeof(*accum->send_proc_ind) );
  accum->local_nbrs = malloc( graph->nedges * sizeof(*accum->local_nbrs) );

  int * accum_ind = malloc( npes * sizeof(*accum_ind) );

#ifdef SORT
  accum->nvals = 0;
  /* Add incident vertex for each edge to accumulator */
  for (pr_int e = 0; e < graph->nedges; e++) {
    pr_accum_add_vtx(accum, graph->nbrs[e]);
  }

  pr_accum_condense(accum);
  pr_accum_local_nbrs(accum, graph);

#else
  accum->nvals = graph->nedges;
  for (int i=0; i<npes+1; i++) {
    accum->bdry[i] = 0;
  }

  int ideal_vtxs = (int) ceil( ((double) graph->tvtxs ) / npes );
  int proc_ind = 0;
  for (pr_int e = 0; e < graph->nedges; e++) {
    while (graph->nbrs[e] >= ideal_vtxs * (proc_ind+1) || graph->nbrs[e] < ideal_vtxs * proc_ind)
    {
      proc_ind = (proc_ind+1)%npes;
    }
    accum->bdry[proc_ind+1]++;
#ifdef PRECOMP_SEND_PROC
    accum->send_proc_ind[e] = proc_ind;
#endif
  }

  for (int i = 0; i<npes; i++) {
    accum->bdry[i+1] += accum->bdry[i];
  }

  for (int i = 0; i<npes; i++) {
    accum_ind[i] = accum->bdry[i];
  }

  for (pr_int e = 0; e < graph->nedges; e++) {
    int ind = graph->nbrs[e] / ideal_vtxs;
    accum->send_ind[ accum_ind[ind] ] = graph->nbrs[e];
    accum_ind[ind]++;
  }
#endif

  free(accum_ind);

  return accum;
}

void pr_accum_condense(
    pr_accum * accum)
{

  if ( accum->nvals == 0 )
    return;  /* Accumulator is already condensed if it is empty */

  pr_int * new_send_ind;

  //radix_sort(accum->send_ind, accum->nvals);
  qsort(accum->send_ind, accum->nvals, sizeof(*accum->send_ind), compfunc);

  /* Counter for number of unique indices in send_ind */
  int count = 1;

  for (pr_int i=1; i<accum->nvals; i++) {
    if (accum->send_ind[i] != accum->send_ind[i-1])
      ++count;
  }

  new_send_ind = malloc( count * sizeof(pr_int) );

  new_send_ind[0] = accum->send_ind[0];

  count = 1;

  /* Create new send_ind array with only unique vertex IDs */
  for (pr_int i=1; i<accum->nvals; i++) {
    if (accum->send_ind[i] != accum->send_ind[i-1]) {
      new_send_ind[count] = accum->send_ind[i];
      count++;
    }
  }

  free(accum->send_ind);
  free(accum->vals);
  accum->send_ind = new_send_ind;
  accum->nvals = count;
  accum->vals = malloc( accum->nvals * sizeof(pr_int) );

}

void pr_accum_local_nbrs(
    pr_accum * accum,
    pr_graph const * const graph)
{
  /* Create array of pointers from a vertex to accum->vals.
   * This array is essentially "graph->nbrs" except pointing to accumulator indices rather than other vertices.
   * This structure can be built before iterations and used throughout the calculation because
   * the communication pattern remains the same across iterations.
   */

  int count = 0;
  for (pr_int v=0; v<graph->nvtxs; ++v) {
    for (pr_int e=graph->xadj[v]; e < graph->xadj[v+1]; ++e) {
      pr_int ind = binary_search(accum->send_ind, accum->nvals, graph->nbrs[e]);

      if (ind == -1) {
        fprintf(stderr, "ERROR: could not locate '%lu' in send_ind array.\n", graph->nbrs[e]);
      }
      else {
        accum->local_nbrs[count++] = ind;
      }
    }
  }
}

void pr_accum_free(
    pr_accum *accum)
{
  free(accum->send_ind);
  free(accum->vals);
  free(accum->bdry);
  free(accum->send_proc_ind);
  free(accum);
}

int compfunc(
    const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}
