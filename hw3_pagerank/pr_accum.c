#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pr_accum.h"
#include "pr_utils.h"
#include "pr_radix_sort.h"

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

  int * accum_ind = malloc( npes * sizeof(*accum_ind) );

  accum->nvals = graph->nedges;
  for (int i=0; i<npes+1; i++) {
    accum->bdry[i] = 0;
  }

  int ideal_vtxs = (int) ceil( ((double) graph->tvtxs ) / npes );
  int proc_ind = 0;
  for (pr_int v = 0; v < graph->nvtxs; v++) {
    for (pr_int e = graph->xadj[v]; e < graph->xadj[v+1]; e++) {
      while (graph->nbrs[e] >= ideal_vtxs * (proc_ind+1) || graph->nbrs[e] < ideal_vtxs * proc_ind)
      {
        proc_ind = (proc_ind+1)%npes;
      }
      accum->bdry[proc_ind+1]++;
      accum->send_proc_ind[e] = proc_ind;
    }
  }

  for (int i = 0; i<npes; i++) {
    accum->bdry[i+1] += accum->bdry[i];
  }

  for (int i = 0; i<npes; i++) {
    accum_ind[i] = accum->bdry[i];
  }

  for (pr_int e = 0; e < graph->nedges; e++) {
    accum->send_ind[ accum_ind[ accum->send_proc_ind[e] ]++ ] = graph->nbrs[e];
  }

  free(accum_ind);

  return accum;
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
