#include <stdlib.h>
#include <stdio.h>

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
    pr_graph const * const graph)
{
  pr_accum * big_accum = malloc( sizeof(pr_accum) );
  big_accum->nvals = 0;
  big_accum->send_ind = malloc( graph->nedges * sizeof(pr_int) );
  big_accum->vals = malloc( graph->nedges * sizeof(pr_int) );

  /* Add incident vertex for each edge to accumulator */
  for (pr_int e = 0; e < graph->nedges; e++) {
    pr_accum_add_vtx(big_accum, graph->nbrs[e]);
  }

  pr_accum * accum = pr_accum_condense(big_accum);

  return accum;
}

pr_accum * pr_accum_condense(
    pr_accum * accum)
{

  pr_accum * new_accum = malloc( sizeof(pr_accum) );
  new_accum->nvals = 0;

  radix_sort(accum->send_ind, accum->nvals);

  /* DOESN'T HANDLE EMPTY ACCUMULATOR CORRECTLY */

  /* Counter for number of unique indices in send_ind */
  pr_int count = 1;

  for (pr_int i=1; i<accum->nvals; i++) {
    if (accum->send_ind[i] != accum->send_ind[i-1])
      ++count;
  }

  new_accum->send_ind = malloc( count * sizeof(pr_int) );

  //new_accum->send_ind[0] = accum->send_ind[0];
  pr_accum_add_vtx( new_accum, accum->send_ind[0]);

  /* Create new send_ind array with only unique vertex IDs */
  for (pr_int i=1; i<accum->nvals; i++) {
    if (accum->send_ind[i] != accum->send_ind[i-1]) {
      //new_accum->send_ind[count] = accum->send_ind[i];
      //count++;
      pr_accum_add_vtx( new_accum, accum->send_ind[i] );
    }
  }

  //new_accum->nvals = count;
  pr_accum_free(accum);
  new_accum->vals = malloc( accum->nvals * sizeof(pr_int) );

  return new_accum;

}

void pr_accum_free(
    pr_accum *accum)
{
  free(accum->send_ind);
  free(accum->vals);
  free(accum);
}
