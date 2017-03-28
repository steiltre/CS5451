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
  accum->nvals++;
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

void pr_accum_condense(
    pr_accum * accum)
{

  radix_sort(accum->send_ind, accum->nvals);

  /* DOESN'T HANDLE EMPTY ACCUMULATOR CORRECTLY */

  /* Counter for number of unique indices in send_ind */
  pr_int count = 1;

  for (pr_int i=1; i<accum->nvals; i++) {
    if (accum->send_ind[i] != accum->send_ind[i-1])
      ++count;
  }

  pr_int * new_send_ind = malloc( count * sizeof(pr_int) );

  count = 1;
  new_send_ind[0] = accum->send_ind[0];

  /* Create new send_ind array with only unique vertex IDs */
  for (pr_int i=1; i<accum->nvals; i++) {
    if (accum->send_ind[i] != accum->send_ind[i-1]) {
      new_send_ind[count] = accum->send_ind[i];
      count++;
    }
  }

  accum->nvals = count;
  free(accum->send_ind);
  accum->send_ind = new_send_ind;
  accum->vals = malloc( accum->nvals * sizeof(pr_int) );

}

void pr_accum_free(
    pr_accum *accum)
{
  free(accum->send_ind);
  free(accum->vals);
  free(accum);
}
