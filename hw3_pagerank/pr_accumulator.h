#ifndef PAGERANK_ACCUMULATOR_H
#define PAGERANK_ACCUMULATOR_H

#include <stdint.h>
#include "mpi.h"

typedef uint64_t pr_int;

/**
 * @brief Structure storing values to send to other processes
 */
typedef struct
{
  /** The number of values stored. */
  pr_int nvals;

  /** Indices of vertices to send values to */
  pr_int * send_ind;

  /** Values to send to vertices. */
  double * vals;
} pr_accum;


/**
 * @brief Add an index to the array of indices
 *
 * @param accum Accumulator to add index to
 * @param vtx The index to be added
 */
void pr_accum_add_vtx(
    pr_accum * accum,
    pr_int const vtx);

/**
 * @brief Add a value to the array of values
 *
 * @param accum Accumulator to add value to
 * @param val Value to add
 * @param vtx Vertex ID corresponding to val
 */
void pr_accum_add_val(
    pr_accum * accum,
    double const val,
    pr_int const vtx);

/**
 * @brief Set values to 0
 *
 * @param accum Accumulator to zero values in
 */
void pr_accum_zero_vals(
    pr_accum * accum);

/**
 * @brief Sort send_ind and remove redundant vertices
 *
 * @param accum Accumulator to condense
 */
void pr_accum_condense(
    pr_accum * accum);

/**
 * @brief Free all memory allocated to accumulator
 *
 * @param accum Accumulator to free
 */
void pr_accum_free(
    pr_accum * accum);

#endif
