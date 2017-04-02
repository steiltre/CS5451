#ifndef PAGERANK_ACCUMULATOR_H
#define PAGERANK_ACCUMULATOR_H

#include <stdint.h>
#include <mpi.h>

#include "pr_graph.h"

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

  /** Process IDs associated to each receiving vertex */
  int * send_proc_ind;

  /** Partitions sent to each process */
  pr_int * bdry;

  /** Values to send to vertices. */
  double * vals;
} pr_accum;


/**
 * @brief Set values to 0
 *
 * @param accum Accumulator to zero values in
 */
void pr_accum_zero_vals(
    pr_accum * accum);

/**
 * @brief Build sparsity structure of accumulator
 *
 * @param graph Graph to get sparsity from
 */
pr_accum * pr_accum_build(
    pr_graph const * const graph,
    const int npes);

/**
 * @brief Sort send_ind and remove redundant vertices
 *
 * @param accum Accumulator to condense
 */
void pr_accum_condense(
    pr_accum * accum);

/**
 * @brief Create local array of vertex neighbors
 *
 * @param accum Accumulator
 * @param graph The graph
 */
void pr_accum_local_nbrs(
    pr_accum * accum,
    pr_graph const * const graph,
    int npes);

/**
 * @brief Free all memory allocated to accumulator
 *
 * @param accum Accumulator to free
 */
void pr_accum_free(
    pr_accum * accum);

/**
 * @brief Comparison function for use in quicksort
 *
 * @param a First value
 * @param b Second value
 *
 * @return
 */
int compfunc(
    const void * a,
    const void * b);

#endif
