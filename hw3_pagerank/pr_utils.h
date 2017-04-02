#ifndef PAGERANK_UTILS_H
#define PAGERANK_UTILS_H

#include <stdint.h>
#include <mpi.h>

typedef uint64_t pr_int;


/**
 * @brief Compute the number of chunks to span a set of items
 *
 * @param chunk_size The size of chunks
 * @param num_items The number of items
 *
 * @return The number of chunks. The last chunk may be smaller than chunk_size
 */
int GetNumChunks(
    int const chunk_size,
    int const num_items );


/**
 * @brief Compute the number of items in current chunk. Last chunk may be smaller than others.
 *
 * @param chunk_id The chunk being processed
 * @param chunk_size The preferred size of chunks
 * @param num_items The total number of items
 *
 * @return Number of items to process in current chunk
 */
int GetChunkSize(
    int const chunk_id,
    int const chunk_size,
    int const num_items );

#endif
