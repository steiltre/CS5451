#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "rs_shared.h"



/**
 * @brief Number of bits to sort with at one time
 */
static int const BIT_CHUNK_SIZE = 4;

/**
 * @brief Width of integers being sorted
 */
static int const NUM_BITS = 32;


/**
 * @brief Dense integer matrix structure
 */
typedef struct
{
    int nrows;
    int ncols;
    int *vals;
} rs_matrix;

/**
 * @brief A wrapper around 'posix_memalign()' to get aligned memory.
 *
 * @param byest How many bytes to allocate
 *
 * @return Allocated memory
 */
void * rs_malloc(
        size_t const bytes)
{
    void *ptr;
    int success = posix_memalign(&ptr, 64, bytes);
    if (success != 0)
    {
        fprintf(stderr, "ERROR: posix_memalign() returned %d.\n", success);
        exit(1);
    }
    return ptr;
}


/**
 * @brief Free memory allocated by 'rs_malloc()'.
 *
 * @param ptr The pointer to free.
 */
void rs_free(
        void *ptr)
{
    free(ptr);
}


/**
 * @brief Allocate a dense matrix
 *
 * @param nrows The number of rows in the matrix
 * @param ncols The number of columns in the matrix
 *
 * @return The allocated matrix
 */
rs_matrix * rs_matrix_alloc(
        int const nrows,
        int const ncols)
{
    rs_matrix * matrix = rs_malloc(sizeof(*matrix));

    matrix->nrows = nrows;
    matrix->ncols = ncols;
    matrix->vals = rs_malloc(nrows * ncols * sizeof(*matrix->vals));

    return matrix;
}

/**
 * @brief Free memory allocated by 'rs_matrix_alloc()'
 *
 * @param matrix The matrix to free
 */
void rs_matrix_free(
        rs_matrix * matrix)
{
    if(matrix == NULL)
    {
        return;
    }

    rs_free(matrix->vals);
    rs_free(matrix);
}

/**
 * @brief Compute the number of chunks to span a set of items
 *
 * @param chunk_size The size of chunks
 * @param num_items The number of items
 *
 * @return The number of chunks. The last chunk may be smaller than chunk_size
 */
static inline int GetNumChunks(
        int const chunk_size,
        int const num_items)
{
    int num_chunks = num_items / chunk_size;
    if (num_items % chunk_size > 0)
    {
        ++num_chunks;
    }
    return num_chunks;
}

/**
 * @brief Compute the number of items in current chunk. Last chunk may be smaller than others.
 *
 * @param chunk_id The chunk being processed
 * @param chunk_size The preferred size of chunks
 * @param num_items The total number of items
 *
 * @return Number of items to process in current chunk
 */
static inline int GetChunkSize(
        int const chunk_id,
        int const chunk_size,
        int const num_items)
{
    int stop = chunk_size;
    if(stop + (chunk_id * chunk_size) > num_items) {
        stop = num_items % chunk_size;
    }
    return stop;
}

/**
 * @brief Compute mask for bitwise & operation
 *
 * @param num_bits Number of bits being sorted at once
 *
 * @return Mask for use with num_bits
 */
static inline int GetMask(
        int const num_bits)
{
    int mask = 0;
    for (int j=0; j<num_bits; j++)
    {
        mask += pow(2,j);
    }

    return mask;
}

/**
 * @brief Sort an array of numbers according to num_bins sequential bits starting at ith bit
 *
 * @param arr Array to be sorted
 * @param i Bit to sort by
 * @param num Number of numbers in arr
 * @param local_num_bins Number of bins to place numbers in
 * @param bin_start Array giving starting location in arr for each bin
 * @param mask Mask to use for bitwise & operation
 * @param nthreads Number of threads
 */
void SortOnRadix(
        uint32_t *arr,
        int i,
        int num,
        int local_num_bins,
        rs_matrix *bin_start,
        int mask)
{
    uint32_t *output;
    int j;

    output = (uint32_t *) malloc(num * sizeof(uint32_t));

    for (j=0; j<num; j++)
    {
        output[j] = 0;
    }

    // Scan array putting values in order within local array for output
    for (j=0; j<num; j++)
    {
        int process_id;
        MPI_Comm_rank(MPI_COMM_WORLD, &process_id);

        for (int k=0; k<local_num_bins; k++)
        {
            if ( ((arr[j] >> i * BIT_CHUNK_SIZE) & mask) == k)
            {
                output[bin_start->vals[thread_id*bin_start->ncols + k]] = arr[j];
                ++bin_start->vals[thread_id*bin_start->ncols+k];  // Increment index for next number in bin
                break;
            }
        }
    }

    // Copy output to original array
    for (j=0; j<num; j++)
    {
        arr[j] = output[j];
    }

    free(output);
}

/**
 * @brief Perform radix sort
 *
 * @param arr Array to sort
 * @param num Number of entries in arr
 */
void RadixSort(
        uint32_t *arr,
        int num)
{

    int npes;
    int num_bits;
    int num_bins;
    int max_bins = pow(2, BIT_CHUNK_SIZE);  // Largest number of bins possible

    MPI_Comm_Size(MPI_COMM_WORLD, &npes);

    rs_matrix *bin_sizes, *bin_start;

    bin_sizes = rs_matrix_alloc(npes, max_bins);
    bin_start = rs_matrix_alloc(npes, max_bins);

    int max_iter = GetNumChunks(BIT_CHUNK_SIZE, NUM_BITS);

    for (int i=0; i<max_iter; i++)
    {

        num_bits = GetChunkSize(i, BIT_CHUNK_SIZE, NUM_BITS);
        num_bins = pow(2, num_bits);
        int mask = GetMask(num_bits);

        // NEED TO FINISH FROM HERE!!!!!!!!!!!!!!!!!!!!
    }
}

int main(
        int argc,
        char *argv[])
{

    return 0;
}
