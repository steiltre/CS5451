#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "rs_shared.h"

/**
 * @brief Number of bits to sort with at one time
 */
static int const BIT_CHUNK_SIZE = 1;

/**
 * @brief Width of integers being sorted
 */
static int const NUM_BITS = 32;

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
        int *bin_start,
        int mask,
        int nthreads)
{
    uint32_t *output;
    int j;

    output = (uint32_t *) malloc(num * sizeof(uint32_t));

    // Scan array putting values in order within local array for output
    #pragma omp parallel for schedule(static)
    for (j=0; j<num; j++)
    {
        int thread_id = omp_get_thread_num();

        for (int k=0; k<local_num_bins; k++)
        {
            if ( ((arr[j] >> i * BIT_CHUNK_SIZE) & mask) == k)
            {
                output[bin_start[thread_id*local_num_bins + k]] = arr[j];
                ++bin_start[thread_id*local_num_bins+k];  // Increment index for next number in bin
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
 * @brief Determine size of bins to use for radix sort
 *
 * @param arr Array being sorted
 * @param i Current bit being sorted
 * @param num Number of numbers in arr
 * @param local_num_bins Number of bins to place numbers into on each thread
 * @param mask Mask for bitwise & operation
 * @param nthreads Number of threads
 * @param[out] bin_sizes Array with size of each bin
 */
void DetermineBinSizes(
        uint32_t *arr,
        int i,
        int num,
        int local_num_bins,
        int mask,
        int nthreads,
        int *bin_sizes)
{
    // Scan array once to find where 0's end and 1's begin
    #pragma omp parallel for schedule(static)
    for (int j=0; j<num; j++)
    {
        int thread_id = omp_get_thread_num();
        for (int k=0; k<local_num_bins; k++)
        {
            if ( ((arr[j] >> i * BIT_CHUNK_SIZE) & mask) == k)
            {
                ++bin_sizes[thread_id*local_num_bins + k];
                #pragma omp atomic
                ++bin_sizes[nthreads*local_num_bins + k];
            }
        }
    }
}

/**
 * @brief Gives bin starting indices from bin sizes array
 *
 * @param bin_sizes Array containing size of each bin
 * @param num_bins Number of bins in arrays
 * @param nthreads Number of threads
 * @param[out] bin_start Array of starting indices for bins
 */
void BinStartIndices(
        int *bin_sizes,
        int local_num_bins,
        int nthreads,
        int *bin_start)
{
    bin_start[0] = 0;

    for (int i=1; i<local_num_bins; i++)
    {
        bin_start[i] = bin_start[i-1] + bin_sizes[nthreads*local_num_bins + i-1];
    }

    #pragma omp parallel for schedule(static)
    for (int j=0; j<local_num_bins; j++)
    {
        for (int k=1; k<nthreads; k++)
        {
            bin_start[k*local_num_bins + j] = bin_start[(k-1)*local_num_bins + j] + bin_sizes[(k-1)*local_num_bins + j];
        }
    }
}

/**
 * @brief Perform radix sort
 *
 * @param arr Array to sort
 * @param num Number of entries in arr
 * @param nthreads Number of threads
 */
void RadixSort(
        uint32_t *arr,
        int num,
        int nthreads)
{
    int num_bits;
    int num_bins;
    int max_bins = pow(2,BIT_CHUNK_SIZE);  // Largest number of bins possible
    int bin_sizes[max_bins * (nthreads + 1)];
    int bin_start[max_bins * nthreads];

    int max_iter = GetNumChunks( BIT_CHUNK_SIZE, NUM_BITS );

    for (int i=0; i<max_iter; i++)
    {

        num_bits = GetChunkSize( i, BIT_CHUNK_SIZE, NUM_BITS );
        num_bins = pow(2,num_bits);
        int mask = GetMask(num_bits);

        for (int j=0; j<num_bins * (nthreads + 1); j++)
        {
            bin_sizes[j] = 0;
        }

        DetermineBinSizes(arr, i, num, num_bins, mask, nthreads, bin_sizes);
        BinStartIndices(bin_sizes, num_bins, nthreads, bin_start);
        SortOnRadix(arr, i, num, num_bins, bin_start, mask, nthreads);
    }
}

/**
 * @brief Read numbers from file
 *
 * @param filename File containing points
 * @param[out] arr Array containing numbers
 * @param[out] num Number of numbers in arr
 */
void ReadFile(
        char const *filename,
        uint32_t **arr,
        int *num)
{
    FILE *fp;
    fp = fopen(filename, "r");
    fscanf(fp, "%d", num);

    (*arr) = (uint32_t *) malloc( *num * sizeof(uint32_t) );

    // Read numbers from data file
    for (int i=0; i<*num; i++)
    {
        fscanf(fp, "%u", arr[0]+i);       //NEED TO UNDERSTAND WHY arr[0]+i IS NEEDED HERE!!!!!
    }

}

int main( int argc, char *argv[] )
{
    char *infile = argv[1];
    int num_threads = atoi(argv[2]);
    char *outfile = argv[3];

    omp_set_num_threads(num_threads);

    int num;   // Number of numbers being sorted

    uint32_t *numbers;   // Array of numbers being sorted

    double start;

    ReadFile( infile, &numbers, &num );

    start = monotonic_seconds();

    RadixSort( numbers, num, num_threads );

    print_time(monotonic_seconds()-start);

    print_numbers( outfile, numbers, num );

    free(numbers);

    return 0;
}

