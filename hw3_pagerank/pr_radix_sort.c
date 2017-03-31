#include <math.h>
#include <stdlib.h>

#include "pr_utils.h"

/**
 * @brief Number of bits to sort with at one time
 */
static pr_int const BIT_CHUNK_SIZE = 4;

/**
 * @brief Width of integers being sorted
 */
static pr_int const NUM_BITS = 32;

pr_int GetMask(
        int const num_bits)
{
    pr_int mask = 0;
    for (int j=0; j<num_bits; j++)
    {
        mask += pow(2,j);
    }

    return mask;
}

void SortOnRadix(
        pr_int *arr,
        int i,
        int num,
        int num_bins,
        int *bin_start,
        pr_int mask)
{
    int * curr_ind;
    pr_int * output;
    int j;

    output = malloc( num * sizeof(*output) );
    curr_ind = malloc( num_bins * sizeof(*curr_ind) );

    // Initialize curr_ind to starting location of bins
    for (j=0; j<num_bins; j++)
    {
        curr_ind[j] = bin_start[j];
    }

    // Scan array putting values in order for output
    for (j=0; j<num; j++)
    {
        for (int k=0; k<num_bins; k++)
        {
            if ( ((arr[j] >> i * BIT_CHUNK_SIZE) & mask) == k)
            {
                output[curr_ind[k]] = arr[j];
                curr_ind[k] += 1;
            }
        }
    }

    // Copy output to original array
    for (j=0; j<num; j++)
    {
        arr[j] = output[j];
    }

    free(output);
    free(curr_ind);
}

void DetermineBinSizes(
        pr_int *arr,
        pr_int i,
        int num,
        int num_bins,
        pr_int mask,
        int *bin_sizes)
{
    // Scan array once to find where 0's end and 1's begin
    for (int j=0; j<num; j++)
    {
        for (int k=0; k<num_bins; k++)
        {
            if ( ((arr[j] >> i * BIT_CHUNK_SIZE) & mask) == k)
            {
                bin_sizes[k] += 1;
            }
        }
    }
}

void BinStartIndices(
        int *bin_sizes,
        int num_bins,
        int *bin_start)
{
    bin_start[0] = 0;

    for (int i=1; i<num_bins; i++)
    {
        bin_start[i] = bin_start[i-1] + bin_sizes[i-1];
    }
}

void radix_sort(pr_int *arr, int num)
{
    int num_bits;
    int num_bins;
    int max_bins = pow(2,BIT_CHUNK_SIZE);  // Largest number of bins possible
    int * bin_sizes;
    int * bin_start;

    bin_sizes = malloc( max_bins * sizeof(*bin_sizes) );
    bin_start = malloc( max_bins * sizeof(*bin_start) );

    int max_iter = GetNumChunks( BIT_CHUNK_SIZE, NUM_BITS );

    for (int i=0; i<max_iter; i++)
    {

        num_bits = GetChunkSize( i, BIT_CHUNK_SIZE, NUM_BITS );
        num_bins = pow(2,num_bits);
        pr_int mask = GetMask(num_bits);

        for (int j=0; j<num_bins; j++)
        {
            bin_sizes[j] = 0;
        }

        DetermineBinSizes(arr, i, num, num_bins, mask, bin_sizes);
        BinStartIndices(bin_sizes, num_bins, bin_start);
        SortOnRadix(arr, i, num, num_bins, bin_start, mask);
    }

    free(bin_sizes);
    free(bin_start);
}
