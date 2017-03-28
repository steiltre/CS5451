#ifndef PAGERANK_RADIX_SORT_H
#define PAGERANK_RADIX_SORT_H

#include <stdint.h>

typedef uint64_t pr_int;


/**
 * @brief Compute mask for bitwise & operation
 *
 * @param num_bits Number of bits being sorted at once
 *
 * @return Mask for use with num_bits
 */
pr_int GetMask(
        pr_int const num_bits);


/**
 * @brief Sort an array of numbers according to num_bins sequential bits starting at ith bit
 *
 * @param arr Array to be sorted
 * @param i Bit to sort by
 * @param num Number of numbers in arr
 * @param num_bins Number of bins to place numbers in
 * @param bin_start Array giving starting location in arr for each bin
 * @param mask Mask to use for bitwise & operation
 */
void SortOnRadix(
        pr_int *arr,
        pr_int i,
        pr_int num,
        pr_int num_bins,
        pr_int *bin_start,
        pr_int mask);


/**
 * @brief Determine size of bins to use for radix sort
 *
 * @param arr Array being sorted
 * @param i Current bit being sorted
 * @param num Number of numbers in arr
 * @param num_bins Number of bins to place numbers into
 * @param mask Mask for bitwise & operation
 * @param[out] bin_sizes Array with size of each bin
 */
void DetermineBinSizes(
        pr_int *arr,
        pr_int i,
        pr_int num,
        pr_int num_bins,
        pr_int mask,
        pr_int *bin_sizes);


/**
 * @brief Gives bin starting indices from bin sizes array
 *
 * @param bin_sizes Array containing size of each bin
 * @param num_bins Number of bins in arrays
 * @param[out] bin_start Array of starting indices for bins
 */
void BinStartIndices(
        pr_int *bin_sizes,
        pr_int num_bins,
        pr_int *bin_start);


/**
 * @brief Perform radix sort
 *
 * @param arr Array to sort
 * @param num Number of entries in arr
 */
void radix_sort(
    pr_int *arr,
    pr_int num);

#endif
