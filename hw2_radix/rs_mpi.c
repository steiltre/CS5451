#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
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
 * @brief Print array to log file for debugging
 *
 * @param arr Array to print
 * @param num Number of entries to print
 */
//void LogArray(
//        uint32_t *arr,
//        int num)
//{
//    int pid;
//    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
//
//    char pid_str[2];
//    snprintf(pid_str, 12, "%d", pid);
//
//    char log_name[20] = "process";
//    char ext_str[4] = ".txt";
//
//    FILE *fp;
//
//    strcat(log_name, pid_str);
//    strcat(log_name, ext_str);
//
//    fp = fopen(log_name, "w");
//
//    for (int i=0; i<num; i++)
//    {
//        fprintf(fp, "%u \n", arr[i]);
//    }
//
//    fclose(fp);
//}


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
        rs_matrix *bin_sizes)
{
    // Scan array once to find where 0's end and 1's begin
    for (int j=0; j<num; j++)
    {
        for (int k=0; k<local_num_bins; k++)
        {
            if ( ((arr[j] >> i * BIT_CHUNK_SIZE) & mask) == k)
            {
                ++bin_sizes->vals[k];
                break;
            }
        }
    }
}


/**
 * @brief Gives bin ending indices from bin sizes array
 *
 * @param bin_sizes Array containing size of each bin
 * @param num_bins Number of bins in arrays
 * @param nthreads Number of threads
 * @param[out] bin_end Array of ending indices for bins
 */
void BinEndIndices(
        rs_matrix *bin_sizes,
        int local_num_bins,
        rs_matrix *global_bin_end,
        rs_matrix *local_bin_end)
{
    int pid;
    int *temp_bin_end;

    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    temp_bin_end = (int *) malloc(local_num_bins * sizeof(int));

    /*
     * Calculate global_bin_end for shuffling values after local sort
     */

    // Use prefix sum to get ending location of local bins in global array
    MPI_Scan(bin_sizes->vals, temp_bin_end, local_num_bins, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // Use Allgather to get ending index of each local bin to every other process (for shuffling values)
    MPI_Allgather(temp_bin_end, local_num_bins, MPI_INT, global_bin_end->vals, local_num_bins, MPI_INT, MPI_COMM_WORLD);

    // Prefix sum isn't quite global indices
    for (int i=1; i<global_bin_end->ncols; i++)
    {
        global_bin_end->vals[(global_bin_end->ncols-1)*global_bin_end->nrows + i] += global_bin_end->vals[(global_bin_end->ncols-1)*global_bin_end->nrows + i-1];
    }

    for (int j=0; j<global_bin_end->nrows-1; j++)
    {
        for (int i=1; i<global_bin_end->ncols; i++)
        {
            global_bin_end->vals[global_bin_end->ncols*j + i] += global_bin_end->vals[(global_bin_end->ncols-1)*global_bin_end->nrows + i-1];
        }
    }

    /*
     * Calculate local_bin_end for sorting on individual process
     */
    local_bin_end->vals[0] = bin_sizes->vals[0];
    for (int i=1; i<local_num_bins; i++)
    {
        local_bin_end->vals[i] = local_bin_end->vals[i-1] + bin_sizes->vals[i];
    }
}

/**
 * @brief Determine where process is going to send its numbers
 *
 * @param bin_end Global ending indices of bins on all processes
 *
 * @return outgoing 1-dimensional matrix of locations to send numbers to
 */
rs_matrix *SendList(
        rs_matrix *bin_end)
{
    int i;
    int send_to[pow(2, BIT_CHUNK_SIZE + 1)];  // Array with size of maximum number of locations to send bits
    int indices[pow(2, BIT_CHUNK_SIZE + 1)];  // Array to store beginning and ending indices for local bins

    int pid;
    Get_Comm_rank(MPI_COMM_WORLD, pid);

    int b = bin_end[ bin_end->ncols * bin_end->nrows - 1 ];  // Use this to give "wrap-around" so first row isn't treated separately

    if (pid == 0)
    {
        indices[i] = 0;
        indices[i] = bin_end->vals[ bin_end->ncols + 1 ] - 1;
        i++;

        for (int j=1; j<bin_end->ncols; j++)
        {
        }
}

/**
 * @brief Determine where process is going to receive its numbers
 *
 * @param bin_end Global ending indices of bins on all processes
 *
 * @return incoming 1-dimensional matrix of locations to receive numbers from
 */
rs_matrix *RecvList(
        rs_matrix *bin_end)
{
}

/**
 * @brief Shuffle values among processes for next sort
 *
 * @param arr Values on local process
 * @param num Number of values in arr
 * @param bin_end Global ending indices of local bins on all processes
 */
void ShuffleValues(
        uint32_t *arr,
        int num,
        rs_matrix *bin_end)
{
}

/**
 * @brief Sort an array of numbers according to num_bins sequential bits starting at ith bit
 *
 * @param arr Array to be sorted
 * @param i Bit to sort by
 * @param num Number of numbers in arr
 * @param local_num_bins Number of bins to place numbers in
 * @param bin_end Array giving ending location in arr for each bin
 * @param mask Mask to use for bitwise & operation
 */
void SortOnRadix(
        uint32_t *arr,
        int i,
        int num,
        int local_num_bins,
        rs_matrix *bin_end,
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
        for (int k=0; k<local_num_bins; k++)
        {
            if ( ((arr[j] >> i * BIT_CHUNK_SIZE) & mask) == k)
            {
                output[bin_end->vals[k]] = arr[j];
                ++bin_end->vals[k];  // Increment index for next number in bin
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

    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    rs_matrix *bin_sizes, *global_bin_end, *local_bin_end;

    bin_sizes = rs_matrix_alloc(1, max_bins);
    global_bin_end = rs_matrix_alloc(npes, max_bins);
    local_bin_end = rs_matrix_alloc(1, max_bins);

    //int max_iter = GetNumChunks(BIT_CHUNK_SIZE, NUM_BITS);
    int max_iter = 1;

    for (int i=0; i<max_iter; i++)
    {

        num_bits = GetChunkSize(i, BIT_CHUNK_SIZE, NUM_BITS);
        num_bins = pow(2, num_bits);
        int mask = GetMask(num_bits);

        for (int j=0; j<max_bins; j++)
        {
            bin_sizes->vals[j] = 0;
        }

        DetermineBinSizes(arr, i, num, num_bins, mask, bin_sizes);
        BinEndIndices(bin_sizes, num_bins, global_bin_end, local_bin_end);
        SortOnRadix(arr, i, num, num_bins, local_bin_end, mask);
        ShuffleValues(arr, num, global_bin_end);
        // NEED TO FINISH FROM HERE!!!!!!!!!!!!!!!!!!!!
    }

    free(bin_sizes);
}

/**
 * @brief Use a single process to read array from file and distribute to other processes
 *
 * @param filename Name of file containing numbers to sort
 * @param arr Array to place numbers into
 * @param num Number of numbers allocated to local process
 */
void ReadFile(
        char const *filename,
        uint32_t **arr,
        int *num)
{
    int pid, npes, chunk_size;
    MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    if (pid == 0)
    {
        int send_num;  // Number of numbers to send to process
        FILE *fp;
        fp = fopen(filename, "r");
        fscanf(fp, "%d", num);

        for (int i=1; i<npes; i++)
        {
            MPI_Send(num, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        *arr = (uint32_t *) malloc(*num *sizeof(uint32_t) );

        chunk_size = (int) ceil( ((double) *num)/npes );

        for (int i=0; i<*num; i++)
        {
            fscanf(fp, "%u", arr[0]+i);
        }

        fclose(fp);

        for (int i=1; i<npes; i++)
        {
            send_num = GetChunkSize(i, chunk_size, *num);
            MPI_Send(*arr+i*chunk_size, send_num, MPI_UINT32_T, i, 0, MPI_COMM_WORLD);
        }

        *num = GetChunkSize(0, chunk_size, *num);

//        LogArray(arr, *num);
    }
    else
    {
        MPI_Recv(num, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        chunk_size = (int) ceil( ((double) *num)/npes );

        *num = GetChunkSize(pid, chunk_size, *num);

        *arr = (uint32_t *) malloc(*num * sizeof(uint32_t) );

        MPI_Recv(*arr, *num, MPI_UINT32_T, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

//        LogArray(arr, *num);
    }

}

int main(
        int argc,
        char *argv[])
{

    MPI_Init(&argc, &argv);

    char *infile = argv[1];
    char *outfile = argv[2];

    uint32_t *numbers;  // Array of numbers being sorted
    int num;  // Number of numbers being sorted

    ReadFile(infile, &numbers, &num);

    RadixSort(numbers, num);

    //free(numbers);

    MPI_Finalize();

    return 0;
}
