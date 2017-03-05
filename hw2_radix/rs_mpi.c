#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
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
 * @brief Implementation of modulo that works for negative dividends
 *
 * @param n Dividend
 * @param d Divisor
 *
 * @return r Remainder
 */
int mod (int n, int d)
{
    int r = n % d;
    if (r < 0)
    {
        r += d;
    }
    return r;
}

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
 * @brief Copy dense matrix
 *
 * @param matrix The matrix to copy
 *
 * @return mat_copy Copy of the original matrix
 */
rs_matrix * rs_matrix_copy(
        rs_matrix * matrix)
{
    rs_matrix * mat_copy;

    mat_copy = rs_matrix_alloc(matrix->nrows, matrix->ncols);

    for (int i=0; i<matrix->nrows * matrix->ncols; i++)
    {
        mat_copy->vals[i] = matrix->vals[i];
    }

    return mat_copy;
}

/**
 * @brief Transpose a dense matrix
 *
 * @param matrix The matrix to transpose
 */
void rs_matrix_transpose(
        rs_matrix * matrix)
{
    rs_matrix * temp_mat;
    temp_mat = rs_matrix_copy(matrix);

    for (int i=0; i<temp_mat->ncols; i++)
    {
        for (int j=0; j<temp_mat->nrows; j++)
        {
            matrix->vals[i*temp_mat->nrows + j] = temp_mat->vals[j*temp_mat->ncols + i];
        }
    }

    matrix->nrows = temp_mat->ncols;
    matrix->ncols = temp_mat->nrows;

    rs_matrix_free(temp_mat);
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

    rs_matrix_transpose(global_bin_end);  // Makes indexing easier in multiple places

    for (int j=1; j<global_bin_end->nrows; j++)
    {
        for (int i=0; i<global_bin_end->ncols; i++)
        {
            global_bin_end->vals[j*global_bin_end->ncols + i] += global_bin_end->vals[j*global_bin_end->ncols - 1];
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

rs_matrix *SendMatrix(
        rs_matrix *bin_end)
{
    int chunk_size;
    int l_index, u_index, l_process, u_process;
    int pid, npes;

    rs_matrix * send_mat;

    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    chunk_size = (int) ceil( ( (double) bin_end->vals[bin_end->ncols*bin_end->nrows-1]) / npes);

    send_mat = rs_matrix_alloc(bin_end->nrows, bin_end->ncols);

    for (int i=0; i<send_mat->nrows*send_mat->ncols; i++)
    {
        send_mat->vals[i] = 0;
    }

    int b = bin_end->vals[ bin_end->ncols * bin_end->nrows - 1 ];  // Use this to give "wrap-around" so first row isn't treated separately

    for (int i=0; i<send_mat->nrows; i++)
    {
        if (bin_end->vals[ i*bin_end->ncols + pid-1 ] == bin_end->vals[ i*bin_end->ncols + pid ])
        {
            continue;
        }

        l_index = mod( bin_end->vals[ mod( i*bin_end->ncols + pid-1, bin_end->ncols*bin_end->nrows ) ], b);
        u_index = bin_end->vals[ i*bin_end->ncols + pid ] - 1;

        if (l_index <= u_index)  // Bin is not empty
        {
            l_process = l_index / chunk_size;
            u_process = u_index / chunk_size;

            send_mat->vals[i*send_mat->ncols + l_process] += (l_process+1) * chunk_size - l_index;
            send_mat->vals[i*send_mat->ncols + u_process] += u_index - (l_process+1)*chunk_size + 1;
        }
    }

    return send_mat;
}

/**
 * @brief Determine where process is going to receive its numbers
 *
 * @param bin_end Global ending indices of bins on all processes
 *
 * @return incoming 1-dimensional matrix of locations to receive numbers from
 */
rs_matrix *RecvMatrix(
        rs_matrix *bin_end)
{
    int sum;
    int pos_sum = 0;  // Used to track change in sum when sum becomes positive
    int ideal_chunk_size, local_chunk_size;
    int pid, npes;

    rs_matrix * recv_mat;

    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    ideal_chunk_size = (int) ceil( ( (double) bin_end->vals[bin_end->ncols*bin_end->nrows-1]) / npes);
    int b = bin_end->vals[ bin_end->ncols * bin_end->nrows - 1];  // Last entry in matrix gives total number of numbers being sorted
    local_chunk_size = GetChunkSize( pid, ideal_chunk_size, b);

    recv_mat = rs_matrix_alloc(bin_end->nrows, bin_end->ncols);

    for (int i = 0; i<recv_mat->nrows*recv_mat->ncols; i++)
    {
        recv_mat->vals[i] = 0;
    }

    sum = -1*pid*ideal_chunk_size;

      for (int i=0; i<recv_mat->nrows; i++)
      {
          for (int j=0; j<recv_mat->ncols; j++)
          {
              sum += bin_end->vals[i *bin_end->ncols +j] - mod( bin_end->vals[ mod( i*bin_end->ncols + j-1, bin_end->ncols*bin_end->nrows ) ], b );

              if (sum > 0)
              {
                  if (sum < local_chunk_size)
                  {
                      recv_mat->vals[i *recv_mat->ncols + j] = sum - pos_sum;
                      pos_sum = sum;
                  }
                  else
                  {
                      recv_mat->vals[i*recv_mat->ncols + j] = local_chunk_size - pos_sum;
                      pos_sum = sum;
                      return recv_mat;
                  }
              }
          }
      }

    return recv_mat;
}

/**
 * @brief Shuffle values among processes for next sort
 *
 * @param arr Values on local process
 * @param num Number of values in arr
 * @param bin_end Global ending indices of local bins on all processes
 * @param local_bin_end Ending indices of bins in local array
 */
void ShuffleValues(
        uint32_t *arr,
        int num,
        rs_matrix *bin_end,
        rs_matrix *local_bin_end)
{
    rs_matrix * send_mat, * recv_mat;
    int npes, pid;
    int send_ind=0, recv_ind=0;
    MPI_Status status;

    uint32_t *temp_arr;

    temp_arr = (uint32_t *) malloc(num * sizeof(uint32_t));

    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    send_mat = SendMatrix(bin_end);
    recv_mat = RecvMatrix(bin_end);

    /* Go through each communication in order of process IDs to make sort stable and prevent deadlocks.
     * Blocking receives will ensure data is communicated in correct order (even when sends are buffered). */

    for (int i=0; i<send_mat->nrows; i++)
    {
        for (int j=0; j<npes; j++)
        {
            if (send_mat->vals[i*send_mat->ncols + j] > 0)  // Send message
            {
                if (j != pid)
                {
                    MPI_Send(&(arr[send_ind]), send_mat->vals[i*send_mat->ncols+j], MPI_UINT32_T, j, 0, MPI_COMM_WORLD);
                    send_ind += send_mat->vals[i*send_mat->ncols+j];
                }
            }

            if (j == pid)  // Receive messages from current bins
            {
                for (int k=0; k<npes; k++)
                {
                    if (recv_mat->vals[i*recv_mat->ncols + k] > 0)  // Receive message
                    {
                        if (k != pid)
                        {
                            MPI_Recv(&(temp_arr[recv_ind]), recv_mat->vals[i*recv_mat->ncols + k], MPI_UINT32_T, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                        }
                        else
                        {
                            memcpy( temp_arr + recv_ind, arr + send_ind, recv_mat->vals[i*recv_mat->ncols+k] * sizeof(uint32_t) );
                            send_ind += send_mat->vals[i*send_mat->ncols+j];
                        }
                        recv_ind += recv_mat->vals[i*recv_mat->ncols + k];
                    }
                }
            }
        }
    }

    /*
    for (int i=0; i<num; i++)
    {
        arr[i] = temp_arr[i];
    } */
    memcpy(arr, temp_arr, num*sizeof(uint32_t));
    free(temp_arr);
    rs_matrix_free(send_mat);
    rs_matrix_free(recv_mat);
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

    rs_matrix *copy_bin_end = rs_matrix_copy(bin_end);

    output = (uint32_t *) malloc(num * sizeof(uint32_t));

    for (j=0; j<num; j++)
    {
        output[j] = 0;
    }

    // Scan array putting values in order within local array for output
    for (j=num-1; j>=0; --j)
    {
        for (int k=0; k<local_num_bins; k++)
        {
            if ( ((arr[j] >> i * BIT_CHUNK_SIZE) & mask) == k)
            {
                --copy_bin_end->vals[k];  // Increment index for next number in bin
                output[copy_bin_end->vals[k]] = arr[j];
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
    rs_matrix_free(copy_bin_end);
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
    local_bin_end = rs_matrix_alloc(1, max_bins);

    int max_iter = GetNumChunks(BIT_CHUNK_SIZE, NUM_BITS);
    //int max_iter = 1;

    for (int i=0; i<max_iter; i++)
    {

        num_bits = GetChunkSize(i, BIT_CHUNK_SIZE, NUM_BITS);
        num_bins = pow(2, num_bits);
        int mask = GetMask(num_bits);
        global_bin_end = rs_matrix_alloc(npes, max_bins);

        for (int j=0; j<max_bins; j++)
        {
            bin_sizes->vals[j] = 0;
        }

        DetermineBinSizes(arr, i, num, num_bins, mask, bin_sizes);
        BinEndIndices(bin_sizes, num_bins, global_bin_end, local_bin_end);
        SortOnRadix(arr, i, num, num_bins, local_bin_end, mask);
        ShuffleValues(arr, num, global_bin_end, local_bin_end);

        rs_matrix_free(global_bin_end);
    }

    rs_matrix_free(bin_sizes);
    //rs_matrix_free(global_bin_end);
    rs_matrix_free(local_bin_end);
}

/**
 * @brief Use a single process to read array from file and distribute to other processes
 *
 * @param filename Name of file containing numbers to sort
 * @param arr Array to place numbers into
 * @param num Number of numbers allocated to local process
 * @param total_num Number of total numbers being sorted
 */
void ReadFile(
        char const *filename,
        uint32_t **arr,
        int *num,
        int *total_num)
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
        fscanf(fp, "%d", total_num);

        *num = *total_num;

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

/**
 * @brief Gather output and print to file
 *
 * @param filename Name of output file
 * @param arr Array of numbers for output
 * @param num Local number of numbers being sorted
 * @param total_num Total number of numbers being sorted
 */
void Output(
        char const *filename,
        uint32_t *arr,
        int num,
        int total_num)
{
    int pid, npes;

    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    if (pid == 0)
    {
        int chunk_size;
        int ideal_chunk_size = (int) ceil( ((double) total_num)/npes );
        MPI_Status status;

        for (int i=1; i<npes; i++)
        {
            chunk_size = GetChunkSize(i, ideal_chunk_size, total_num);

            MPI_Recv(&(arr[i*(ideal_chunk_size)]), chunk_size, MPI_UINT32_T, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        }

        print_numbers(filename, arr, total_num);
    }
    else
    {
        MPI_Send(arr, num, MPI_UINT32_T, 0, 0, MPI_COMM_WORLD);
    }
}

int main(
        int argc,
        char *argv[])
{

    MPI_Init(&argc, &argv);

    char *infile = argv[1];
    char *outfile = argv[2];

    double start, end;

    uint32_t *numbers;  // Array of numbers being sorted
    int num;  // Number of numbers being sorted locally
    int total_num;  // Number of numbers being sorted globally
    int pid;

    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    ReadFile(infile, &numbers, &num, &total_num);

    start = MPI_Wtime();

    RadixSort(numbers, num);

    end = MPI_Wtime();

    if (pid == 0)
        print_time(end-start);

    Output(outfile, numbers, num, total_num);

    free(numbers);

    MPI_Finalize();

    return 0;
}
