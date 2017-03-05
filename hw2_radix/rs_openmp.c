#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

/**
 * @brief Number of bits to sort with at one time
 */
static int const BIT_CHUNK_SIZE = 4;

/**
 * @brief Width of integers being sorted
 */
static int const NUM_BITS = 32;

/**
 * @brief Length of rows with padding
 */
static int const PAD_WIDTH = 64;

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
 * * @brief Write an array of integers to a file.
 * *
 * * @param filename The name of the file to write to.
 * * @param numbers The array of numbers.
 * * @param nnumbers How many numbers to write.
 * */
static void print_numbers(
            char const * const filename,
                uint32_t const * const numbers,
                    uint32_t const nnumbers)
{
      FILE * fout;

        /* open file */
        if((fout = fopen(filename, "w")) == NULL) {
                fprintf(stderr, "error opening '%s'\n", filename);
                    abort();
                      }

          /* write the header */
          fprintf(fout, "%d\n", nnumbers);

            /* write numbers to fout */
            for(uint32_t i = 0; i < nnumbers; ++i) {
                    fprintf(fout, "%d\n", numbers[i]);
                      }

              fclose(fout);
}

/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 200809L
#include <time.h>

/* OSX timer includes */
#ifdef __MACH__
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#endif

/**
 * * @brief Return the number of seconds since an unspecified time (e.g., Unix
 * *        epoch). This is accomplished with a high-resolution monotonic timer,
 * *        suitable for performance timing.
 * *
 * * @return The number of seconds.
 * */
static inline double monotonic_seconds()
{
#ifdef __MACH__
      /* OSX */
      static mach_timebase_info_data_t info;
        static double seconds_per_unit;
          if(seconds_per_unit == 0) {
                  mach_timebase_info(&info);
                      seconds_per_unit = (info.numer / info.denom) / 1e9;
                        }
            return seconds_per_unit * mach_absolute_time();
#else
              /* Linux systems */
              struct timespec ts;
                clock_gettime(CLOCK_MONOTONIC, &ts);
                  return ts.tv_sec + ts.tv_nsec * 1e-9;
#endif
}

/**
 * * @brief Output the seconds elapsed while sorting. This excludes input and
 * *        output time. This should be wallclock time, not CPU time.
 * *
 * * @param seconds Seconds spent sorting.
 * */
static void print_time(
            double const seconds)
{
      printf("Sort Time: %0.04fs\n", seconds);
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

    #pragma omp parallel for schedule(static)
    for (j=0; j<num; j++)
    {
        output[j] = 0;
    }

    // Scan array putting values in order within local array for output
    #pragma omp parallel for schedule(static)
    for (j=0; j<num; j++)
    {
        int thread_id = omp_get_thread_num();

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
    #pragma omp parallel for schedule(static)
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
        rs_matrix *bin_sizes)
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
                ++bin_sizes->vals[thread_id*bin_sizes->ncols + k];
                break;
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
        rs_matrix *bin_sizes,
        int local_num_bins,
        rs_matrix *bin_start)
{

    for (int i=0; i<local_num_bins; i++)
    {
        for (int j=0; j<bin_start->nrows; j++)
        {
            bin_start->vals[j*bin_start->ncols + i] = 0;
        }
    }

    for (int i=0; i<local_num_bins; i++)
    {
        if (i == 0)
        {
            bin_start->vals[i] = 0;
        }
        else
        {
            bin_start->vals[i] = bin_start->vals[(bin_start->nrows-1)*bin_start->ncols + i-1] + bin_sizes->vals[(bin_sizes->nrows-1)*bin_sizes->ncols + i-1];
        }

        for (int j=1; j<bin_sizes->nrows; j++)
        {
            bin_start->vals[j*bin_start->ncols + i] += bin_start->vals[(j-1)*bin_start->ncols + i] + bin_sizes->vals[(j-1)*bin_sizes->ncols + i];
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

    rs_matrix *bin_sizes, *bin_start;

    // Use num_bins columns if num_bins > PAD_WIDTH, otherwise pad with extra columns
    if (max_bins > PAD_WIDTH)
    {
        bin_sizes = rs_matrix_alloc(nthreads, max_bins);
        bin_start = rs_matrix_alloc(nthreads, max_bins);
    }
    else
    {
        bin_sizes = rs_matrix_alloc(nthreads, PAD_WIDTH);
        bin_start = rs_matrix_alloc(nthreads, PAD_WIDTH);
    }

    int max_iter = GetNumChunks( BIT_CHUNK_SIZE, NUM_BITS );

    for (int i=0; i<max_iter; i++)
    {

        num_bits = GetChunkSize( i, BIT_CHUNK_SIZE, NUM_BITS );
        num_bins = pow(2,num_bits);
        int mask = GetMask(num_bits);

        #pragma omp parallel for schedule(static)
        for (int j=0; j<bin_sizes->nrows; j++)
        {
            for (int k=0; k<max_bins; k++)
            {
                bin_sizes->vals[j*bin_sizes->ncols + k] = 0;
            }
        }

        DetermineBinSizes(arr, i, num, num_bins, mask, bin_sizes);
        BinStartIndices(bin_sizes, num_bins, bin_start);
        SortOnRadix(arr, i, num, num_bins, bin_start, mask);
    }

    free(bin_sizes);
    free(bin_start);
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

    #pragma omp parallel for schedule(static)
    for (int i=0; i<*num; i++)
    {
        (*arr)[i] = 0;
    }

    // Read numbers from data file
    for (int i=0; i<*num; i++)
    {
        fscanf(fp, "%u", arr[0]+i);       //NEED TO UNDERSTAND WHY arr[0]+i IS NEEDED HERE!!!!!
    }

    fclose(fp);

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

