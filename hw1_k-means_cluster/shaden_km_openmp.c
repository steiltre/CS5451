/**
* @file km_openmp.c
* @brief A multi-threaded k-means clustering code using OpenMP. 
* @author Shaden Smith <shaden@cs.umn.edu>
* @version 0.1.0
* @date 2017-01-31
*/


/******************************************************************************
 * Includes
 *****************************************************************************/

/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define NDEBUG /* disable `assert()` for performance runs */
#include <assert.h>

/* OSX timer includes */
#ifdef __MACH__
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#endif

/**
* @brief When assigning a point, how many centroids do we examine at a time?
*/
static int const CLUST_CHUNK_SIZE = 256;


/**
* @brief How many points do we assign at a time?
*/
static int const POINT_CHUNK_SIZE = 2048;




/******************************************************************************
 * Types
 *****************************************************************************/


/**
* @brief Just a dense matrix structure.
*/
typedef struct
{
  int nrows;
  int ncols;
  double * vals;
} km_matrix;






/******************************************************************************
 * Utility functions
 *****************************************************************************/


/**
* @brief Return the number of seconds since an unspecified time (e.g., Unix
*        epoch). This is accomplished with a high-resolution monotonic timer,
*        suitable for performance timing.
*
* @return The number of seconds.
*/
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
* @brief A wrapper around `posix_memalign()` to get aligned memory.
*
* @param bytes How many bytes to allocate?
*
* @return Allocated memory.
*/
void * km_malloc(
    size_t const bytes)
{
  void * ptr;
  int success = posix_memalign(&ptr, 64, bytes);
  if(success != 0) {
    fprintf(stderr, "ERROR: posix_memalign() returned %d.\n", success);
    exit(1);
  }
  return ptr;
}


/**
* @brief Free memory allocated by `km_malloc()`.
*
* @param ptr The pointer to free.
*/
void km_free(
    void * ptr)
{
  free(ptr);
}



/**
* @brief Allocate a dense matrix.
*
* @param nrows The number of rows in the matrix.
* @param ncols The number of columns in the matrix.
*
* @return The allocated matrix.
*/
km_matrix * km_matrix_alloc(
    int const nrows,
    int const ncols)
{
  km_matrix * matrix = km_malloc(sizeof(*matrix));

  matrix->nrows = nrows;
  matrix->ncols = ncols;
  matrix->vals = km_malloc(nrows * ncols * sizeof(*matrix->vals));

  return matrix;
}


/**
* @brief Free memory allocated by `km_matrix_alloc()`.
*
* @param matrix The matrix to free.
*/
void km_matrix_free(
    km_matrix * matrix)
{
  if(matrix == NULL) {
    return;
  }

  km_free(matrix->vals);
  km_free(matrix);
}


/**
* @brief Compute the number of chunks to span a set of items.
*
* @param chunk_size The size of the chunks.
* @param num_items The number of items.
*
* @return  The number of chunks. The last chunk may be smaller than chunk_size.
*/
static inline int get_num_chunks(
    int const chunk_size,
    int const num_items)
{
  int num_chunks = num_items / chunk_size;
  if(num_items % chunk_size > 0) {
    ++num_chunks;
  }
  return num_chunks;
}



/**
* @brief Compute the number of rows to process in a chunk. The last chunk may
*        be smaller than the rest.
*
* @param chunk_id The chunk we are processing.
* @param chunk_size The ideal size of the chunks.
* @param num_items The number of total items.
*
* @return How many items in the chunk to process.
*/
static inline int get_chunk_end(
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
* @brief Write a matrix to a file.
*
* @param matrix The matrix to write.
* @param fname The filename to write to.
*/
void write_matrix(
    km_matrix const * const matrix,
    char const * const fname)
{
  FILE * fout;
  if((fout = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "unable to open '%s' for writing.\n", fname);
    exit(EXIT_FAILURE);
  }

  fprintf(fout, "%d %d\n", matrix->nrows, matrix->ncols);
  for(int i=0; i < matrix->nrows; ++i) {
    for(int j=0; j < matrix->ncols; ++j) {
      fprintf(fout, "%0.3f ", matrix->vals[j + (i*matrix->ncols)]);
    }
    fprintf(fout, "\n");
  }
  fclose(fout);
}



/**
* @brief Parse a file into a row-major matrix of doubles. The first two entries
*        in the file should specify the number of rows and columns,
*        respectively.
*
* @param fname The name of the file to read.
*
* @return The matrix data.
*/

km_matrix * read_matrix(
    char const * const fname)
{
  /* Open file */
  FILE * fin;
  if((fin = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "unable to open '%s' for reading.\n", fname);
    exit(EXIT_FAILURE);
  }

  int nrows;
  int ncols;

  /* Read dimensions.*/
  fscanf(fin, "%d", &nrows);
  fscanf(fin, "%d", &ncols);

  km_matrix * matrix = km_matrix_alloc(nrows, ncols);

  /* first-touch rows -- useful for NUMA systems */
  #pragma omp parallel for schedule(static)
  for(int i=0; i < nrows; ++i) {
    for(int j=0; j < ncols; ++j) {
      matrix->vals[j + (i*ncols)] = 0.;
    }
  }

  /* Read in row-major matrix. */
  for(size_t i=0; i < nrows * ncols; ++i) {
    fscanf(fin, "%lf", &(matrix->vals[i]));
  }

  fclose(fin);

  return matrix;
}


int rand_int()
{
  return abs((rand() + RAND_MAX) + rand());
}




/******************************************************************************
 * Clustering functions
 *****************************************************************************/

/**
* @brief Compute the Euclidean distance between two vectors.
*
* @param x The first vector.
* @param y The second vector.
* @param N The lengths of the vectors.
*
* @return The distance between x and y.
*/
static inline double vector_distance(
    double const * const restrict x,
    double const * const restrict y,
    int const N)
{
  double dist = 0.;
  for(int i=0; i < N; ++i) {
    dist += (y[i] - x[i]) * (y[i] - x[i]);
  }
  return sqrt(dist);
}




/**
* @brief Initialize centroids for k-means clustering. Each centroid is just a
*        randomly-chosen point in the data.
*
*        XXX: This function selects points *with* replacement. If you are
*        unlucky, or want many clusters, you will have identical centroids.
*
* @param points The points to cluster. We will select 'nclusters' of these.
* @param nclusters The number of clusters to initialize.
*
* @return  An (nclusters x dim) matrix of initial centroids.
*/
km_matrix * init_centroids(
    km_matrix const * const points,
    int const nclusters)
{
  km_matrix * centroids = km_matrix_alloc(nclusters, points->ncols);

  /* save some typing */
  int const npoints = points->nrows; /* only for random assignment */
  int const dim     = points->ncols;

  for(int c=0; c < nclusters; ++c) {
//#define SELECT_RAND_POINTS
#ifdef SELECT_RAND_POINTS
    int const point_id = rand_int() % npoints;
#else
    int const point_id = c;
#endif

    /* Copy row 'point_id' into centroid matrix. */
    double       * const restrict centroid_row = centroids->vals + (c * dim);
    double const * const restrict point_row = points->vals + (point_id * dim);
    for(int d=0; d < dim; ++d) {
      centroid_row[d] = point_row[d];
    }
  }

  return centroids;
}



/**
* @brief Compute the assignment of each point to a cluster. This performs the
* assignments in a parallel, cache-friendly way. When the centroid matrix does
* not fit in cache, assigning just one point at a time will force the entire
* centroids matrix to be fetched from RAM each point. Instead, we break the
* data points and centroids up into chunks, and we process chunks at a time.
* This allows us to tune how much of the data and centroid matrices we touch at
* a time.
*
* @param points The points to cluster.
* @param centroids The centroids representing each cluster.
* @param[out] assignments The array of cluster assignments.
*
* @return How many points changed cluster assignment.
*/
int compute_assignments(
    km_matrix const * const points,
    km_matrix const * const centroids,
    int * const restrict assignments)
{
  double const start = monotonic_seconds();

  /* save some typing */
  int const dim = points->ncols;

  int    min_cent[POINT_CHUNK_SIZE];
  double min_dist[POINT_CHUNK_SIZE];

  /* get number of point and cluster chunks */
  int const num_point_chunks = get_num_chunks(POINT_CHUNK_SIZE, points->nrows);
  int const num_clust_chunks = get_num_chunks(CLUST_CHUNK_SIZE, centroids->nrows);

  /* the number of points which changed cluster */
  int nchanged = 0;

  #pragma omp parallel reduction(+: nchanged)
  for(int pc=0; pc < num_point_chunks; ++pc) {
    int const i_stop = get_chunk_end(pc, POINT_CHUNK_SIZE, points->nrows);

    /* initialize chunk data */
    #pragma omp for schedule(static)
    for(int i=0; i < i_stop; ++i) {
      min_cent[i] = -1;
      min_dist[i] = DBL_MAX;
    }

    /* foreach cluster chunk */
    for(int cc=0; cc < num_clust_chunks; ++cc) {
      int const c_stop = get_chunk_end(cc, CLUST_CHUNK_SIZE, centroids->nrows);

      /* foreach point in chunk */
      #pragma omp for schedule(static)
      for(int i=0; i < i_stop; ++i) {
        int const i_global = i + (pc * POINT_CHUNK_SIZE);
        assert(i_global < points->nrows);
        double const * const point_row = points->vals + (i_global * dim);

        /* foreach centroid in chunk */
        for(int c=0; c < c_stop; ++c) {
          int const c_global = c + (cc * CLUST_CHUNK_SIZE);
          assert(c_global < centroids->nrows);
          double const * const cent_row = centroids->vals + (c_global * dim);

          double const dist = vector_distance(point_row, cent_row, dim);
          if(dist < min_dist[i]) {
            min_dist[i] = dist;
            min_cent[i] = c_global;
          }
        } /* foreach centroid in chunk */
      } /* foreach point in chunk */
    } /* foreach centroid chunk */

    /* update cluster assignment */
    #pragma omp for schedule(static)
    for(int i=0; i < i_stop; ++i) {
      int const i_global = i + (pc * POINT_CHUNK_SIZE);
      /* if a change was made, we haven't converged */
      if(min_cent[i] != assignments[i_global]) {
        ++nchanged;
      }
      assignments[i_global] = min_cent[i];
    }
  } /* foreach point chunk */

  double const elapsed = monotonic_seconds() - start;
  printf("    assignments: %0.3fs\n", elapsed);

  return nchanged;
}



/**
* @brief Update the centroids representing each cluster. This function relies
* on atomic hardware to work efficiently. Another option is to allocate
* thread-local centroids, compute new centroids totally in parallel, then
* perform a parallel reduction.
*
* @param points The points we are clustering.
* @param assignments The assignments of points to clusters.
* @param[out] centroids The centroids to update.
*/
void update_centroids(
    km_matrix const * const points,
    int const * const restrict assignments,
    km_matrix * const centroids)
{
  double const start =  monotonic_seconds();

  /* reset centroids */
  #pragma omp parallel for schedule(static)
  for(int x=0; x < centroids->nrows * centroids->ncols; ++x) {
    centroids->vals[x] = 0.;
  }

  int * cluster_sizes = km_malloc(centroids->nrows * sizeof(*cluster_sizes));
  #pragma omp parallel for schedule(static)
  for(int c=0; c < centroids->nrows; ++c) {
    cluster_sizes[c] = 0;
  }

  int const dim = centroids->ncols;

  /* push updates to centroids */
  #pragma omp parallel for schedule(static)
  for(int i=0; i < points->nrows; ++i) {
    int const c_id = assignments[i];
    #pragma omp atomic
    ++cluster_sizes[c_id];

    double const * const restrict point_row = points->vals + (i * dim);
    double       * const restrict cent_row = centroids->vals + (c_id * dim);

    /* add point to centroid */
    for(int d=0; d < dim; ++d) {
      #pragma omp atomic
      cent_row[d] += point_row[d];
    }
  }

  /* average points to get new centroids */
  #pragma omp parallel for schedule(static)
  for(int c=0; c < centroids->nrows; ++c) {
    double * const restrict cent_row = centroids->vals + (c * dim);
    for(int d=0; d < dim; ++d) {
      cent_row[d] /= cluster_sizes[c];
    }
  }

  double const elapsed = monotonic_seconds() - start;
  printf("    update: %0.3fs\n", elapsed);

  km_free(cluster_sizes);
}


/**
* @brief Perform k-means clustering.
*
* @param points
* @param npoints
* @param dim
* @param nclusters
*
* @return
*/
int * kmeans(
    km_matrix const * const points,
    int const nclusters)
{
  int const npoints = points->nrows;

  int * assignments = km_malloc(npoints * sizeof(*assignments));
  km_matrix * centroids = init_centroids(points, nclusters);

  double time_start = monotonic_seconds();

  /* compute initial assignment of points */
  #pragma omp parallel for schedule(static)
  for(int p=0; p < npoints; ++p) {
    assignments[p] = 0;
  }
  compute_assignments(points, centroids, assignments);

  int num_its = 0;
  for(int i=0; i < 20; ++i) {
    double start = monotonic_seconds();

    update_centroids(points, assignments, centroids);
    int nchanged = compute_assignments(points, centroids, assignments);

    ++num_its;
    double iter_time = monotonic_seconds() - start;
    printf("  iter: %d (%d moved) (%0.3fs)\n", num_its, nchanged, iter_time);

    if(nchanged == 0) {
      printf("Clustering converged in %d iterations.\n", num_its);
      break;
    }
  }

  double elapsed_time = monotonic_seconds() - time_start;

  printf("k-means: (%0.3fs)\n", elapsed_time);

  write_matrix(centroids, "centroids.txt");

  km_matrix_free(centroids);
  return assignments;
}


/******************************************************************************
 * Entry point
 *****************************************************************************/

int main(
    int argc,
    char * * argv)
{
  if(argc < 4) {
    printf("usage: %s <data> <num_clusters> <num_threads>\n", argv[0]);
    return EXIT_FAILURE;
  }

  char * data_fname = argv[1];

  /* Parse input */
  km_matrix * points = read_matrix(data_fname);
  int nclusters = atoi(argv[2]);

  /* set number of threads */
  int nthreads = atoi(argv[3]);
  omp_set_num_threads(nthreads);

  int * assignments = kmeans(points, nclusters);

  /* optionally write assignments */
  FILE * fout;
  if((fout = fopen("clusters.txt", "w")) == NULL) {
    fprintf(stderr, "unable to open 'clusters.txt' for writing.\n");
    exit(EXIT_FAILURE);
  }

  /* write cluster assignments */
  for(int i=0; i < points->nrows; ++i) {
    fprintf(fout, "%d\n", assignments[i]);
  }

  fclose(fout);

  /* cleanup */
  km_matrix_free(points);
  km_free(assignments);

  return EXIT_SUCCESS;
}



