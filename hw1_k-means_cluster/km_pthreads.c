#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "timing.h"

void *FindClosestCentroid( void *input_arg );
void *DetermineCentroids( void *input_arg );
void WriteOutput( int *closest_centroid, double *centroids, int dim, int num_clusters, int num_points );
static void print_time( double const seconds );
pthread_mutex_t closest_centroid_lock;
pthread_mutex_t determine_centroid_lock;
pthread_mutex_t cluster_population_lock;

static int max_trials = 20;

struct find_centroid_data
{
    int num_points;
    int dim;
    int num_clusters;
    int *cluster_population;
    double *points;
    double *centroids;
    int *closest_centroid;
    int *update_flag;
};

struct determine_centroid_data
{
    int num_points;
    int dim;
    int num_clusters;
    int *cluster_population;
    double *points;
    double *centroids;
    int *closest_centroid;
};

/* Determines centroids of data points using k-means clustering algorithm.
 * Requires command line arguments for text file containing data points,
 * number of clusters, and number of threads.
 */
int main(int argc,char *argv[])
{
    double start,end;
    int i,j,k;
    int updated;
    updated = 1;

    // Get command line arguments
    char *filename = argv[1];
    int num_clusters = atoi(argv[2]);
    int num_threads = atoi(argv[3]);

    pthread_t threads[num_threads];

    int num_points;
    int dim;
    FILE *fp;
    fp = fopen(filename, "r");
    fscanf(fp, "%d", &num_points);
    fscanf(fp, "%d", &dim);

    int *thread_points;  // Array for storing beginning index of points allocated to each thread
    thread_points = (int *) malloc((num_threads + 1) * sizeof(int));


    for (i=0; i<num_threads+1; i++)
    {
        *(thread_points + i) = (int) (floor( ((double) i) / num_threads * num_points ));
    }

    double *centroids; // Stores centroid locations
    centroids = (double *) malloc(num_clusters * dim * sizeof(double));

    double *points; // Stores data point locations
    points = (double *) malloc(num_points * dim * sizeof(double));

    int *closest_centroid; // Stores index of nearest centroid for each point
    closest_centroid = (int *) malloc(num_points * sizeof(int));

    int *cluster_populations; // Stores number of points closest to each centroid
    cluster_populations = (int *) malloc(num_points * dim *sizeof(int));

    // Create mutual exclusion lock for updating ID of closest centroid to points and for updating centroids
    pthread_mutex_init(&closest_centroid_lock, NULL);
    pthread_mutex_init(&determine_centroid_lock, NULL);
    pthread_mutex_init(&cluster_population_lock, NULL);

    // Read points from file
    for (i=0; i<num_points; i++) {
        for (j=0; j<dim; j++) {
            fscanf(fp, "%lf", &points[dim * i + j]);
        }
    }

    fclose(fp);

    start = monotonic_seconds();

    // Initialize centroids and clusters
    for (i=0; i<num_clusters; i++) {
        cluster_populations[i] = 0;
        for (j=0; j<dim; j++) {
            centroids[dim * i + j] = points[dim * i + j];
        }
    }

    for (i=0; i<num_points; i++) {
        closest_centroid[i] = 0;
    }

    // Create array for passing necessary arguments to each thread for finding closest centroid to each point
    struct find_centroid_data *closest_centroid_thread_data;
    closest_centroid_thread_data = (struct find_centroid_data *) malloc(num_threads * sizeof( struct find_centroid_data ));

    // Find initial closest centroid to each data point
    for (i=0; i<num_threads; i++)
    {
        closest_centroid_thread_data[i].num_points = thread_points[i+1] - thread_points[i];
        closest_centroid_thread_data[i].dim = dim;
        closest_centroid_thread_data[i].num_clusters = num_clusters;
        closest_centroid_thread_data[i].cluster_population = cluster_populations;
        closest_centroid_thread_data[i].points = &points[ thread_points[i] * dim ];
        closest_centroid_thread_data[i].centroids = centroids;
        closest_centroid_thread_data[i].closest_centroid = &closest_centroid[ thread_points[i] ];
        closest_centroid_thread_data[i].update_flag = &updated;

        pthread_create(&threads[i], NULL, FindClosestCentroid, (void *) (closest_centroid_thread_data+i));
    }

    for (i=0; i<num_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }

    // Create array for passing necessary arguments to each thread for updating centroid locations
    struct determine_centroid_data *determine_centroid_thread_data;
    determine_centroid_thread_data = (struct determine_centroid_data *) malloc(num_threads * sizeof( struct determine_centroid_data ));

    i=0;
    updated = 1; // For tracking number of data points changing clusters in each iteration

    while( updated>0 && i<max_trials )
    {

        updated = 0;

        // Set each centroid location to origin to be updated (only use nonempty centroids)
        for (j=0; j<num_clusters; j++)
        {
            if (cluster_populations[j] > 0)
            {
                for (k=0; k<dim; k++)
                {
                    centroids[j*dim + k] = 0;
                }
            }
        }

        /*
         * Find the centroid of each cluster
         * Each thread computes contribution to averages from the subset of points allocated to it (approx. n / num_threads)
         */
        for (j=0; j<num_threads; j++)
        {
            determine_centroid_thread_data[j].num_points = thread_points[j+1] - thread_points[j];
            determine_centroid_thread_data[j].dim = dim;
            determine_centroid_thread_data[j].num_clusters = num_clusters;
            determine_centroid_thread_data[j].cluster_population = cluster_populations;
            determine_centroid_thread_data[j].points = &points[ thread_points[j] * dim];
            determine_centroid_thread_data[j].centroids = centroids;
            determine_centroid_thread_data[j].closest_centroid = &closest_centroid[ thread_points[j] ];

            pthread_create(&threads[j], NULL, DetermineCentroids, (void *) (determine_centroid_thread_data+j));
        }

        for (j=0; j<num_threads; j++)
        {
            pthread_join(threads[j], NULL);
        }

        for (j=0; j<num_clusters; j++) {
            cluster_populations[j] = 0;
        }

        /*
         * Determine closest centroid to each data point
         * Each thread computes closest centroid for each point in subset allocated to it (approx. n / num_threads)
         */
        for (j=0; j<num_threads; j++)
        {
            closest_centroid_thread_data[j].num_points = thread_points[j+1] - thread_points[j];
            closest_centroid_thread_data[j].dim = dim;
            closest_centroid_thread_data[j].num_clusters = num_clusters;
            closest_centroid_thread_data[j].cluster_population = cluster_populations;
            closest_centroid_thread_data[j].points = &points[ thread_points[j] * dim ];
            closest_centroid_thread_data[j].centroids = centroids;
            closest_centroid_thread_data[j].closest_centroid = &closest_centroid[ thread_points[j] ];
            closest_centroid_thread_data[j].update_flag = &updated;

            pthread_create(&threads[j], NULL, FindClosestCentroid, (void *) (closest_centroid_thread_data+j));
        }

        for (j=0; j<num_threads; j++)
        {
            pthread_join(threads[j], NULL);
        }

        i++;

    }

    end = monotonic_seconds();

    print_time(end-start);

    WriteOutput(closest_centroid, centroids, dim, num_clusters, num_points);

}

/*
 * Computes closest centroid to each given point.
 * Distance to every centroid is calculated. Keeps track of closest.
 * Closest centroid information is kept in local variables and updated
 * all at once after all calculations.
 */
void *FindClosestCentroid( void *input_arg )
{
    struct find_centroid_data *local_data;
    local_data = (struct find_centroid_data*) input_arg;

    int num_local_points = local_data->num_points;
    int dim = local_data->dim;
    int num_clusters = local_data->num_clusters;
    int *global_cluster_population = local_data->cluster_population;
    double *points = local_data->points;
    double *centroids = local_data->centroids;
    int *global_closest_cent = local_data->closest_centroid;  // Global variable. Need lock before updating.
    int *global_update = local_data->update_flag;

    double min_dist, dist;
    int local_closest_cent[num_local_points];
    int local_cluster_population[num_clusters];

    int closest_cent;

    int i,j,k;

    for (i=0; i<num_clusters; i++)
    {
        local_cluster_population[i] = 0;
    }

    for (i=0; i<num_local_points; i++)
    {
        dist = 0;
        for (j=0; j<dim; j++)
        {
            dist += pow( points[i*dim + j] - centroids[j], 2 );
        }
        min_dist = dist;
        closest_cent = 0;

        for (k=1; k<num_clusters; k++)
        {
            dist = 0;
            for (j=0; j<dim; j++)
            {
                dist += pow( points[i*dim + j] - centroids[k*dim + j], 2 );
            }

            if (dist < min_dist)
            {
                closest_cent = k;
                min_dist = dist;
            }
        }
        local_closest_cent[i] = closest_cent; // Update closest centroid to data point
        local_cluster_population[closest_cent] += 1; // Update number of data points associated to chosen centroid
    }

    // Lock global variables and update from local variables
    pthread_mutex_lock(&closest_centroid_lock);
    for (i=0; i<num_local_points; i++)
    {
        if (global_closest_cent[i] != local_closest_cent[i])
        {
            *global_update += 1; // Update number of data points changing clusters
            global_closest_cent[i] = local_closest_cent[i];
        }
    }
    pthread_mutex_unlock(&closest_centroid_lock);

    pthread_mutex_lock(&cluster_population_lock);
    for (i=0; i<num_clusters; i++)
    {
        global_cluster_population[i] += local_cluster_population[i];
    }
    pthread_mutex_unlock(&cluster_population_lock);

    pthread_exit(0);
}

/*
 * Finds the centroids for each cluster.
 * Set to find contribution from a subset of points by dividing by the
 * population for whole clusters, rather than population of cluster contained in
 * the given subset of points. Centroid location contributions are stored in local
 * variables and updated to global variables after calculation.
 */
void *DetermineCentroids( void *input_arg )
{
    struct determine_centroid_data *local_data;
    local_data = (struct determine_centroid_data*) input_arg;

    int num_local_points = local_data->num_points;
    int dim = local_data->dim;
    int num_clusters = local_data->num_clusters;
    int *cluster_population = local_data->cluster_population;
    double *points = local_data->points;
    double *centroids = local_data->centroids;
    int *closest_centroid = local_data->closest_centroid;

    double local_cluster_sum[num_clusters * dim]; // Local contribution to centroid coordinates

    int i,j;
    int cluster_id;

    for (i=0; i<num_clusters; i++) {
        for (j=0; j<dim; j++) {
            local_cluster_sum[i*dim + j] = 0;
        }
    }

    for (i=0; i<num_local_points; i++)
    {
        cluster_id = closest_centroid[i];

        for (j=0; j<dim; j++)
        {
            local_cluster_sum[cluster_id*dim + j] += points[i*dim + j];
        }
    }

    // Lock global variables and update from local variables
    pthread_mutex_lock(&determine_centroid_lock);

    for (i=0; i<num_clusters; i++)
    {
        for (j=0; j<dim; j++)
        {
            if (cluster_population[i] != 0)
            {
                *(centroids + i*dim + j) += local_cluster_sum[i*dim +j] / cluster_population[i];
            }
        }
    }

    pthread_mutex_unlock(&determine_centroid_lock);

    pthread_exit(0);

}

void WriteOutput( int *closest_centroid, double *centroids, int dim, int num_clusters, int num_points )
{
    int i,j;
    FILE *fp;
    fp = fopen("clusters.txt", "w");

    // Write cluster of each point to file
    for (i=0; i<num_points; i++)
    {
        fprintf(fp, "%d\n", closest_centroid[i]);
    }

    fclose(fp);

    fp = fopen("centroids.txt", "w");

    fprintf(fp, "%d %d\n", num_clusters, dim);

    // Print centroids without trailing space at end of each line
    for (i=0; i<num_clusters; i++)
    {
        for (j=0; j<dim; j++)
        {
            fprintf(fp, "%.3lf ", centroids[i*dim+j]);
        }
        //fprintf(fp, "%.3lf", centroids[i*dim+dim-1]);
        fprintf(fp, "\n");
    }

    fclose(fp);
}

/**
 * * @brief Output the seconds elapsed while clustering.
 * *
 * * @param seconds Seconds spent on k-means clustering, excluding IO.
 * */
static void print_time(double const seconds)
{
      printf("k-means clustering time: %0.04fs\n", seconds);
}
