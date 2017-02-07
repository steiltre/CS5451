#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "timing.h"

void FindClosestCentroid( int num_points, int dim, int num_clusters, int *cluster_population, double *points, double *centroids, int *closest_centroid, int *update_flag );
void DetermineCentroids( int num_points, int dim, int num_clusters, int *cluster_population, double *points, double *centroids, int *closest_centroid );
void WriteOutput( int *closest_centroid, double *centroids, int dim, int num_clusters, int num_points );
static void print_time( double const seconds );

static int max_trials = 20;

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
    int err;

    // Get command line arguments
    char *filename = argv[1];
    int num_clusters = atoi(argv[2]);
    int num_threads = atoi(argv[3]);

    omp_set_num_threads(num_threads);

    int num_points;
    int dim;
    FILE *fp;
    fp = fopen(filename, "r");
    err = fscanf(fp, "%d", &num_points);
    err = fscanf(fp, "%d", &dim);

    double *centroids; // Stores centroid locations
    centroids = (double *) malloc(num_clusters * dim * sizeof(double));

    double *points; // Stores data point locations
    points = (double *) malloc(num_points * dim * sizeof(double));

    int *closest_centroid; // Stores index of nearest centroid for each point
    closest_centroid = (int *) malloc(num_points * sizeof(int));

    int *cluster_populations; // Stores number of points closest to each centroid
    cluster_populations = (int *) malloc(num_points * dim *sizeof(int));

    // Read points from file
    for (i=0; i<num_points; i++) {
        for (j=0; j<dim; j++) {
            err = fscanf(fp, "%lf", &points[dim * i + j]);
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

    FindClosestCentroid( num_points, dim, num_clusters, cluster_populations, points, centroids, closest_centroid, &updated );

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

        DetermineCentroids(num_points, dim, num_clusters, cluster_populations, points, centroids, closest_centroid);

        for (j=0; j<num_clusters; j++) {
            cluster_populations[j] = 0;
        }

        FindClosestCentroid( num_points, dim, num_clusters, cluster_populations, points, centroids, closest_centroid, &updated );

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
void FindClosestCentroid( int num_points, int dim, int num_clusters, int *cluster_population, double *points, double *centroids, int *closest_centroid, int *update_flag)
{

    #pragma omp parallel
    {
        // Initialize local variables
        int local_updates = 0;
        int local_closest_cent[num_points];
        int local_cluster_population[num_clusters];
        int i;
        for (i=0; i<num_clusters; i++)
        {
            local_cluster_population[i] = 0;
        }

        for (i=0; i<num_points; i++)
        {
            local_closest_cent[i] = -1; // Set to unused value to easily tell if it's been updated
        }

        #pragma omp for schedule(static)
        for (i=0; i<num_points; i++)
        {
            double dist, min_dist;
            int j,k;
            int closest_cent;
            dist = 0;
            for (j=0; j<dim; j++)
            {
                dist += pow( points[i*dim+j] - centroids[j], 2);
            }
            min_dist = dist;
            closest_cent = 0;

            for (k=1; k<num_clusters; k++)
            {
                dist = 0;
                for (j=0; j<dim; j++)
                {
                    dist += pow( points[i*dim+j] - centroids[k*dim+j], 2);
                }

                if (dist < min_dist)
                {
                    closest_cent = k;
                    min_dist = dist;
                }
            }

            local_closest_cent[i] = closest_cent;
            local_cluster_population[closest_cent] += 1;

           }

        // Update global variables with local values
        #pragma omp critical (closest_centroid)
        {
            for(i=0; i<num_points; i++)
            {
                if ((local_closest_cent[i] != -1) & (local_closest_cent[i] != closest_centroid[i]))
                {
                    local_updates += 1; // Track number of points that changed cluster
                    closest_centroid[i] = local_closest_cent[i];
                }
            }
        }

        #pragma omp critical (cluster_population)
        {
            for (i=0; i<num_points; i++)
            {
                cluster_population[i] += local_cluster_population[i];
            }
        }

        #pragma omp critical (update_flag)
        {
            *update_flag += local_updates;
        }
    }
}

/*
 * Finds the centroids for each cluster.
 * Set to find contribution from a subset of points by dividing by the
 * population for whole clusters, rather than population of cluster contained in
 * the given subset of points. Centroid location contributions are stored in local
 * variables and updated to global variables after calculation.
 */
void DetermineCentroids( int num_points, int dim, int num_clusters, int *cluster_population, double *points, double *centroids, int *closest_centroid )
{

    int i,j;
    double local_cluster_sum[num_clusters * dim]; // Local contribution to centroid coordinates
    double cluster_weights[num_clusters]; // Weight to use for averaging in each cluster

    for (i=0; i<num_clusters; i++) {
        for (j=0; j<dim; j++) {
            local_cluster_sum[i*dim + j] = 0;
        }

        if (cluster_population[i] > 0)
        {
            cluster_weights[i] = ( (double) 1 ) / cluster_population[i];
        }
        else
        {
            cluster_weights[i] = 0;
        }
    }

    #pragma omp parallel firstprivate(local_cluster_sum)
    {
        int j;
        #pragma omp for schedule(static)
        for (i=0; i<num_points; i++)
        {
            int cluster_id;
            cluster_id = closest_centroid[i];

            // Add points weighted coordinates to cluster centroid
            for (j=0; j<dim; j++)
            {
                local_cluster_sum[cluster_id*dim+j] += points[i*dim+j] * cluster_weights[cluster_id];
            }
        }

        // Update global values with local contributions
        #pragma omp critical (determine_centroid)
        {
            for (i=0; i<num_clusters; i++)
            {
                for (j=0; j<dim; j++)
                {
                    centroids[i*dim+j] += local_cluster_sum[i*dim+j];
                }
            }
        }
    }
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
