#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void *FindClosestCentroid( void *input_arg );
void *DetermineCentroids( void *input_arg );
void WriteOutput( int *closest_centroid, double *centroids, int dim, int num_clusters, int num_points );
pthread_mutex_t closest_centroid_lock;
pthread_mutex_t determine_centroid_lock;
pthread_mutex_t cluster_population_lock;

struct find_centroid_data
{
    int num_points;
    int dim;
    int num_clusters;
    int min_global_point_id;  // minimum point ID in global numbering
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

int main(int argc,char *argv[])
{
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

    int points_per_thread = num_points / num_threads;
    int *thread_points;
    thread_points = (int *) malloc((num_threads + 1) * sizeof(int));


    for (i=0; i<num_threads+1; i++)
    {
        *(thread_points + i) = (int) (floor( ((double) i) / num_threads * num_points ));
    }

    double *centroids;
    centroids = (double *) malloc(num_clusters * dim * sizeof(double));

    double *points;
    points = (double *) malloc(num_points * dim * sizeof(double));

    int *closest_centroid;
    closest_centroid = (int *) malloc(num_points * sizeof(int));

    int *cluster_populations;
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

    struct find_centroid_data *closest_centroid_thread_data;
    closest_centroid_thread_data = (struct find_centroid_data *) malloc(num_threads * sizeof( struct find_centroid_data ));

    for (i=0; i<num_threads; i++)
    {
        closest_centroid_thread_data[i].num_points = thread_points[i+1] - thread_points[i];
        closest_centroid_thread_data[i].dim = dim;
        closest_centroid_thread_data[i].num_clusters = num_clusters;
        closest_centroid_thread_data[i].min_global_point_id = thread_points[i];
        closest_centroid_thread_data[i].cluster_population = cluster_populations;
        closest_centroid_thread_data[i].points = &points[ thread_points[i] * dim ];
        closest_centroid_thread_data[i].centroids = centroids;
        closest_centroid_thread_data[i].closest_centroid = closest_centroid;
        closest_centroid_thread_data[i].update_flag = &updated;

        pthread_create(&threads[i], NULL, FindClosestCentroid, (void *) (closest_centroid_thread_data+i));
    }

    for (i=0; i<num_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }

    struct determine_centroid_data *determine_centroid_thread_data;
    determine_centroid_thread_data = (struct determine_centroid_data *) malloc(num_threads * sizeof( struct determine_centroid_data ));

    i=0;
    updated = 1;

    while( updated>0 && i<20 )
    {

        updated = 0;

        for (j=0; j<num_threads; j++)
        {
            for (k=0; k<dim; k++)
            {
                centroids[j*dim + k] = 0;
            }
        }

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

        for (j=0; j<num_threads; j++)
        {
            closest_centroid_thread_data[j].num_points = thread_points[j+1] - thread_points[j];
            closest_centroid_thread_data[j].dim = dim;
            closest_centroid_thread_data[j].num_clusters = num_clusters;
            closest_centroid_thread_data[j].min_global_point_id = thread_points[j];
            closest_centroid_thread_data[j].cluster_population = cluster_populations;
            closest_centroid_thread_data[j].points = &points[ thread_points[j] * dim ];
            closest_centroid_thread_data[j].centroids = centroids;
            closest_centroid_thread_data[j].closest_centroid = closest_centroid;
            closest_centroid_thread_data[j].update_flag = &updated;

            pthread_create(&threads[j], NULL, FindClosestCentroid, (void *) (closest_centroid_thread_data+j));
        }

        for (j=0; j<num_threads; j++)
        {
            pthread_join(threads[j], NULL);
        }

        i++;

    }

    WriteOutput(closest_centroid, centroids, dim, num_clusters, num_points);

}

void *FindClosestCentroid( void *input_arg )
{
    struct find_centroid_data *local_data;
    local_data = (struct find_centroid_data*) input_arg;

    int num_local_points = local_data->num_points;
    int dim = local_data->dim;
    int num_clusters = local_data->num_clusters;
    int min_point_id = local_data->min_global_point_id;
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
            }
        }
        local_closest_cent[i] = closest_cent;
        local_cluster_population[closest_cent] += 1;
    }

    pthread_mutex_lock(&closest_centroid_lock);
    for (i=0; i<num_local_points; i++)
    {
        if (global_closest_cent[min_point_id + i] != local_closest_cent[i])
        {
            *global_update += 1;
            global_closest_cent[min_point_id + i] = local_closest_cent[i];
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

    double local_cluster_sum[num_clusters * dim];

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
            if (cluster_population[cluster_id] > 0) {
                local_cluster_sum[cluster_id*dim + j] += points[i*dim + j] / cluster_population[cluster_id];
            }
        }
    }

    pthread_mutex_lock(&determine_centroid_lock);

    for (i=0; i<num_clusters; i++)
    {
        for (j=0; j<dim; j++)
        {
            *(centroids + i*dim + j) += local_cluster_sum[i*dim +j];
        }
    }

    pthread_mutex_unlock(&determine_centroid_lock);

    pthread_exit(0);

}

void WriteOutput( int *closest_centroid, double *centroids, int dim, int num_clusters, int num_points )
{
    int i,j;
    FILE *fp1, *fp2;
    fp1 = fopen("clusters.txt", "w");

    for (i=0; i<num_points; i++)
    {
        fprintf(fp1, "%d\n", closest_centroid[i]);
    }

    fclose(fp1);

    fp2 = fopen("centroids.txt", "w");

    for (i=0; i<num_clusters; i++)
    {
        for (j=0; j<dim; j++)
        {
            fprintf(fp2, "%.3lf ", centroids[i*dim+j]);
        }
        fprintf(fp2, "\n");
    }

    fclose(fp2);
}
