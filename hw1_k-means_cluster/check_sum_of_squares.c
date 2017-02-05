#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char *argv[])
{
    int i,j;
    int cluster_id;

    int num_points;
    int dim;
    int num_clusters;
    double *points;
    double *centroids;
    int *closest_centroid;

    FILE *fp;
    char *filename = argv[1];

    centroids = (double *) malloc(num_clusters * dim * sizeof(double));
    points = (double *) malloc(num_points * dim * sizeof(double));
    closest_centroid = (int *) malloc(num_points * sizeof(int));

    // Read data points from file
    fp = fopen(filename, "r");
    fscanf(fp, "%d", &num_points);
    fscanf(fp, "%d", &dim);

    for (i=0; i<num_points; i++) {
        for (j=0; j<dim; j++) {
            fscanf(fp, "%lf", &points[dim * i +j]);
        }
    }

    // Read centroids from file
    filename = argv[2];
    fp = fopen(filename, "r");
    fscanf(fp, "%d", &num_clusters);
    fscanf(fp, "%d", &dim);

    for (i=0; i<num_clusters; i++) {
        for (j=0; j<dim; j++) {
            fscanf(fp, "%lf", &centroids[dim*i+j]);
        }
    }

    // Read closest centroid from file
    filename = argv[3];
    fp = fopen(filename, "r");

    for (i=0; i<num_points; i++) {
        fscanf(fp, "%d", &closest_centroid[i]);
    }

    double sum_of_dist = 0;

    for (i=0; i<num_points; i++)
    {
        cluster_id = closest_centroid[i];
        for (j=0; j<dim; j++)
        {
            sum_of_dist += pow(points[i*dim+j] - centroids[cluster_id*dim+j],2);
        }
    }

    printf("The sum of squared distances is %lf\n", sum_of_dist);

}
