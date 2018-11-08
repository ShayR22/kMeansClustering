
#ifndef CLUSTER_HEADER
#define CLUSTER_HEADER

#include "point_header.h"
#include "vector_header.h"

#define NUMBER_OF_ELEMENTS_IN_CLUSTER 3
#define NUMBER_OF_ELEMENTS_IN_CLUSTER_HELPER 3

typedef struct CLUSTER
{
	//TODO currently id is not needed check 
	int id;
	vector_t center;
	double maxDistancePoint;

}cluster_t;


typedef struct CLUSTER_HELPER
{
	double maxDistancePoint;
	vector_t sumPointLocation;
	int numOfPoints;
	int *pointsIDs;
}clusterHelper_t;


void cluster_init(cluster_t *cluster, vector_t *center);
void clusterHelper_allocate(clusterHelper_t *helpers, int helpersCount ,int numPoints);

void cluster_add_to_each_point_its_nearest_cluster(cluster_t *clusters, int k, point_t *points, int pointsCount);
void cluster_add_to_each_clusterHelper_its_points(clusterHelper_t *helpers, int k, point_t *pointsWorkLoad, int pointCount);

void cluster_calc_clusterHelper_points_location_sum(cluster_t *clusters, clusterHelper_t *clusterHelpers, int k, point_t *points);
void cluster_reset_clusterHelpers(clusterHelper_t *helpers, int helpersCount);

int cluster_is_center_changes(cluster_t *clusters, clusterHelper_t *helpers, int k);

void cluster_copy_clusterHelpers(clusterHelper_t *copyFrom, clusterHelper_t *copyTo, int numOfHelpers);

void cluster_calculate_helpers_diameter(clusterHelper_t *helpers, int helpersCount, point_t *points, int pointsCount);

double cluster_calculate_quality_measure(cluster_t *clusters, int clustersCount);


MPI_Datatype cluster_create_cluster_mpi_struct();
MPI_Datatype cluster_create_clusterHelper_mpi_struct();


void cluster_print_clusters(cluster_t *clusters, int clustersCount);
void cluster_print_clusterHelpers(clusterHelper_t *helpers, int helpersCount);


#endif // !CLUSTER_HEADER
