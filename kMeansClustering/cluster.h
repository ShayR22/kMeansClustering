
#ifndef CLUSTER_HEADER
#define CLUSTER_HEADER

#include "point.h"
#include "vector.h"

#define NUMBER_OF_ELEMENTS_IN_CLUSTER 2
#define NUMBER_OF_ELEMENTS_IN_CLUSTER_HELPER 3

typedef struct CLUSTER
{
	vector_t center;
	double diameter;

}cluster_t;


typedef struct CLUSTER_HELPER
{
	double diameter;
	vector_t sumPointLocation;
	int numOfPoints;
	int *pointsIDs;
}clusterHelper_t;

/*
copy center into the cluster and assign the cluster an id 

Parameters:
- cluster: pointer to a cluster
- center: pointer to a center vecotr_t

Returns: void

*/
void cluster_init(cluster_t *cluster, vector_t *center);

/*
allocate for each clusterHelper in the helpers array its array of pointsIDs

Parameters:
- helpers: pointer to a clusteHelpers arrays
- helpersCount: length of helpers array
- numPoints: amount of ids to allocate for each clusterHelper

Returns: void

*/
void clusterHelper_allocate(clusterHelper_t *helpers, int helpersCount ,int numPoints);

/*
Assign to each point in the points array its nearest cluster from the clusters array based on euclidean geometry

Parameters:
- clusters: pointer to clusters array
- k: length of clusters array.
- points: pointer to points array.
- pointsCount: length of points array.

Returns: void
*/
void cluster_add_to_each_point_its_nearest_cluster(cluster_t *clusters, int k, point_t *points, int pointsCount);

/*
Assign to each pointsIDs ((int*) points id array) in each clusterHelper (helpers array), its
points ids from pointsWorkLoad (points array).

Parameters:
- helpers: pointer to clusterHelpers array
- k: length of clusterHelpers array.
- pointsWorkLoad: pointer to points array.
- pointsCount: length of points array.

Returns: void
*/
void cluster_add_to_each_clusterHelper_its_points(clusterHelper_t *helpers, int k, point_t *pointsWorkLoad, int pointCount);

/*
sum each helper's pointsIDs ((int*) points ids array).

Parameters:
- helpers: pointer to clusterHelpers array
- helpersCount: length of clusterHelpers array.
- points: pointer to array of points

Returns: void
*/
void cluster_calc_clusterHelper_points_location_sum(clusterHelper_t *helpers, int k, point_t *points);

/*
reset each clusterHelper in the array by zeroing its:
- diameter
- sumPointLocation
- numOfPoints

Parameters:
- clusterHelpers: pointer to clusterHelpers array
- clusterHelperCount: length of clusterHelpers array.

Returns: void
*/
void cluster_reset_clusterHelpers(clusterHelper_t *helpers, int helpersCount);

/*
compare each cluster center to its corrosponding (array index) clusterHelper center and update the cluster center
if they are different

Parameters:
- clusters: pointer to clusters array
- helpers: pointer to clusterHelpers array
- k: length of clusters array and clusterHelpers array.

Returns:
- 1 there was 1 pair that didnt have the same center
- 0 all pair centers are the same
*/
int cluster_is_center_changes(cluster_t *clusters, clusterHelper_t *helpers, int k);

/*
copy from each clusterHelper in "copyFrom" to its corresponding clusterHelper in "copyTo" its:
- diameter
- sumPointLocation
- numOfPoints

Parameters:
- copyFrom: pointer to clusterHelpers array
- copyTo: pointer to clusterHelpers array
- numOfHelpers: length of helpers to copy

Returns: void
*/
void cluster_copy_clusterHelpers(clusterHelper_t *copyFrom, clusterHelper_t *copyTo, int numOfHelpers);

/*
calculate each helper's diameter

Parameters:
- helpers: pointer to clusterHelpers array
- helpersCount: length of clusterHelpers array.
- points: pointer to array of points
- pointsCount: length of points array

Returns: void
*/
void cluster_calculate_helpers_diameter(clusterHelper_t *helpers, int helpersCount, point_t *points, int pointsCount);

/*
calculate the quality measure of an array of clusters.
calculation is done by the formula:
q = (di+ dj/Dij)/ (clustersCount * clusterCount -1).
di represents the diameter of cluster number i.
Dij represents the distance between the 2 centers of cluster number i and cluster number j.

Parameters:
- clusters: pointer to clusters array
- clustersCount: length of clusters array and clusterHelpers array.

Returns:
- double representing the clusters quality meassure
*/
double cluster_calculate_quality_measure(cluster_t *clusters, int clustersCount);

/*
create and commit cluster mpi sturct for later use

Parameters: void

Returns:
- MPI_DataType
*/
MPI_Datatype cluster_create_cluster_mpi_struct();

/*
free array of clusterHelpers

Parameters:
- helpers: pointer to clusterHelpers array
- helpersCount: length of clusterHelpers array.

Returns: void
*/
void cluster_free_clusterHelpers_array(clusterHelper_t *helpers, int helpersCount);

/*
print to console array of clusters

Parameters:
- clusters: pointer to clusters array
- clustersCount: length of clusters array.

Returns: void
*/
void cluster_print_clusters(cluster_t *clusters, int clustersCount);

#endif // !CLUSTER_HEADER
