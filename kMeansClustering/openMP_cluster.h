
#ifndef OPENMP_CLUSTER_HEADER
#define OPENMP_CLUSTER_HEADER

#include "cluster.h"
#include "point.h"


/* 
Assign to each point in the points array its nearest cluster from the clusters array based on euclidean geometry

Parameters:
- clusters: pointer to clusters array
- k: length of clusters array.
- points: pointer to points array.
- pointsCount: length of points array.

Returns: void
*/
void openMP_add_to_each_point_its_nearest_cluster(cluster_t *clusters, int k, point_t *points, int pointsCount);


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
void openMP_add_to_each_clusterHelper_its_points(clusterHelper_t *helpers, int k, point_t *pointsWorkLoad, int pointCount);

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
void openMP_reset_clusterHelpers(clusterHelper_t *clusterHelpers, int clusterHelperCount);

/*
sum each helper's pointsIDs ((int*) points ids array).

Parameters:
- helpers: pointer to clusterHelpers array
- helpersCount: length of clusterHelpers array.
- points: pointer to array of points 

Returns: void
*/
void openMP_calc_clusterHelpers_points_location_sum(clusterHelper_t *helpers, int helpersCount, point_t *points);

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
int openMP_is_center_changes(cluster_t *clusters, clusterHelper_t *helpers, int k);


/*
calculate each helper's diameter

Parameters:
- helpers: pointer to clusterHelpers array
- helpersCount: length of clusterHelpers array.
- points: pointer to array of points

Returns: void
*/
void openMP_calc_clusterHelpers_diamaters(clusterHelper_t *helpers, int helpersCount, point_t *allPoints);


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
double openMP_calc_QM(cluster_t *clusters, int clustersCount);

#endif // !OPENMP_CLUSTER_HEADER

