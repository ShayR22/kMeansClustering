
#ifndef CUDA_HEADER
#define CUDA_HEADER
#include "cluster.h"
#include "point.h"

/*
	Increments an array of points by a deltT.
	this increment is based on the formula: newLocation = oldLocation + speed * deltaT.
	
	Parameters:
	- points: pointer to array of point_t.
    - pointsCouns:length of points array.
	- dletaT: amount of time to advance the points's location.

	Returns:
	 - 1: for success.
	 - 0: for failure.
*/
int cuda_increment_points(point_t *points, int pointsCount, double deltaT);

/*
classify for each points the cluster to which its belongs
this classification is oclidian based

Parameters:
- points: pointer to points array
- pointsCount:length of points array.
- clusters: pointer to cluster array.
- clustersCount: length of the clusters array

Returns:
- 1: for success.
- 0: for failure.
*/
int cuda_add_to_each_point_its_nearest_cluster(point_t *points, int pointsCount, cluster_t *clusters, int clustersCount);

/*
classify for each points the cluster to which its belongs
this classification is oclidian based

Parameters:
- points: pointer to points array
- pointsCount:length of points array.
- clusters: pointer to clusterHelpers array.
- clustersCount: length of the clusterHelpers array

Returns:
- 1: for success.
- 0: for failure.
*/
int cuda_calc_each_helper_diamater(point_t *points, int pointsCount, clusterHelper_t *helpers, int helpersCount);

#endif //!CUDA_HEADER