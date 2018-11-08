
#ifndef CUDA_HEADER
#include "cluster_header.h"
#include "point_header.h"



#define CUDA_HEADER

//TODO add calc QM method

int cuda_increment_points(point_t *points, int pointsCount, double deltaT);
int cuda_add_to_each_point_its_nearest_cluster(point_t *points, int pointsCount, cluster_t *clusters, int clustersCount);
int cuda_calc_each_helper_diamater(point_t *points, int pointsCount, clusterHelper_t *helpers, int helpersCount);

//TODO this method isn't viable because copying such large data (points array too) takes too long
//int cuda_add_to_each_clusterHelper_its_points(clusterHelper_t *helpers, int helpersCount, point_t *points, int pointsCount);

#endif