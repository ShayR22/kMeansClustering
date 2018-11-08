
#ifndef OPENMP_CLUSTER_HEADER
#define OPENMP_CLUSTER_HEADER

#include "cluster_header.h"
#include "point_header.h"

void openMP_add_to_each_point_its_nearest_cluster(cluster_t *clusters, int k, point_t *points, int pointsCount);
void openMP_add_to_each_clusterHelper_its_points(clusterHelper_t *helpers, int k, point_t *pointsWorkLoad, int pointCount);

void openMP_reset_clusterHelpers(clusterHelper_t *clusterHelpers, int clusterHelperCount);
void openMP_calc_clusterHelpers_points_location_sum(clusterHelper_t *helpers, int helpersCount, point_t *points);

int openMP_is_center_changes(cluster_t *clusters, clusterHelper_t *helpers, int k);

void openMP_calc_clusterHelpers_diamaters(clusterHelper_t *helpers, int helpersCount, point_t *allPoints);

double openMP_calc_QM(cluster_t *clusters, int clustersCount);

#endif // !OPENMP_CLUSTER_HEADER

