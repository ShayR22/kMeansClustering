
#ifndef MASTER_HEADER
#define MASTER_HEADER
#include "mpi_master_slave_header.h"
#include "cluster_header.h"
#include "point_header.h"
#include <mpi.h>

// this method returns the number of points the master need to work on
void master_distribute_points(MPI_Datatype *pointDataType,point_t *points, int pointsCount, int numProcs,
						   point_t **pointsWorkLoad, int *numPoints);

void master_recv_clusterHelpers(clusterHelper_t *helpers, int helpersCount, int numProcs, int numOfPoints, MPI_Datatype *clusterHelperDataType);
void master_fold_clusterHelpers(clusterHelper_t *helpers, int helperCount, int numProcs);

void master_send_clusterHelpers(clusterHelper_t *helpers, int helpersCount, int numProcs, int *helpersToWorkOn, MPI_Datatype *clusterHelperDataType);

void master_copy_diameters_from_its_own_helpers(cluster_t *clusters, clusterHelper_t *helpers, int helpersCount);


//recv helpers results directly into the clusters
void master_recv_clusters_diameters(cluster_t *clusters, int helpersCoumt, int numProcs, int numOfHelpersToWorkOn, double *diameterBuffer);
#endif // !MASTER_HEADER
