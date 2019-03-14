#ifndef MASTER_HEADER
#define MASTER_HEADER
#include "mpi_master_slave.h"
#include "cluster.h"
#include "point.h"
#include <mpi.h>

/*
This method require all machine/process have the same points array
Distribute the points across the proccess that are available (including to the master itself).
Distribution is done by sending the 2 ints:
 - first: where to start working inside the points array
 - second: the amount of points to work on, within the array

Parameters:
- pointDataType: mpi structure representing a point_t sturct.
- points: pointer to points array.
- pointsCount: length of points array.
- numProcs: number of machines/processes.
- pointsWorkLoad: a pointer to be initialize to the future work starting point.
- numPoins: a pointer to be initialize the number of points to work on.

Returns: void
*/
void master_distribute_points(MPI_Datatype *pointDataType, point_t *points, int pointsCount, int numProcs,
						      point_t **pointsWorkLoad, int *numPoints);

/*
Assumption: helpers length is at least helpersCount*numProcs
Recevie helpers from other machines/processes

Parameters:
- helpers: pointer to clusterHelpers array.
- helperCount: length of helpers array.
- numProcs: number of machines/processes.

Returns: void
*/
void master_recv_clusterHelpers(clusterHelper_t *helpers, int helpersCount, int numProcs);

/*
Assumption: helpers length is at least helpersCount*numProcs
Fold the results across the sets of helpers (helpersCount size per set) into the first set.

Parameters:
- helpers: pointer to clusterHelpers array.
- helperCount: length of helpers array.
- numProcs: number of machines/processes.

Returns: void
*/
void master_fold_clusterHelpers(clusterHelper_t *helpers, int helperCount, int numProcs);

/*
Scatter helpers across machines/processes by calculating how many helpers per machine/process and send them to it

Parameters:
- helpers: pointer to clusterHelpers array.
- helperCount: length of helpers array.
- numProcs: number of machines/processes.
- helpersToWorkOn: pointer to be initialize to the future work starting index.

Returns: void
*/
void master_send_clusterHelpers(clusterHelper_t *helpers, int helpersCount, int numProcs, int *helpersToWorkOn);

/*
Scatter helpers across machines/processes

Parameters:
- helpers: pointer to clusterHelpers array.
- helperCount: length of helpers array.
- numProcs: number of machines/processes.
- helpersToWorkOn: pointer to be initialize to the future work starting index.

Returns: void
*/
void master_copy_diameters_from_its_own_helpers(cluster_t *clusters, clusterHelper_t *helpers, int helpersCount);


/*
Scatter helpers across machines/processes

Parameters:
- clusters: pointer to clusters array.
- clustersCount: length of clusters array.
- numProcs: number of machines/processes.
- numOfDiametersToRecvFromSlave: number of diameters to recevie from each slave. (except last slave might
 have more according to clustersCount/numProcs (remainder)).

Returns: void
*/
void master_recv_clusters_diameters(cluster_t *clusters, int clustersCount, int numProcs, int numOfDiametersToRecvFromSlave, double *diameterBuffer);

#endif // !MASTER_HEADER
