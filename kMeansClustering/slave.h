
#ifndef SLAVE_HEADER
#define SLAVE_HEADER
#include "mpi_master_slave.h"
#include "cluster.h"
#include "point.h"
#include <mpi.h>

/*
Assumption: master and slave has the same points array
Recevie where to work on inside the points array and how many points to work on

Parameters:
- pointDataType: mpi structure representing a point_t sturct.
- points: pointer to points array.
- pointsCount: length of points array.
- pointsWorkLoad: a pointer to be initialize the future work starting point.
- pointsCount: a pointer to be initialize the number of points to work on.

Returns: void
*/
void slaveRecvPoints(MPI_Datatype *pointDataType, point_t *points,  point_t **pointsWorkLoad, int *pointsCount);

/*
Send an array of helpers to master (master is initialize to id 0).

Parameters:
- helpers: pointer to clusterHelpers array.
- helpersCount: length of helpers array.

Returns: void
*/
void slaveSendHelpers(clusterHelper_t *helpers, int helpersCount);

/*
Assumption number of helpers that will be recvied wont be larger than helpers size
Recv amount of helpers to be recevied and then recv those amount of helpers.

Parameters:
- numOfHelpersToWorkOn: pointer to be initialize with the number of helpers to work on and recvied in the method.
- helpers: pointer to clusterHelpers array.

Returns: void.
*/
void slave_recv_clusterHelpers(clusterHelper_t *helpers, int *numOfHelpersToWorkOn);


/*
Fill diameter buffer with the helpers diameters and send it to master (master is initialize to id 0).

Parameters:
- helpers: pointer to clusterHelpers array.
- helpersCount: length of helpers array.
- diameterBuffer: pointer that would conatin all the helpers diameters

Returns: void.
*/
void slave_send_clusterhelpers_diameters(clusterHelper_t *helpers, int helpersCount, double *diameterBuffer);

#endif // !SLAVE_HEADER
