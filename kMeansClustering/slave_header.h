
#ifndef SLAVE_HEADER
#define SLAVE_HEADER
#include "mpi_master_slave_header.h"
#include "cluster_header.h"
#include <mpi.h>
#include "point_header.h"

void slaveRecvPoints(MPI_Datatype *pointDataType,point_t *points,  point_t **pointsWorkLoad, int *pointsCount);

void slaveSendHelpers(MPI_Datatype *clusterHelperType, clusterHelper_t *helpers, int helpersCount);

void slave_recv_clusterHelpers(MPI_Datatype *clusterHelperType, clusterHelper_t *helpers, int *numOfHelpersToWorkOn);

void slave_send_clusterhelpers_diameters(clusterHelper_t *helpers, int helpersCount, double *diameterBuffer);

#endif // !SLAVE_HEADER
