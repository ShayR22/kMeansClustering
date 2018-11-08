#include "mpi_master_slave_header.h"
#include "slave_header.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void slaveRecvPoints(MPI_Datatype *pointDataType,point_t *points, point_t **pointsWorkLoad, int *pointsCount)
{
	MPI_Status status;

	//receive where to start from in the points array;
	int startFrom;
	MPI_Recv(&startFrom, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);
	// change pointsWorkLoad to the correct place
	*pointsWorkLoad = &(points[startFrom]);

	//recevie the number of points to work on
	MPI_Recv(pointsCount, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);
}

void slaveSendHelpers(MPI_Datatype *clusterHelperType, clusterHelper_t *helpers, int helpersCount)
{
	//SEND the helpers and after that send all the arrays one after the other
	
	int argForSend; // Arg for MPI API
	MPI_Send(&argForSend, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);

	// send the amount of points 
	vector_t *currentSumPointsLocation;
	int *currentNumOfPoints;
	int *currentPointsIDs;

	int i;
	for (i = 0; i < helpersCount; i++)
	{

		currentSumPointsLocation = &(helpers[i].sumPointLocation);
		currentNumOfPoints = &(helpers[i].numOfPoints);
		currentPointsIDs = helpers[i].pointsIDs;

		MPI_Send(currentSumPointsLocation, NUM_OF_ELEMENTS_IN_VECTOR, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD);
		MPI_Send(currentNumOfPoints, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
		if ((*currentNumOfPoints) != 0)
		{
			MPI_Send(currentPointsIDs, (*currentNumOfPoints), MPI_INT, MASTER, 0, MPI_COMM_WORLD);
		}
	}
}

//TODO FIX method recving helper override IDS allocating
void slave_recv_clusterHelpers(MPI_Datatype *clusterHelperType, clusterHelper_t *helpers, int *numOfHelpersToWorkOn)
{
	int numOfHelpersToWork;
	MPI_Status status;
	MPI_Recv(&numOfHelpersToWork, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);
	
//	printf("slave: recv numOfhelperToWorkOn = %d\n", numOfHelpersToWork);
//	fflush(stdout);

	*numOfHelpersToWorkOn = numOfHelpersToWork;

	if (numOfHelpersToWork != NO_DATA)
	{
		double *currentMaxDistance;
		vector_t *currentSumPointsLocation;
		int *currentNumOfPoints;
		int *currentPointsIDs;

		int i;
		for (i = 0; i < *numOfHelpersToWorkOn; i++)
		{
			currentMaxDistance = &(helpers[i].maxDistancePoint);
			currentSumPointsLocation = &(helpers[i].sumPointLocation);
			currentNumOfPoints = &(helpers[i].numOfPoints);
			currentPointsIDs = helpers[i].pointsIDs;
		//	printf("slave before recv: diameter: %lf, numPoints: %d\n", *currentMaxDistance, *currentNumOfPoints);
		//	fflush(stdout);

			MPI_Recv(currentMaxDistance, 1, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD, &status);
			//printf("slave: recv maxDistance = %lf\n", currentMaxDistance);
			//fflush(stdout);

			MPI_Recv(currentSumPointsLocation, NUM_OF_ELEMENTS_IN_VECTOR, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD, &status);
			//printf("slave: recv sumPointsLocation \n");
			//fflush(stdout);
			
			
			MPI_Recv(currentNumOfPoints, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);

			if ((*currentNumOfPoints) != 0)
			{
				MPI_Recv(currentPointsIDs, (*currentNumOfPoints), MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);
			}

		}
	}
	else
	{
	//	printf("slave got 0 helpers\n");
	//	fflush(stdout);
	}

}

void slave_send_clusterhelpers_diameters(clusterHelper_t *helpers, int helpersCount, double *diameterBuffer)
{

	if (helpersCount > 0)
	{
		//printf("slave send clusterHelpers diameters, helpersCount = %d\n", helpersCount);
		//fflush(stdout);
		int i;
		for (i = 0; i < helpersCount; i++)
		{
			diameterBuffer[i] = helpers[i].maxDistancePoint;
		}
		int argForSend;
		MPI_Send(&argForSend, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
		MPI_Send(diameterBuffer, helpersCount, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD);
	}
	else
	{
		//printf("slave not sending anything\n");
		//fflush(stdout);
	}

}
