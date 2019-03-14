#include "mpi_master_slave.h"
#include "slave.h"
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

void slaveSendHelpers(clusterHelper_t *helpers, int helpersCount)
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

		if ((*currentNumOfPoints) == 0)
		{
			int noPoints = -1;
			MPI_Send(&noPoints, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
		}
		else
		{
			MPI_Send(currentNumOfPoints, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
			MPI_Send(currentSumPointsLocation, NUM_OF_ELEMENTS_IN_VECTOR, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD);
			MPI_Send(currentPointsIDs, (*currentNumOfPoints), MPI_INT, MASTER, 0, MPI_COMM_WORLD);
		}
	}
}

void slave_recv_clusterHelpers(clusterHelper_t *helpers, int *numOfHelpersToWorkOn)
{
	int numOfHelpersToWork;
	MPI_Status status;
	MPI_Recv(&numOfHelpersToWork, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);
	
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
			currentMaxDistance = &(helpers[i].diameter);
			currentSumPointsLocation = &(helpers[i].sumPointLocation);
			currentNumOfPoints = &(helpers[i].numOfPoints);
			currentPointsIDs = helpers[i].pointsIDs;

			MPI_Recv(currentMaxDistance, 1, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(currentSumPointsLocation, NUM_OF_ELEMENTS_IN_VECTOR, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD, &status);	
			MPI_Recv(currentNumOfPoints, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);

			if ((*currentNumOfPoints) != 0)
			{
				MPI_Recv(currentPointsIDs, (*currentNumOfPoints), MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);
			}
		}
	}
}

void slave_send_clusterhelpers_diameters(clusterHelper_t *helpers, int helpersCount, double *diameterBuffer)
{
	if (helpersCount > 0)
	{
		int i;
		for (i = 0; i < helpersCount; i++)
		{
			diameterBuffer[i] = helpers[i].diameter;
		}
		int argForSend;
		MPI_Send(&argForSend, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
		MPI_Send(diameterBuffer, helpersCount, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD);
	}
}
