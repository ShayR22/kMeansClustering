#include "mpi_master_slave_header.h"
#include "master_header.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

static void recv_clusterHelpers(clusterHelper_t *helpers, int helpersCount, MPI_Datatype *clusterHelperDataType);
static void send_helpers(int sendTo, int numOfHelpersToWorkOn, clusterHelper_t *helpers, MPI_Datatype *clusterHelperDataType);
static void recv_diameters(cluster_t *clusters, int numProcs, int helpersCount, int recvAmount, double *diameterBuffer);

void master_distribute_points(MPI_Datatype *pointDataType, point_t *points, int pointsCount, int numProcs,
	point_t **pointsWorkLoad, int *numPoints)
{
	int numOfPointsToWorkOn = pointsCount / numProcs;
	*numPoints = numOfPointsToWorkOn;
	*pointsWorkLoad = points; // master will work from the start

	int i;
	for (i = 1; i < numProcs; i++) //i start from 1 to exclude master
	{
		int startFrom = i * numOfPointsToWorkOn;
		if (i == numProcs - 1)
		{
			//add remainder of points to last process/machine
			numOfPointsToWorkOn += pointsCount - numProcs*numOfPointsToWorkOn;
		}
		// send where to start and how many points to work on
		MPI_Send(&startFrom, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Send(&numOfPointsToWorkOn, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	}
}


void master_recv_clusterHelpers(clusterHelper_t *helpers, int helpersCount, int numProcs, int numOfPoints, MPI_Datatype *clusterHelperDataType)
{
	int i;
	for (i = 1; i < numProcs; i++)// starts from 1 to not include master
	{
		recv_clusterHelpers( helpers, helpersCount, clusterHelperDataType);
	}

	//TEST START
	//printf("master: printing all helpers after recv\n");
	int j;
	for (i = 0; i < helpersCount * numProcs; i++)
	{
	//	printf("helper number %d has %d points\n", i, helpers[i].numOfPoints);
	//	fflush(stdout);
		for (j = 0; j < helpers[i].numOfPoints; j++)
		{
			//printf("%d, ", helpers[i].pointsIDs[j]);
			//fflush(stdout);
		}
	//	printf("\n");
	//	fflush(stdout);
	}
	//TEST END
}

static void recv_clusterHelpers(clusterHelper_t *helpers, int helpersCount, MPI_Datatype *clusterHelperDataType)
{
	int sentFrom;
	MPI_Status status;

	int argForRecv; // argument only for MPI API
	MPI_Recv(&argForRecv, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
	sentFrom = status.MPI_SOURCE;
	clusterHelper_t *locationToAdd = &(helpers[sentFrom * helpersCount]);
	
	vector_t *currentSumPointsLocation;
	int *currentNumOfPoints;
	int *currentPointsIDs;

	int i;
	for (i = 0; i < helpersCount; i++)
	{
		currentSumPointsLocation = &(locationToAdd[i].sumPointLocation);
		currentNumOfPoints = &(locationToAdd[i].numOfPoints);
		currentPointsIDs = locationToAdd[i].pointsIDs;

		MPI_Recv(currentSumPointsLocation, NUM_OF_ELEMENTS_IN_VECTOR, MPI_DOUBLE, sentFrom, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(currentNumOfPoints, 1, MPI_INT, sentFrom, 0, MPI_COMM_WORLD, &status);
		if (*currentNumOfPoints != 0)
		{
			MPI_Recv(currentPointsIDs, *currentNumOfPoints, MPI_DOUBLE, sentFrom, 0, MPI_COMM_WORLD, &status);
		}
	}

//	printf("\nmaster recv helpers from %d\n", sentFrom);
//	fflush(stdout);
	for (i = 0; i < helpersCount; i++)
	{
	//	printf("\nbuffer in %d:\n", i);
	//	fflush(stdout);
	//	printf("numPoints = %d, diameter = %lf sumLocation is: \n", locationToAdd[i].numOfPoints, locationToAdd[i].maxDistancePoint);
	//	fflush(stdout);
	//	vector_print_vector(&(locationToAdd[i].sumPointLocation));
	//	printf("\n");
	//	fflush(stdout);
	}

}

void master_fold_clusterHelpers(clusterHelper_t *helpers, int helperCount, int numProcs)
{
	clusterHelper_t *masterHelper;
	clusterHelper_t *currentFold;
	int masterCurrentNumPoints;
	int *pointsLastLocation;
	int numberOfBytesToCopy;


	// go over indexes and sum them rather summing each set
	/*
		we have 2 procs and 2 clusters meaning master has 0,1,2,3
		we start saying masterHelper = 0;
		current fold = 1*2 + 0 = 2
	*/

	int i, j;
	for (i = 0; i < helperCount; i++)
	{
		masterHelper = &(helpers[i]);
		for (j = 1; j < numProcs; j++)
		{
			currentFold = &(helpers[j*helperCount + i]);

			//add all the currentFold pointsID to the master's pointsID
			masterCurrentNumPoints = masterHelper->numOfPoints;
			pointsLastLocation = &(masterHelper->pointsIDs[masterCurrentNumPoints]);
			numberOfBytesToCopy = (currentFold->numOfPoints) * sizeof(int);

			memcpy(pointsLastLocation, currentFold->pointsIDs, numberOfBytesToCopy);

			// update the new number of points
			masterHelper->numOfPoints += currentFold->numOfPoints;

			//add the sum of the vectors
			vector_add_vector(&(masterHelper->sumPointLocation), &(currentFold->sumPointLocation));

		}
		//devide the sum of the vector in the number of points
		vector_divide_vector(&(masterHelper->sumPointLocation), masterHelper->numOfPoints);
	}
}

void master_send_clusterHelpers(clusterHelper_t *helpers, int helpersCount, int numProcs, int *helpersToWorkOn, MPI_Datatype *clusterHelperDataType)
{
	int numOfHelpersToWorkOn = helpersCount / numProcs;

	*helpersToWorkOn = numOfHelpersToWorkOn;

	clusterHelper_t *helpersToSend = NULL;
	int i;
	for (i = 1; i < numProcs; i++) //i starts from 1 to exclude master
	{
		helpersToSend = &(helpers[i * numOfHelpersToWorkOn]);
		if (i == numProcs - 1)
		{
			numOfHelpersToWorkOn += helpersCount - (numOfHelpersToWorkOn * numProcs); // last process calc remainder
		}
		send_helpers(i, numOfHelpersToWorkOn, helpersToSend, clusterHelperDataType);
	}
}

static void send_helpers(int sendTo, int numOfHelpersToWorkOn, clusterHelper_t *helpers, MPI_Datatype *clusterHelperDataType)
{
	/*
	if number of helpers to send is 0, a noData(-1) will be sent instead.
	in a for loop helpers will be sent:
		maxDistance
		sumPointLocation
		numOfPoints
		pointsIDS
	special case numOfPoints equal 0 no pointsIDS will not be sent
	*/


	int noData = -1;
	if (numOfHelpersToWorkOn == 0)
	{
		MPI_Send(&noData, 1, MPI_INT, sendTo, 0, MPI_COMM_WORLD);
	}
	else
	{
		MPI_Send(&numOfHelpersToWorkOn, 1, MPI_INT, sendTo, 0, MPI_COMM_WORLD);

		double *currentMaxDistance;
		vector_t *currentSumPointsLocation;
		int *currentNumOfPoints;
		int *currentPointsIDs;

		int i;
		for (i = 0; i < numOfHelpersToWorkOn; i++)
		{
			currentMaxDistance = &(helpers[i].maxDistancePoint);
			currentSumPointsLocation = &(helpers[i].sumPointLocation);
			currentNumOfPoints = &(helpers[i].numOfPoints);
			currentPointsIDs = helpers[i].pointsIDs;

			MPI_Send(currentMaxDistance, 1, MPI_DOUBLE, sendTo, 0, MPI_COMM_WORLD);
			MPI_Send(currentSumPointsLocation, NUM_OF_ELEMENTS_IN_VECTOR, MPI_DOUBLE, sendTo, 0, MPI_COMM_WORLD);
			MPI_Send(currentNumOfPoints, 1, MPI_INT, sendTo, 0, MPI_COMM_WORLD);
			if (*currentNumOfPoints != 0)
			{
				MPI_Send(currentPointsIDs, (*currentNumOfPoints), MPI_INT, sendTo, 0, MPI_COMM_WORLD);
			}
		}
	}
}

void master_copy_diameters_from_its_own_helpers(cluster_t *clusters, clusterHelper_t *helpers, int helpersCount)
{
	int i;
	for (i = 0; i < helpersCount; i++)
	{
		clusters[i].maxDistancePoint = helpers[i].maxDistancePoint;
	}
}

void master_recv_clusters_diameters(cluster_t *clusters, int clustersCount, int numProcs, int numOfHelpersToWorkOn, double *diameterBuffer)
{	

	// special case last proc/machine calc all the remainder (remainder is the whole)
	if (numOfHelpersToWorkOn == 0)
	{
		recv_diameters(clusters, numProcs, clustersCount, numOfHelpersToWorkOn, diameterBuffer);
	}
	else
	{
		int i;
		for (i = 1; i < numProcs; i++) // i = 1 because master doesnt require to recv from himself
		{
			recv_diameters(clusters, numProcs ,clustersCount, numOfHelpersToWorkOn, diameterBuffer);
		}
	}
}

static void recv_diameters(cluster_t *clusters, int numProcs ,int helpersCount, int recvAmount, double *diameterBuffer)
{

	MPI_Status status;

	// understand who to get data from
	int argForRecv = 0; // argument only for MPI API
	MPI_Recv(&argForRecv, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
	int recvFrom = status.MPI_SOURCE;
	
	//calculate how much data to recv (last proc/machine should have the remainder)
	int bufferRecv = recvAmount;
	if (recvFrom == numProcs - 1)
	{
		// if 3 procs and 10 helpers, 2 first procs recv 3 helpers and the 3rd will recv 4 helpers
		bufferRecv += helpersCount - numProcs * recvAmount; 
	}

	MPI_Recv(diameterBuffer, bufferRecv, MPI_DOUBLE, recvFrom, 0, MPI_COMM_WORLD, &status);

	//calculate where in the cluster array to recv from
	int initPos = recvAmount * recvFrom;

	// special case where one proc recv all then it should start from the starts
	if (helpersCount == recvAmount)
		initPos = 0;

	int i;
	for (i = 0; i < bufferRecv; i++)
	{
		clusters[initPos + i].maxDistancePoint = diameterBuffer[i];
	}
}
