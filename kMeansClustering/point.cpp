#include "point.h"
#include "openMP_general.h"
#include "openMP_point.h"
#include "cuda.h"
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stddef.h>
#include <omp.h>

void point_advance_points_by_dt(point_t*points, int pointsCount, double deltaT)
{
	double static cudaPointDevisionClusterToPoint = 0.9;
	int numOfCudaPoints = (int)(pointsCount * cudaPointDevisionClusterToPoint);
	if (numOfCudaPoints == 0)
	{
		numOfCudaPoints++;
	}
	else if (numOfCudaPoints == pointsCount)
	{
		numOfCudaPoints--;
	}
	int numOfOpenMPpoints = pointsCount - numOfCudaPoints;
	point_t *cudaPoints = points; // give a different name for easy understanding in later code
	point_t *openMPpoints = &(points[numOfCudaPoints]);

	double openMPTimeStart;
	double openMPTimeEnd;
	double cudaTimeStart;
	double cudaTimeEnd;

#pragma omp parallel sections
	{
		#pragma omp section
		{
			openMPTimeStart = MPI_Wtime();
			openMP_advance_points_by_deltaT(openMPpoints, numOfOpenMPpoints, deltaT);
			openMPTimeEnd = MPI_Wtime();
		}
		#pragma omp section
		{
			cudaTimeStart = MPI_Wtime();
			if (!cuda_increment_points(cudaPoints, numOfCudaPoints, deltaT))
			{
				printf("CUDA FAILED in distributing to each point its cluster \n");
				fflush(stdout);
			}
			cudaTimeEnd = MPI_Wtime();
		}
	}
	openMP_calc_new_deviation_based_on_time(cudaTimeStart, cudaTimeEnd, openMPTimeStart, openMPTimeEnd, &cudaPointDevisionClusterToPoint);
}

void point_advance_point_by_dt(point_t *point, double deltaT)
{
	int i;
	for (i = 0; i < NUM_OF_ELEMENTS_IN_VECTOR; i++)
	{
		point->location.axis[i] += point->speed.axis[i] * deltaT;
	}
}

void point_add_point_location(point_t *p1, point_t *p2)
{
	vector_add_vector(&(p1->location), &(p2->location));
}

double point_get_distance(point_t *p1, point_t *p2)
{
	return vector_get_distance(&(p1->location), &(p2->location));
}


// ======================================= PRINTING METHODS ===========================================
void point_print_points(point_t *points, int pointsCount)
{
	int i;
	for (i = 0; i < pointsCount; i++)
	{
		point_print_point(&(points[i]));
	}
}

void point_print_point(point_t *p)
{
	char axis[] = { 'X','Y','Z' };
	printf("points %d belongs to: %d and its location:\n", p->id, p->clusterBelongTo);
	int i;
	for (i = 0; i < NUM_OF_ELEMENTS_IN_VECTOR; i++)
	{
		printf("%c: %lf\n", axis[i], p->location.axis[i]);
		fflush(stdout);
	}
}


MPI_Datatype point_create_mpi_struct()
{
	// id + 2 vectors (location and speed) + cluster belongs too (int) + distance from cluster (double) 
	int numOfItems = NUM_OF_ELEMENTS_IN_POINT;
	int blocklengths[NUM_OF_ELEMENTS_IN_POINT] = { 1, NUM_OF_ELEMENTS_IN_VECTOR, NUM_OF_ELEMENTS_IN_VECTOR, 1};
	MPI_Datatype types[NUM_OF_ELEMENTS_IN_POINT] = { MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT };

	MPI_Aint offsets[NUM_OF_ELEMENTS_IN_POINT];
	offsets[0] = offsetof(point_t, id);
	offsets[1] = offsetof(point_t, location);
	offsets[2] = offsetof(point_t, speed);
	offsets[3] = offsetof(point_t, clusterBelongTo);
	
	MPI_Datatype mpiPointType;
	MPI_Type_create_struct(numOfItems, blocklengths, offsets, types, &mpiPointType);
	MPI_Type_commit(&mpiPointType);

	return mpiPointType;
}

