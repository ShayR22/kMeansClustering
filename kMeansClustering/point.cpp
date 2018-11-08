
#include "point_header.h"
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stddef.h>

static void incrementLocation(point_t *point, double deltaT);
static void advancePointInDt(point_t *point, double deltaT);


void point_incrementByDT(point_t *points, int pointsCount)
{

}


void point_increment_location(point_t *points, int pointsCount, double deltaT)
{
	point_t *currentPoint;
	int i;
	for (i = 0; i < pointsCount; i++)
	{
		currentPoint = &(points[i]);
		incrementLocation(currentPoint, deltaT);
	}
}

static void incrementLocation(point_t *p, double deltaT)
{
	int i;
	for (i = 0; i < NUM_OF_ELEMENTS_IN_VECTOR; i++)
	{
		p->location.axis[i] += p->speed.axis[i] * deltaT;
	}
}


void point_advance_points_by_dt(point_t*points, int pointsCount, double deltaT)
{
	int i;
	for (i = 0; i < pointsCount; i++)
	{
		advancePointInDt(&(points[i]), deltaT);
	}
}

static void advancePointInDt(point_t *point, double deltaT)
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

