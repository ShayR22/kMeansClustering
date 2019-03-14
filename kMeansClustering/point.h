#ifndef POINT_HEADER
#define POINT_HEADER

#include "vector.h"
#include <mpi.h>

#define NUM_OF_ELEMENTS_IN_POINT 4

typedef struct
{
	int id;
	vector_t location; 
	vector_t speed;
	int clusterBelongTo;
} point_t;


void point_advance_points_by_dt(point_t*points, int pointsCount, double deltaT);
void point_advance_point_by_dt(point_t *point, double deltaT);

void point_add_point_location(point_t *p1, point_t *p2);
double point_get_distance(point_t *p1, point_t *p2);

/*
create and commit point_t mpi sturct for later use

Parameters: void

Returns:
- MPI_DataType
*/
MPI_Datatype point_create_mpi_struct();

void point_print_points(point_t *points, int pointsCount);
void point_print_point(point_t *p);



#endif // ! POINT_HEADER

