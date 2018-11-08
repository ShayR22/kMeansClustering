
#include "vector_header.h"
#include <math.h>
#include <stdio.h>

#define DOUBLE_EPSILON_PRECISION 0.0001

double vector_get_distance(vector_t *v1, vector_t *v2)
{
	double sumOfDeltas = 0;
	double currentDelta;
	int i;
	for (i = 0; i < NUM_OF_ELEMENTS_IN_VECTOR; i++)
	{
		currentDelta = v1->axis[i] - v2->axis[i];
		sumOfDeltas += currentDelta * currentDelta;
	}

	return sqrt(sumOfDeltas);
}

void vector_add_vector(vector_t *v, vector_t *toBeAdd)
{
	int i;
	for (i = 0; i < NUM_OF_ELEMENTS_IN_VECTOR; i++)
	{
		v->axis[i] += toBeAdd->axis[i];
	}
}

void vector_divide_vector(vector_t *v, double divide)
{
	int i;
	for (i = 0; i < NUM_OF_ELEMENTS_IN_VECTOR; i++)
	{
		v->axis[i] = v->axis[i] / divide;
	}
}

int vector_is_equal(vector_t *v1, vector_t *v2)
{
	int i;
	for (i = 0; i < NUM_OF_ELEMENTS_IN_VECTOR; i++)
	{
		if (fabs(v1->axis[i] -v2->axis[i]) > DOUBLE_EPSILON_PRECISION)
		{
			return 0;
		}
	}
	return 1;
}



void vector_print_vector(vector_t *vector)
{
	printf("x: %lf, y: %lf, z: %lf\n", vector->axis[0], vector->axis[1], vector->axis[2]);
	fflush(stdout);
}