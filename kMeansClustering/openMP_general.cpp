#include "openMP_general.h"
#include <math.h>
#include <omp.h>

#include <stdio.h>

#define DEVISION_COEFFICIENT 10.0

static void calc_new_deviation_based_on_time_helper(double *currentDevision, double deviation);

void openMP_calc_pragma_range(int count, int *threadWorking, int *range, int *id, int *start, int *end)
{
	*threadWorking = omp_get_num_threads();
	*range = count / *threadWorking;
	*id = omp_get_thread_num();
	*start = (*id) * (*range);
	*end = *start + *range;
	if (*id == *threadWorking - 1) // last thread
	{
		*end += (count - *end); // add remainder
	}
}

void openMP_calc_new_deviation_based_on_time(double t1Start, double t1End, double t2Start, double t2End, double *t1currentDevision)
{
	double t1 = t1End - t1Start;
	double t2 = t2End - t2Start;
	double totalTime = t1+t2;
	double t1Precent = t1 / totalTime;
	calc_new_deviation_based_on_time_helper(t1currentDevision, t1Precent);
}

static void calc_new_deviation_based_on_time_helper(double *currentDevision, double deviation)
{
	double absDelta = fabs(0.5 - deviation) / DEVISION_COEFFICIENT;
	if (deviation > 0.5)
	{
		*currentDevision -= absDelta;
	}
	else
	{
		*currentDevision += absDelta;
	}
}
