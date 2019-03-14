#include "openMP_point.h"
#include "openMP_general.h"
#include <omp.h>


void openMP_advance_points_by_deltaT(point_t *points, int pointsCount, double deltaT)
{
	#pragma omp parallel
	{
		int threadWorking, range, id, start, end;
		openMP_calc_pragma_range(pointsCount, &threadWorking, &range, &id, &start, &end);
		int i;
		for (i = start; i < end; i++)
		{
			point_advance_point_by_dt(&(points[i]), deltaT);
		}
	}
}