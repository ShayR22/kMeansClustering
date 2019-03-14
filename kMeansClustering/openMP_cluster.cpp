#include "openMP_cluster.h"
#include "openMP_general.h"
#include <omp.h>
#include <stdlib.h>

static void openMpaddToEachPointItsNearestCluster(point_t *point, cluster_t *clusters, int clusterCount);
static void openMPAddToEachHelperItsPoints(clusterHelper_t *helper, int helperID, point_t *points, int pointsCount);

static void calculateClusterHelperPointsLocationSum(clusterHelper_t *helper, point_t *points);
static int updateCenter(cluster_t *cluster, clusterHelper_t *helper);
static void calc_helper_diameter(clusterHelper_t *helper, point_t *points);

static double calcClusterQM(cluster_t *clusters, int clusterAt);

void openMP_add_to_each_point_its_nearest_cluster(cluster_t *clusters, int k, point_t *points, int pointsCount)
{

	#pragma omp parallel
	{
		int threadWorking, range, id, start, end;
		openMP_calc_pragma_range(pointsCount, &threadWorking, &range, &id, &start, &end);
		int i;
		for (i = start; i < end; i++)
		{
			openMpaddToEachPointItsNearestCluster(&(points[i]), clusters, k);
		}
	}
}

static void openMpaddToEachPointItsNearestCluster(point_t *point, cluster_t *clusters, int clusterCount)
{
	cluster_t *nearestCluster = NULL;
	double minDistance = INT_MAX;

	double currentDistance;
	vector_t *currentClusterCenter;

	int belongToo;
	int i;
	for (i = 0; i < clusterCount; i++)
	{
		currentClusterCenter = &(clusters[i].center);
		currentDistance = vector_get_distance(&(point->location), currentClusterCenter);
		if (currentDistance < minDistance)
		{
			minDistance = currentDistance;
			nearestCluster = &(clusters[i]);
			belongToo = i;
		}
	}
	point->clusterBelongTo = belongToo;
}

void openMP_add_to_each_clusterHelper_its_points(clusterHelper_t *helpers, int k, point_t *pointsWorkLoad, int pointCount)
{
	#pragma omp parallel
	{
		int threadWorking, range, id, start, end;
		openMP_calc_pragma_range(k, &threadWorking, &range, &id, &start, &end);
		int i;
		for (i = start; i < end; i++)
		{
			openMPAddToEachHelperItsPoints(&(helpers[i]), i, pointsWorkLoad, pointCount);
		}
	}
}

static void openMPAddToEachHelperItsPoints(clusterHelper_t *helper, int helperID, point_t *points, int pointsCount)
{

	int *helperNumPoints = &(helper->numOfPoints);
	int *IDS = helper->pointsIDs;

	int i;
	for (i = 0; i < pointsCount; i++)
	{
		if (points[i].clusterBelongTo == helperID)
		{
			IDS[(*helperNumPoints)++] = points[i].id;
		}
	}
}

void openMP_reset_clusterHelpers(clusterHelper_t *clusterHelpers, int clusterHelperCount)
{
	#pragma omp parallel
	{
		int threadWorking, range, id, start, end;
		openMP_calc_pragma_range(clusterHelperCount, &threadWorking, &range, &id, &start, &end);
		int i;
		for (i = start; i < end; i++)
		{
			clusterHelpers[i].numOfPoints = 0;
			clusterHelpers[i].sumPointLocation = { 0,0,0 };
			clusterHelpers[i].diameter = 0;
		}
	}
}

void openMP_calc_clusterHelpers_points_location_sum(clusterHelper_t *helpers, int helpersCount, point_t *points)
{
	#pragma omp parallel
	{
		int threadWorking, range, id, start, end;
		openMP_calc_pragma_range(helpersCount, &threadWorking, &range, &id, &start, &end);
		int i;
		for (i = start; i < end; i++)
		{
			calculateClusterHelperPointsLocationSum(&(helpers[i]), points);
		}
	}
}

static void calculateClusterHelperPointsLocationSum(clusterHelper_t *helper, point_t *points)
{
	int i;
	vector_t *sumLocation = &(helper->sumPointLocation);
	vector_t *currentPointLocation;
	int currentPointID;
	for (i = 0; i < helper->numOfPoints; i++)
	{
		//extract the i'th place id from the clusterHelper and add the point's location of that id
		currentPointID = helper->pointsIDs[i];
		currentPointLocation = &(points[currentPointID].location);
		vector_add_vector(sumLocation, currentPointLocation);
	}
}

int openMP_is_center_changes(cluster_t *clusters, clusterHelper_t *helpers, int k)
{
	int anyChange = 0;
	#pragma omp parallel
	{
		int threadWorking, range, id, start, end;
		openMP_calc_pragma_range(k, &threadWorking, &range, &id, &start, &end);
		int i;
		for (i = start; i < end; i++)
		{
			if (updateCenter(&(clusters[i]), &(helpers[i])))
			{
				anyChange = 1;		
			}
		}
	}
	return anyChange;
}

static int updateCenter(cluster_t *cluster, clusterHelper_t *helper)
{
	vector_t *clusterLocation = &(cluster->center);
	vector_t *helperLocation = &(helper->sumPointLocation);
	if (!vector_is_equal(clusterLocation, helperLocation)) // if not equal update center and return true
	{
		*clusterLocation = *helperLocation;
		return 1;
	}
	return 0;
}

void openMP_calc_clusterHelpers_diamaters(clusterHelper_t *helpers, int helpersCount, point_t *allPoints)
{
	
	#pragma omp parallel
	{
		int threadWorking, range, id, start, end;
		clusterHelper_t *currentHelper;
		openMP_calc_pragma_range(helpersCount, &threadWorking, &range, &id, &start, &end);
		int i;
		for (i = start; i < end; i++)
		{
			currentHelper = &(helpers[i]);
			calc_helper_diameter(currentHelper, allPoints);
		}

	}
}

static void calc_helper_diameter(clusterHelper_t *helper, point_t *points)
{
	double diameter = 0;
	double currentDiameter;

	vector_t *compareA;
	vector_t *compareB;

	//loop on all points and compare each point to all the other points with a second loop.
	int i, j;
	for (i = 0; i < helper->numOfPoints - 1; i++)
	{
		compareA = &(points[helper->pointsIDs[i]].location);
		for (j = i + 1; j < helper->numOfPoints; j++)
		{
			compareB = &(points[helper->pointsIDs[j]].location);
			currentDiameter = vector_get_distance(compareA, compareB);
			if (currentDiameter > diameter)
			{
				diameter = currentDiameter;
			}
		}
	}
	helper->diameter = diameter;
}

double openMP_calc_QM(cluster_t *clusters, int clustersCount)
{
	int threadsNumber = omp_get_max_threads();
	double *qResults = (double*)calloc(threadsNumber, sizeof(double));
	

	#pragma omp parallel
	{
		int threadWorking, range, id, start, end;
		openMP_calc_pragma_range(clustersCount, &threadWorking, &range, &id, &start, &end);
		
		int i;
		for (i = start; i < end; i++)
		{
			qResults[id] += calcClusterQM(clusters, i);
		}
	}
	
	double q = 0;
	int i;
	for (i = 0; i < threadsNumber; i++)
	{	
		q += qResults[i];
	}
	free(qResults);
	q /= clustersCount * (clustersCount - 1);


	return q;

}

static double calcClusterQM(cluster_t *clusters, int clusterPosCalc)
{
	double clusterPosCalcDiameter = clusters[clusterPosCalc].diameter;
	vector_t *clusterPosCalcCenter = &(clusters[clusterPosCalc].center);

	double currentDiameter;
	vector_t *currentCenter;
	double currentClustersDistance;

	double qm = 0;
	int i;
	for (i = 0; i < clusterPosCalc; i++)
	{
		currentDiameter = clusters[i].diameter;
		currentCenter = &(clusters[i].center);
		currentClustersDistance = vector_get_distance(clusterPosCalcCenter, currentCenter);

		qm += (clusterPosCalcDiameter + currentDiameter) / currentClustersDistance;
	}
	return qm;
}


