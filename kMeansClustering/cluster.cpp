#include "cuda.h"
#include "cluster.h"
#include "vector.h"
#include "openMP_cluster.h"
#include "openMP_general.h"
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static int get_number_to_work_on_without_zero_or_count(int count, double devision);

static void print_cluster(cluster_t *cluster);

static void copy_clusterHelper(clusterHelper_t *copyFrom, clusterHelper_t *copyTo);

static void calc_helper_diameter(clusterHelper_t *helper, point_t *points);

// Public methods
void cluster_init(cluster_t *cluster, vector_t *center)
{
	cluster->center = *center;
}

void clusterHelper_allocate(clusterHelper_t *helpers, int helpersCount, int numPoints)
{
	int i;
	for (i = 0; i < helpersCount; i++)
	{
		helpers[i].pointsIDs = (int*)calloc(numPoints, sizeof(int));
	}
}

void cluster_reset_clusterHelpers(clusterHelper_t *helpers, int helpersCount)
{
	openMP_reset_clusterHelpers(helpers, helpersCount);
}

void cluster_calc_clusterHelper_points_location_sum(clusterHelper_t *clusterHelpers, int k, point_t *points)
{
	openMP_calc_clusterHelpers_points_location_sum(clusterHelpers, k, points);
}

void cluster_add_to_each_point_its_nearest_cluster(cluster_t *clusters, int k, point_t *points, int pointsCount)
{
	double static cudaPointDevisionClusterToPoint = 0.5;
	int numOfCudaPoints = get_number_to_work_on_without_zero_or_count(pointsCount , cudaPointDevisionClusterToPoint);
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
			openMP_add_to_each_point_its_nearest_cluster(clusters, k, openMPpoints, numOfOpenMPpoints);
			openMPTimeEnd = MPI_Wtime();
		}
		#pragma omp section
		{	
			cudaTimeStart = MPI_Wtime();
			if (!cuda_add_to_each_point_its_nearest_cluster(cudaPoints, numOfCudaPoints, clusters, k))
			{
				printf("CUDA FAILED in distributing to each point its cluster \n");
				fflush(stdout);
			}
			cudaTimeEnd = MPI_Wtime();
		}
	}
	openMP_calc_new_deviation_based_on_time(cudaTimeStart, cudaTimeEnd, openMPTimeStart, openMPTimeEnd, &cudaPointDevisionClusterToPoint);
}

void cluster_add_to_each_clusterHelper_its_points(clusterHelper_t *helpers, int k, point_t *pointsWorkLoad, int pointCount)
{
	openMP_add_to_each_clusterHelper_its_points(helpers, k, pointsWorkLoad, pointCount);
}

int cluster_is_center_changes(cluster_t *clusters, clusterHelper_t *helpers, int k)
{
	return openMP_is_center_changes(clusters, helpers, k);
}


void cluster_copy_clusterHelpers(clusterHelper_t *copyFrom, clusterHelper_t *copyTo, int numOfHelpers)
{
	clusterHelper_t *from;
	clusterHelper_t *to;
	int i;
	for (i = 0; i < numOfHelpers; i++)
	{
		from = &(copyFrom[i]);
		to = &(copyTo[i]);
		copy_clusterHelper(from, to);
	}
}

static void copy_clusterHelper(clusterHelper_t *copyFrom, clusterHelper_t *copyTo)
{
	copyTo->diameter = copyFrom->diameter;
	copyTo->numOfPoints = copyFrom->numOfPoints;
	copyTo->sumPointLocation = copyFrom->sumPointLocation;
}

void cluster_calculate_helpers_diameter(clusterHelper_t *helpers, int helpersCount, point_t *points, int pointsCount)
{
	//cuda_calc_each_helper_diamater(points, pointsCount, helpers, helpersCount);
	//openMP_calc_clusterHelpers_diamaters(helpers, helpersCount, points);

	if (helpersCount < 2)
	{
		cuda_calc_each_helper_diamater(points, pointsCount, helpers, helpersCount);
	}
	else
	{
		double static cudaHelperDevisionDiameter = 0.5;
		int numOfCudaHelpers = get_number_to_work_on_without_zero_or_count(helpersCount, cudaHelperDevisionDiameter);
		int numOfOpenMPHelpers = helpersCount - numOfCudaHelpers;

		clusterHelper_t *cudaHelpers = helpers; // give a different name for easy understanding in later code
		clusterHelper_t *openMPHelpers = &(helpers[numOfCudaHelpers]);

		double openMPTimeStart;
		double openMPTimeEnd;
		double cudaTimeStart;
		double cudaTimeEnd;

		#pragma omp parallel sections
		{
			#pragma omp section
			{
				openMPTimeStart = MPI_Wtime();
				openMP_calc_clusterHelpers_diamaters(openMPHelpers, numOfOpenMPHelpers, points);
				openMPTimeEnd = MPI_Wtime();
			}
			#pragma omp section
			{
				cudaTimeStart = MPI_Wtime();
				cuda_calc_each_helper_diamater(points, pointsCount, cudaHelpers, numOfCudaHelpers);
				cudaTimeEnd = MPI_Wtime();
			}
		}
		openMP_calc_new_deviation_based_on_time(cudaTimeStart, cudaTimeEnd, openMPTimeStart, openMPTimeEnd, &cudaHelperDevisionDiameter);
	}

}

double cluster_calculate_quality_measure(cluster_t *clusters, int clustersCount)
{
	return openMP_calc_QM(clusters, clustersCount);
}

static int get_number_to_work_on_without_zero_or_count(int count, double devision)
{
	int numToWorkOn = (int)(count * devision);
	if (numToWorkOn == 0)
	{
		numToWorkOn++;
	}
	else if (numToWorkOn == count) {
		numToWorkOn--;
	}
	return numToWorkOn;
}

void cluster_free_clusterHelpers_array(clusterHelper_t *helpers, int helpersCount)
{
	int i;
	for (i = 0; i < helpersCount; i++)
	{
		free(helpers[i].pointsIDs);
	}
	free(helpers);
}

// ================================ PRINTING METHODS BELOW ========================================

void cluster_print_clusters(cluster_t *clusters, int clustersCount)
{
	int i;
	for (i = 0; i < clustersCount; i++)
	{
		print_cluster(&(clusters[i]));
	}
	printf("\n");
}

static void print_cluster(cluster_t *cluster)
{
	vector_print_vector(&(cluster->center));
	printf("cluster Diameter is %lf\n\n", cluster->diameter);
	fflush(stdout);
}


MPI_Datatype cluster_create_cluster_mpi_struct()
{
	// 2 vectors (location and speed) + cluster belongs too (int) + distance from cluster (double) 
	int numOfItems = NUMBER_OF_ELEMENTS_IN_CLUSTER;
	int blocklengths[NUMBER_OF_ELEMENTS_IN_CLUSTER] = { NUM_OF_ELEMENTS_IN_VECTOR,1 };
	MPI_Datatype types[NUMBER_OF_ELEMENTS_IN_CLUSTER] = { MPI_DOUBLE, MPI_DOUBLE};

	MPI_Aint offsets[NUMBER_OF_ELEMENTS_IN_CLUSTER];
	offsets[0] = offsetof(cluster_t, center);
	offsets[1] = offsetof(cluster_t, diameter);
	

	MPI_Datatype mpiClusterType;
	MPI_Type_create_struct(numOfItems, blocklengths, offsets, types, &mpiClusterType);
	MPI_Type_commit(&mpiClusterType);

	return mpiClusterType;
}
