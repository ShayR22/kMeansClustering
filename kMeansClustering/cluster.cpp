#include "cudaHeader.h"
#include "cluster_header.h"
#include "vector_header.h"
#include "openMP_cluster_header.h"
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>


static void printCluster(cluster_t *cluster);
static void printClusterHelper(clusterHelper_t *helper);

static void copy_clusterHelper(clusterHelper_t *copyFrom, clusterHelper_t *copyTo);

static void calc_helper_diameter(clusterHelper_t *helper, point_t *points);

// Public methods
void cluster_init(cluster_t *cluster, vector_t *center)
{
	static int ID = 0;
	cluster->center = *center;
	cluster->id = ID++;
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

void cluster_calc_clusterHelper_points_location_sum(cluster_t *clusters, clusterHelper_t *clusterHelpers, int k, point_t *points)
{
	openMP_calc_clusterHelpers_points_location_sum(clusterHelpers, k, points);
}


void cluster_add_to_each_point_its_nearest_cluster(cluster_t *clusters, int k, point_t *points, int pointsCount)
{
	int numOfCudaPoints = pointsCount / 2;
	int numOfOpenMPpoints = pointsCount - numOfCudaPoints;
	point_t *cudaPoints = points; // give a different name for easy understanding in later code
	point_t *openMPpoints = &(points[numOfOpenMPpoints]);
	
	//TODO resume pragma sections
	openMP_add_to_each_point_its_nearest_cluster(clusters, k, points, pointsCount);

	/*#pragma omp parallel sections
	{
		#pragma omp section
		{
			openMP_add_to_each_point_its_nearest_cluster(clusters, k, openMPpoints, numOfOpenMPpoints);
		}
		#pragma omp section
		{
			if (!cuda_add_to_each_point_its_nearest_cluster(cudaPoints, numOfCudaPoints, clusters, k))
			{
				printf("CUDA FAILED in distributing to each point its cluster");
				fflush(stdout);
			}
		}
	}*/
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
	copyTo->maxDistancePoint = copyFrom->maxDistancePoint;
	copyTo->numOfPoints = copyFrom->numOfPoints;
	copyTo->sumPointLocation = copyFrom->sumPointLocation;
}


//TODO make this method with openMP and CUDA
void cluster_calculate_helpers_diameter(clusterHelper_t *helpers, int helpersCount, point_t *points, int pointsCount)
{

	openMP_calc_clusterHelpers_diamaters(helpers, helpersCount, points);
	//if (helpersCount < 2)
	//{
	//	openMP_calc_clusterHelpers_diamaters(helpers, helpersCount, points);
	//}
	//else
	//{
	//	int numOfCudaHelpers = helpersCount / 2;
	//	int numOfOpenMPHelpers = helpersCount - numOfCudaHelpers;
	//	clusterHelper_t *cudaHelpers = helpers; // give a different name for easy understanding in later code
	//	clusterHelper_t *openMPHelpers = &(helpers[numOfOpenMPHelpers]);

	//	#pragma omp parallel sections
	//	{
	//		#pragma omp section
	//		{
	//			openMP_calc_clusterHelpers_diamaters(openMPHelpers, numOfOpenMPHelpers, points);
	//		}
	//		#pragma omp section
	//		{
	//			cuda_calc_each_helper_diamater(points, pointsCount, cudaHelpers, numOfCudaHelpers);
	//		}
	//	}
	//}

}

double cluster_calculate_quality_measure(cluster_t *clusters, int clustersCount)
{
	return openMP_calc_QM(clusters, clustersCount);
}


void cluster_print_clusters(cluster_t *clusters, int clustersCount)
{
	int i;
	for (i = 0; i < clustersCount; i++)
	{
		printCluster(&(clusters[i]));
	}
	printf("\n");
}

static void printCluster(cluster_t *cluster)
{
	vector_print_vector(&(cluster->center));
	printf("cluster Diameter is %lf\n\n", cluster->maxDistancePoint);
	fflush(stdout);
}

void cluster_print_clusterHelpers(clusterHelper_t *helpers, int helpersCount)
{
	int i;
	for (i = 0; i < helpersCount; i++)
	{
		printf("helper number %d\n", i);
		fflush(stdout);
		printClusterHelper(&(helpers[i]));
	}

}

static void printClusterHelper(clusterHelper_t *helper)
{
	int i;
	printf("helper has %d points and there are:\n", helper->numOfPoints);
	fflush(stdout);
	for (i = 0; i < helper->numOfPoints; i++)
	{
		printf("%d, ", helper->pointsIDs[i]);
		fflush(stdout);
	}
	printf("\nsumOfVecotrs is:\n");
	fflush(stdout);
	vector_print_vector(&(helper->sumPointLocation));
	printf("maxDistance = %lf\n\n", helper->maxDistancePoint);
	fflush(stdout);
}

MPI_Datatype cluster_create_cluster_mpi_struct()
{
	// id + 2 vectors (location and speed) + cluster belongs too (int) + distance from cluster (double) 
	int numOfItems = NUMBER_OF_ELEMENTS_IN_CLUSTER;
	int blocklengths[NUMBER_OF_ELEMENTS_IN_CLUSTER] = { 1, NUM_OF_ELEMENTS_IN_VECTOR,1 };
	MPI_Datatype types[NUM_OF_ELEMENTS_IN_POINT] = { MPI_INT, MPI_DOUBLE, MPI_DOUBLE};

	MPI_Aint offsets[NUM_OF_ELEMENTS_IN_POINT];
	offsets[0] = offsetof(cluster_t, id);
	offsets[1] = offsetof(cluster_t, center);
	offsets[2] = offsetof(cluster_t, maxDistancePoint);
	

	MPI_Datatype mpiClusterType;
	MPI_Type_create_struct(numOfItems, blocklengths, offsets, types, &mpiClusterType);
	MPI_Type_commit(&mpiClusterType);

	return mpiClusterType;
}

MPI_Datatype cluster_create_clusterHelper_mpi_struct()
{
	// id + 2 vectors (location and speed) + cluster belongs too (int) + distance from cluster (double) 
	int numOfItems = NUMBER_OF_ELEMENTS_IN_CLUSTER_HELPER;
	int blocklengths[NUMBER_OF_ELEMENTS_IN_CLUSTER_HELPER] = {1,NUM_OF_ELEMENTS_IN_VECTOR, 1};
	MPI_Datatype types[NUMBER_OF_ELEMENTS_IN_CLUSTER_HELPER] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};

	MPI_Aint offsets[NUMBER_OF_ELEMENTS_IN_CLUSTER_HELPER];
	offsets[0] = offsetof(clusterHelper_t, maxDistancePoint);
	offsets[1] = offsetof(clusterHelper_t, sumPointLocation);
	offsets[2] = offsetof(clusterHelper_t, numOfPoints);

	MPI_Datatype mpiClusterHelperType;
	MPI_Type_create_struct(numOfItems, blocklengths, offsets, types, &mpiClusterHelperType);
	MPI_Type_commit(&mpiClusterHelperType);

	return mpiClusterHelperType;
}



