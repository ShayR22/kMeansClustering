
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cudaHeader.h"
#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define CPU_TO_GPU cudaMemcpyHostToDevice
#define GPU_TO_CPU cudaMemcpyDeviceToHost
#define WRAP_SIZE 32
#define NUM_OF_THREADS_PER_BLOCK 256

#define RETURN_IF_NOT_TRUE(truth)if(!truth) {return 0;}

static int setDevice();

static int allocate_GPU_space(void** data, int dataCount, int dataTypeSize);
static int free_GPU_space(void* data);

static int copy_to_GPU_memory(void* copyTo, void* copyFrom, int count, int typeSize);
static int copy_from_GPU_memory(void* copyTo, void* copyFrom, int count, int typeSize);

static int checkError(cudaError error);

__device__ double getDistance(vector_t *v1, vector_t *v2)
{
	int i;
	double sumOfDeltaSquares = 0;
	double delta;
	for (i = 0; i < NUM_OF_ELEMENTS_IN_VECTOR; i++)
	{
		delta = v1->axis[i] - v2->axis[i];
		sumOfDeltaSquares += delta * delta;
	}
	return sqrt(sumOfDeltaSquares);
}


//increment points START
__device__ void advance_point_in_deltaT(point_t *p, double deltaT)
{
	int i;
	for (i = 0; i < NUM_OF_ELEMENTS_IN_VECTOR; i++)
	{
		p->location.axis[i] += p->speed.axis[i] * deltaT;
	}
}

__global__ void increment_points(point_t *devPoints, int pointsCount, double deltaT)
{
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	if (index < pointsCount)
	{
		point_t *point = &(devPoints[index]);
		advance_point_in_deltaT(point, deltaT);
	}
}


int cuda_increment_points(point_t *points, int pointsCount, double deltaT)
{
	point_t *devPoints;
	setDevice();
	// allocate space and then copy points in to it
	RETURN_IF_NOT_TRUE(allocate_GPU_space((void**)&devPoints, pointsCount, sizeof(point_t)));
	RETURN_IF_NOT_TRUE(copy_to_GPU_memory(devPoints, points, pointsCount, sizeof(point_t)));

	int numOfBlocks = pointsCount / NUM_OF_THREADS_PER_BLOCK + 1;
	increment_points <<<numOfBlocks, NUM_OF_THREADS_PER_BLOCK >>> (devPoints, pointsCount, deltaT);
	RETURN_IF_NOT_TRUE(copy_from_GPU_memory(points, devPoints, pointsCount, sizeof(point_t)));

	RETURN_IF_NOT_TRUE(free_GPU_space(devPoints));
	return 1;
}
//increment points END


// calcPointsNearestCluster START
__global__ void calcPointNearestCluster(point_t *devPoints, int pointsCount, cluster_t *devClusters, int clustersCount)
{
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	if (index < pointsCount)
	{
		vector_t *pointLocation = &(devPoints[index].location);
		vector_t *currentClusterLocation;
		cluster_t *nearestCluster = NULL;

		double distance = INT_MAX;
		double currentDistance;

		int i;
		for (i = 0; i < clustersCount; i++)
		{
			currentClusterLocation = &(devClusters[i].center);
			currentDistance = getDistance(pointLocation, currentClusterLocation);
			if(currentDistance < distance)
			{
				distance = currentDistance;
				nearestCluster = &(devClusters[i]);
			}
		}
		devPoints[index].clusterBelongTo = nearestCluster->id;
	}
}

int cuda_add_to_each_point_its_nearest_cluster(point_t *points, int pointsCount, cluster_t *clusters, int clustersCount)
{
	setDevice();
	
	cluster_t *devClusters;
	point_t *devPoints;
		
	// allocate GPU space
	RETURN_IF_NOT_TRUE(allocate_GPU_space((void**)&devClusters, clustersCount, sizeof(cluster_t)));
	RETURN_IF_NOT_TRUE(allocate_GPU_space((void**)&devPoints, pointsCount, sizeof(point_t)));
	
	//copy clusters and point to GPU's memory
	RETURN_IF_NOT_TRUE(copy_to_GPU_memory(devClusters, clusters, clustersCount, sizeof(cluster_t)));
	RETURN_IF_NOT_TRUE(copy_to_GPU_memory(devPoints, points, pointsCount, sizeof(point_t)));
		
	int numBlocks = pointsCount / NUM_OF_THREADS_PER_BLOCK + 1;
	calcPointNearestCluster <<<numBlocks, NUM_OF_THREADS_PER_BLOCK>>> (devPoints, pointsCount, devClusters, clustersCount);
		
	RETURN_IF_NOT_TRUE(copy_from_GPU_memory(points, devPoints, pointsCount, sizeof(point_t)));
		
	RETURN_IF_NOT_TRUE(free_GPU_space(devClusters));
	RETURN_IF_NOT_TRUE(free_GPU_space(devPoints));
	
	return 1;
}
// calcPointsNearestCluster END

// calcDiamaters START

__global__ void calc_diameters(clusterHelper_t *devHelpers, int *devHelpersIDs ,point_t *devPoints)
{
	int helperNumber = blockDim.x * blockIdx.x;
	clusterHelper_t *helper = &(devHelpers[helperNumber]);
	int numPoints = helper->numOfPoints;

	vector_t *v1;
	vector_t *v2;
	point_t *p1;
	point_t *p2;
	int *currentIDS;
	int sumIDs = 0;

	double maxDistance = 0;
	double currentDistance;

	int i, j;
	for (i = 0; i < numPoints - 1; i++)
	{
		currentIDS = &(devHelpersIDs[sumIDs]);
		p1 = &(devPoints[currentIDS[i]]);
		v1 = &(p1->location);
		for (j = i + 1; j < numPoints; j++)
		{
			p2 = &(devPoints[currentIDS[j]]);
			v2 = &(p2->location);
			currentDistance = getDistance(v1, v2);

			if (currentDistance > maxDistance)
			{
				maxDistance = currentDistance;
			}
		}
		sumIDs += devHelpers[i].numOfPoints;
	}
	helper->maxDistancePoint = maxDistance;
}

//TODO this work with small sample data but not with big data
int cuda_calc_each_helper_diamater(point_t *points, int pointsCount, clusterHelper_t *helpers, int helpersCount)
{
	setDevice();
	
	clusterHelper_t *devHelpers;
	point_t *devPoints;
	int *devHelpersIDs;
	
	// allocate GPU space
	RETURN_IF_NOT_TRUE(allocate_GPU_space((void**)&devHelpers, helpersCount, sizeof(clusterHelper_t)));
	RETURN_IF_NOT_TRUE(allocate_GPU_space((void**)&devPoints, pointsCount, sizeof(point_t)));
	
	//copy helpers and point to GPU's memory
	RETURN_IF_NOT_TRUE(copy_to_GPU_memory(devPoints, points, pointsCount, sizeof(point_t)));
	RETURN_IF_NOT_TRUE(copy_to_GPU_memory(devHelpers, helpers, helpersCount, sizeof(clusterHelper_t)));

	RETURN_IF_NOT_TRUE(allocate_GPU_space((void**)&devHelpersIDs, pointsCount, sizeof(int)));


	//copy to devHelpersIDS all points IDS from the helpers;
	int k;
	int numPoints;
	int currentSum = 0;
	clusterHelper_t *currentHelper;
	for (k = 0; k < helpersCount; k++)
	{
		currentHelper = &(helpers[k]);
		numPoints = currentHelper->numOfPoints;
		RETURN_IF_NOT_TRUE(copy_to_GPU_memory(&(devHelpersIDs[currentSum]), currentHelper->pointsIDs, numPoints, sizeof(int)));		
		currentSum += numPoints;
	}

	int numBlocks = helpersCount;
	calc_diameters << <numBlocks, 1 >> > (devHelpers, devHelpersIDs ,devPoints);

	RETURN_IF_NOT_TRUE(copy_from_GPU_memory(helpers, devHelpers, helpersCount, sizeof(clusterHelper_t)));

	RETURN_IF_NOT_TRUE(free_GPU_space(devHelpersIDs));
	RETURN_IF_NOT_TRUE(free_GPU_space(devHelpers));
	RETURN_IF_NOT_TRUE(free_GPU_space(devPoints));

	return 1;
}

// calcDiamaters END

static int setDevice()
{
	int device = 0;
	cudaError_t error = cudaSetDevice(device);
	return checkError(error);
}

static int allocate_GPU_space(void** data, int dataCount, int dataTypeSize)
{
	cudaError error = cudaMalloc(data, dataCount * dataTypeSize);
	return checkError(error);
}

static int free_GPU_space(void* data)
{
	cudaError error = cudaFree(data);
	return checkError(error);
}

static int copy_to_GPU_memory(void* copyTo, void* copyFrom, int count, int typeSize)
{
	cudaError error = cudaMemcpy(copyTo, copyFrom, count*typeSize, CPU_TO_GPU);
	return checkError(error);
}

static int copy_from_GPU_memory(void *copyTo, void* copyFrom, int count, int typeSize)
{
	cudaError error = cudaMemcpy(copyTo, copyFrom, count*typeSize, GPU_TO_CPU);
	return checkError(error);
}


static int checkError(cudaError error)
{
	if (error != cudaSuccess)
	{
		printf("error number %d occured\n", error);
		fflush(stdout);
		return 0;
	}
	return 1;
}

// calcHelperPoints START
//__global__ void calcHelpersPoints(clusterHelper_t *devHelpers, int *devHelpersPointArray, point_t *devPoints, int pointsCount)
//{
//	int helperNumber = blockDim.x * blockIdx.x;
//	int *numOfPoints = &(devHelpersPointArray[helperNumber * pointsCount]);
//
//	int i;
//	for (i = 0; i < pointsCount; i++)
//	{
//		if (devPoints[i].clusterBelongTo == helperNumber)
//		{
//			numOfPoints[i] = devPoints[i].clusterBelongTo;
//			devHelpers[helperNumber].numOfPoints++;
//		}
//	}
//}
//
//
//int cuda_add_to_each_clusterHelper_its_points(clusterHelper_t *helpers, int helpersCount, point_t *points, int pointsCount)
//{
//	setDevice();
//
//	clusterHelper_t *devHelpers;
//	int *devHelpersPointsArray;
//	point_t *devPoints;
//
//	// allocate GPU space
//	RETURN_IF_NOT_TRUE(allocate_GPU_space((void**)&devHelpers, helpersCount, sizeof(clusterHelper_t)));
//	RETURN_IF_NOT_TRUE(allocate_GPU_space((void**)&devHelpersPointsArray, helpersCount * pointsCount, sizeof(int)));
//	RETURN_IF_NOT_TRUE(allocate_GPU_space((void**)&devPoints, pointsCount, sizeof(point_t)));
//
//	//copy clusters and point to GPU's memory
//	RETURN_IF_NOT_TRUE(copy_to_GPU_memory(devPoints, points, pointsCount, sizeof(point_t)));
//	RETURN_IF_NOT_TRUE(copy_to_GPU_memory(devHelpers, helpers, helpersCount, sizeof(clusterHelper_t)));
//
//
//	int numBlocks = helpersCount;
//	calcHelpersPoints <<<numBlocks, 1>>> (devHelpers, devHelpersPointsArray ,devPoints, pointsCount);
//
//
//
//	RETURN_IF_NOT_TRUE(copy_from_GPU_memory(helpers, devHelpers, helpersCount, sizeof(clusterHelper_t)));
//	int i;
//	for (i = 0; i < helpersCount; i++)
//	{
//		RETURN_IF_NOT_TRUE(copy_from_GPU_memory(helpers[i].pointsIDs, &(devHelpersPointsArray[i*pointsCount]), pointsCount, sizeof(int)));
//	}
//
//	RETURN_IF_NOT_TRUE(free_GPU_space(devHelpers));
//	RETURN_IF_NOT_TRUE(free_GPU_space(devPoints));
//	RETURN_IF_NOT_TRUE(free_GPU_space(devHelpersPointsArray));
//
//	return 1;
//
//}
// calcHelperPoints END


//__global__ void calcPointCluster(cluster_t *devClusters, int clusterCount, point_t *devPoints, int pointCount)
//{
//	int offset = blockIdx.x * blockDim.x + threadIdx.x;
//	if (offset < pointCount) // pointCount doesnt have to be divisble by 256
//	{
//		point_t *point = &(devPoints[offset]);
//		int clusterNumber = -1;
//		double minDistance = INT_MAX;
//		double currentDistance;
//		vector_t *currentClusterCenter;
//
//		//iterate on all cluster for the point and save the cluster number and its distance from the point.
//		int i, j;
//		double sumOfDeltas = 0;
//		double currentDelta;
//		for (i = 0; i < clusterCount; i++)
//		{
//			currentClusterCenter = &(devClusters[i].center);
//
//			//TODO put in method
//			//calculate the distance between the point location and the cluster center
//			sumOfDeltas = 0;
//			for (j = 0; j < NUM_OF_ELEMENTS_IN_VECTOR; j++)
//			{
//				currentDelta = currentClusterCenter->axis[j] - point->location.axis[j];
//				sumOfDeltas += currentDelta * currentDelta;
//			}
//			currentDistance = sqrt(sumOfDeltas);
//
//			//update relevant data for the point if 
//			if (currentDistance < minDistance)
//			{
//				minDistance = currentDistance;
//				clusterNumber = i;
//			}
//		}
//
//		point->clusterNumber = clusterNumber;
//		point->distanceFromCluster = minDistance;
//	}
//}
//
//
//
//int cuda_add_nearest_point_to_cluster(cluster_t *clusters, int clusterCount, point_t *points, int pointCount)
//{
//	//TODO fix docu below
//	/*
//	 need to go over each point and assign to the right cluster
//    each thread can represent a point, in order for this to work
//	beside the clusters and point we will need an array of int in the size of pointCount.
//	the array will represent to which cluster the point in the id of it "i" belongs meaning
//	int belongCluster[pointCount];
//	belongCluster[0] = 4 => means point 0 belong to cluster 4
//	another array of double in the size of pointCount which will represent the "i" point distance from
//	its corrosponding cluster
//
//	 note: altohugh the threads dont need the whole cluster and only its center the entire cluster will be copy
//	as it easier and faster than extrapolating thier centers and then copy them (neglegible clusterCount should be
//	fairly small).
//
//	 number of points can be to very large therefore the division will be with 256 threads per block 
//	256 was chosen as it is divisible by the wrap size (32) 
//	
//	
//	*/
//
//	setDevice();
//
//	cluster_t *devClusters;
//	point_t *devPoints;
//	
//	// allocate GPU space
//	RETURN_IF_NOT_TRUE(allocate_GPU_space((void**)&devClusters, clusterCount, sizeof(cluster_t)));
//	RETURN_IF_NOT_TRUE(allocate_GPU_space((void**)&devPoints, pointCount, sizeof(point_t)));
//
//	//copy clusters and point to GPU's memory
//	RETURN_IF_NOT_TRUE(copy_to_GPU_memory(devClusters, clusters, clusterCount, sizeof(cluster_t)));
//	RETURN_IF_NOT_TRUE(copy_to_GPU_memory(devPoints, points, pointCount, sizeof(point_t)));
//	
//	int numBlocks = pointCount / NUM_OF_THREADS_PER_BLOCK + 1;
//	calcPointCluster <<<numBlocks, NUM_OF_THREADS_PER_BLOCK>>> (devClusters, clusterCount, devPoints, pointCount);
//	
//	RETURN_IF_NOT_TRUE(copy_from_GPU_memory(points, devPoints, pointCount, sizeof(point_t)));
//	
//	RETURN_IF_NOT_TRUE(free_GPU_space(devClusters));
//	RETURN_IF_NOT_TRUE(free_GPU_space(devPoints));
//
//	return 1;
//}
//
////int cuda_calc_new_center(cluster_t *clusters, int clusterCount)
////{
////
////	return 1;
////}
//
//
//
//__global__ void calcSumOfHelpersFromPoints(clusterHelper_t *devHelpers, point_t *devPoints, int numPoints)
//{
//	int offset = blockIdx.x * blockDim.x + threadIdx.x;
//	if (offset < numPoints)
//	{
//		point_t *point = &(devPoints[offset]);
//		clusterHelper_t *helper = &(devHelpers[point->clusterNumber]);
//		
//		double *helperVectorArg;
//		double pointVectorArg;
//		int j;
//		for (j = 0; j < NUM_OF_ELEMENTS_IN_VECTOR; j++)
//		{
//			helperVectorArg = &(helper->sumPointLocation.axis[j]);
//			pointVectorArg = point->location.axis[j];
//			atomicAdd((float*)helperVectorArg, (float)pointVectorArg);
//		}
//	}
//}
//
//int cuda_calc_partialSums(clusterHelper_t *clusterHelpers, int numOfHelpers, point_t *points, int pointCount)
//{
//
//	setDevice();
//	int numOfBlocks = pointCount / NUM_OF_THREADS_PER_BLOCK;
//
//	point_t *devPoints;
//	clusterHelper_t *devClusterHelpers;
//
//	//allocate gpu space
//	RETURN_IF_NOT_TRUE(allocate_GPU_space((void**)&devPoints, pointCount, sizeof(point_t)));
//	RETURN_IF_NOT_TRUE(allocate_GPU_space((void**)&devClusterHelpers, numOfHelpers, sizeof(clusterHelper_t)));
//	
//	//copy points to gpu
//	RETURN_IF_NOT_TRUE(copy_to_GPU_memory(devPoints, points, pointCount, sizeof(point_t)));
//
//	calcSumOfHelpersFromPoints <<<numOfBlocks, NUM_OF_THREADS_PER_BLOCK>>>(devClusterHelpers, devPoints, pointCount);
//
//	RETURN_IF_NOT_TRUE(copy_from_GPU_memory(clusterHelpers, devClusterHelpers, numOfHelpers, sizeof(clusterHelper_t)));
//
//	RETURN_IF_NOT_TRUE(free_GPU_space(devClusterHelpers));
//	RETURN_IF_NOT_TRUE(free_GPU_space(devPoints));
//
//
//	return 1;
//}