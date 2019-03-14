#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda.h" // my header
#include <cuda.h> // official api header
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define CPU_TO_GPU cudaMemcpyHostToDevice
#define GPU_TO_CPU cudaMemcpyDeviceToHost
#define WRAP_SIZE 32
#define NUM_OF_THREADS_PER_BLOCK 256
#define NUM_OF_POINTS_TO_GO_OVER_IN_CALC_DIAMETER_PER_THREAD 256
#define FIND_MAX_RUN_PER_POINT 10

#define RETURN_IF_NOT_TRUE(truth)if(!truth) {return 0;}

static void copy_ids_from_helpers_to_array(clusterHelper_t *helpers, int helpersCount, int *IDs);
static int find_2_power_of(int numPoints);
static void find_max_distance_for_each_point_in_helper(clusterHelper_t *helper, point_t *devPoints, double *devDiameters, int *devHelpersIDs, int *offset);
static void zero_to_the_power_of_two_complement(clusterHelper_t *helper, double *devDiameters, int *goOver);
static void zero_devArray(double *devDiameters, int arraySize);

static double find_max_from_devDiameters(double *devDiameters, int goOver);

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

__global__ void calcPointNearestCluster(point_t *devPoints, int pointsCount, cluster_t *devClusters, int clustersCount)
{
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	if (index < pointsCount)
	{
		vector_t *pointLocation = &(devPoints[index].location);
		vector_t *currentClusterLocation;

		double distance = INT_MAX;
		double currentDistance;
		int belongToo;
		int i;
		for (i = 0; i < clustersCount; i++)
		{
			currentClusterLocation = &(devClusters[i].center);
			currentDistance = getDistance(pointLocation, currentClusterLocation);
			if(currentDistance < distance)
			{
				distance = currentDistance;
				belongToo = i;
			}
		}
		devPoints[index].clusterBelongTo = belongToo;
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

__global__ void calcDiaPerOneHelperIteration(point_t *devPoints, double *devDiameters, int *devHelpersIDs, int devHelperIDsOffset, int devChunkSize, int start, int jumpSize)
{
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (offset + start < devChunkSize)
	{
		int *IDS = &(devHelpersIDs[devHelperIDsOffset]);
		vector_t *v1 = &(devPoints[IDS[offset]].location); // extract point in the offset location
		vector_t *v2;
		double maxDiameterPointOffset = 0;
		double currentDist;

		int startLoopFrom = offset + start + 1;
		int endLoopIn = startLoopFrom + jumpSize;
		if (endLoopIn > devChunkSize)
		{
			endLoopIn = devChunkSize;
		}

		//iterate to the right of the chunkSize and get maxValue for given point
		int i;
		for (i = startLoopFrom; i < endLoopIn; i++)
		{
			v2 = &(devPoints[IDS[i]].location);

			currentDist = getDistance(v1, v2);
			if (currentDist > maxDiameterPointOffset)
			{
				maxDiameterPointOffset = currentDist;
			}
		}
		if (maxDiameterPointOffset > devDiameters[offset])
		{
			devDiameters[offset] = maxDiameterPointOffset;
		}
	}
}

__global__ void find_max_iteration_2_power_n(double *devDiameters, int halfCompares)
{
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	if (index < halfCompares)
	{
		if (devDiameters[index] < devDiameters[2 * halfCompares - index - 1])
		{
			devDiameters[index] = devDiameters[2 * halfCompares - index - 1];
		}
	}
}

__global__ void zeroLocation(double *devDiameters, int devDiametersCount, int goOver)
{
	int index = blockDim.x * blockIdx.x + threadIdx.x + devDiametersCount;
	if (index < goOver + devDiametersCount)
	{
		devDiameters[index] = 0;
	}
}

int cuda_calc_each_helper_diamater(point_t *points, int pointsCount, clusterHelper_t *helpers, int helpersCount)
{
		setDevice();
		point_t *devPoints;
		int *allHelpersIDs;
		int *devHelpersIDs;
		double *devDiameters;
		int devDiametersCompPower2 = find_2_power_of(pointsCount);
			
		//allocate diameters in GPU
		RETURN_IF_NOT_TRUE(allocate_GPU_space((void**)&devDiameters, devDiametersCompPower2, sizeof(double)));
		//allocate points in GPU
		RETURN_IF_NOT_TRUE(allocate_GPU_space((void**)&devPoints, pointsCount, sizeof(point_t)));
		//copy points to GPU
		RETURN_IF_NOT_TRUE(copy_to_GPU_memory(devPoints, points, pointsCount, sizeof(point_t)));
		//allocate CPU array for all helpers IDs
		RETURN_IF_NOT_TRUE((allHelpersIDs = (int*)malloc(sizeof(int) * pointsCount)) != NULL);
		//copy each helper's ids into one array
		copy_ids_from_helpers_to_array(helpers, helpersCount, allHelpersIDs);
		//allocate helperIDs in GPU
		RETURN_IF_NOT_TRUE(allocate_GPU_space((void**)&devHelpersIDs, pointsCount, sizeof(int)));
		//copy helpersIDs to GPU
		RETURN_IF_NOT_TRUE(copy_to_GPU_memory(devHelpersIDs, allHelpersIDs, pointsCount, sizeof(int)));
		
		int goOver;
		int offset = 0;
		int i;
		for (i = 0; i < helpersCount; i++)
		{
			zero_devArray(devDiameters,helpers[i].numOfPoints);
			find_max_distance_for_each_point_in_helper(&(helpers[i]), devPoints, devDiameters, devHelpersIDs, &offset);

			zero_to_the_power_of_two_complement(&(helpers[i]), devDiameters, &goOver);

			helpers[i].diameter = find_max_from_devDiameters(devDiameters, goOver);
		}

		//free CPU
		free(allHelpersIDs);
		//free CUDA 
		RETURN_IF_NOT_TRUE(free_GPU_space(devDiameters)); // free diameters
		RETURN_IF_NOT_TRUE(free_GPU_space(devPoints)); // free points
		RETURN_IF_NOT_TRUE(free_GPU_space(devHelpersIDs)); // free IDs array	
		return 1;
}

static void copy_ids_from_helpers_to_array(clusterHelper_t *helpers, int helpersCount, int *IDs)
{
	int offSetIDs = 0;
	int i;
	for (i = 0; i < helpersCount; i++)
	{
		memcpy(&(IDs[offSetIDs]), helpers[i].pointsIDs, helpers[i].numOfPoints * sizeof(int));
		offSetIDs += helpers[i].numOfPoints;
	}
}

static void find_max_distance_for_each_point_in_helper(clusterHelper_t *helper, point_t *devPoints, double *devDiameters, int *devHelpersIDs ,int *offset)
{
	int numBlocksCalcDiameters;
	int maxIteration = (helper->numOfPoints + FIND_MAX_RUN_PER_POINT - 1) / FIND_MAX_RUN_PER_POINT;
	int start = 0;
	int i;
	for (i = 0; i < maxIteration; i++)
	{
		numBlocksCalcDiameters = ((helper->numOfPoints + NUM_OF_THREADS_PER_BLOCK - 1 - start) / NUM_OF_THREADS_PER_BLOCK);
		calcDiaPerOneHelperIteration<<<numBlocksCalcDiameters, NUM_OF_THREADS_PER_BLOCK >>>(devPoints, devDiameters, devHelpersIDs, *offset, helper->numOfPoints, start, FIND_MAX_RUN_PER_POINT);
		start += FIND_MAX_RUN_PER_POINT;
	}
	*offset += helper->numOfPoints;

}

static void zero_to_the_power_of_two_complement(clusterHelper_t *helper, double *devDiameters, int *goOver)
{
	*goOver = find_2_power_of(helper->numOfPoints);
	int numPointsToZero = *goOver - helper->numOfPoints;
	if (numPointsToZero != 0)
	{
		int numBlocksForZero = (numPointsToZero + NUM_OF_THREADS_PER_BLOCK - 1) / NUM_OF_THREADS_PER_BLOCK;
		zeroLocation << <numBlocksForZero, NUM_OF_THREADS_PER_BLOCK >> > (devDiameters, helper->numOfPoints, *goOver); //O(n /(k * cudaThreads)) ~= O(1)
	}
}

static void zero_devArray(double *devDiameters,int arraySize)
{
	int numBlocksForZero = (arraySize + NUM_OF_THREADS_PER_BLOCK - 1) / NUM_OF_THREADS_PER_BLOCK;
	zeroLocation << <numBlocksForZero, NUM_OF_THREADS_PER_BLOCK >> > (devDiameters, 0, arraySize); //O(n /(k * cudaThreads)) ~= O(1)
}

static double find_max_from_devDiameters(double *devDiameters, int goOver)
{
	double result;
	int i;
	for (i = goOver / 2; i >= 1; i= i / 2)
	{
		int numBlocksForFindMax = (i + NUM_OF_THREADS_PER_BLOCK - 1) / NUM_OF_THREADS_PER_BLOCK;
		if (numBlocksForFindMax != 0)
		{
			find_max_iteration_2_power_n << <numBlocksForFindMax, NUM_OF_THREADS_PER_BLOCK >> > (devDiameters, i);
		}
		else
		{
			find_max_iteration_2_power_n << <1, NUM_OF_THREADS_PER_BLOCK >> > (devDiameters, i);
		}
	}

	copy_from_GPU_memory(&result, devDiameters, 1, sizeof(double));
	return result;
}

static int find_2_power_of(int numPoints)
{
	int x = 1;
	while (x < numPoints)
	{
		x *= 2;
	}
	return x;
}

static int setDevice()
{
	int device = 0;
	cudaSetDeviceFlags(cudaDeviceBlockingSync);
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
		printf(cudaGetErrorName(error));
		printf("\n");
		fflush(stdout);
		return 0;
	}
	return 1;
}
