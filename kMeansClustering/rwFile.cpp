#define _CRT_SECURE_NO_WARNINGS

#include "rwFileHeader.h"
#include "cluster_header.h"
#include "point_header.h"
#include "vector_header.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>


#define INIT_LINE_SIZE_BYTES (3*sizeof(int) + 3*sizeof(double))

static void readPoints(FILE *file, point_t *points, int numPoints);
static void readPoint(FILE *file, point_t *p);


void file_read_file(cluster_t  **emptyClusters, int *k, point_t **emptyPoints, int *pointsCount,
			  double *T, double *dt, int *LIMIT, double *QM, char* fileName)
{
	// The first line of the file contains   N    K    T   dT   LIMIT   QM.
	// Next lines are Initial Positions and Velocities of the points(xi, yi, zi, vxi, vyi, vzi)
	FILE *file;
	file = fopen(fileName, "r");
	if (!file)
	{
		printf("file open failed\n");
	}
	fscanf(file, "%d %d %lf %lf %d %lf", pointsCount, k, T, dt, LIMIT, QM);



	*emptyPoints = (point_t*)calloc(*pointsCount, sizeof(point_t));
	readPoints(file, *emptyPoints, *pointsCount);

	*emptyClusters = (cluster_t*)calloc(*k, sizeof(cluster_t));
	cluster_t *clusters = *emptyClusters;
	int i;
	for(i = 0; i < *k; i++)
	{
		cluster_init(&(clusters[i]), &((*emptyPoints)[i].location));
	}
}

static void readPoints(FILE *file, point_t *points, int numPoints)
{
	//regular read
	int i;
	for (i = 0; i < numPoints; i++)
	{
		readPoint(file, &(points[i]));
		points[i].id = i;
	}

	
	//TODO maybe fix parallel reading
	//point_t *pointsss = (point_t*)calloc(numPoints, sizeof(point_t));
	//#pragma omp parallel
	//{
	//	//compute start location and endlocation and read points 
	//	// start location is seek location in the file
	//	int threadsUsed = omp_get_num_threads();
	//	int numPointsToRead = numPoints / threadsUsed;
	//	int id = omp_get_thread_num();
	//	int start = id * numPointsToRead;
	//	int end = start + numPointsToRead;
	//	if (id == threadsUsed - 1) // last thread
	//	{
	//		end += (numPoints - end);
	//	}	
	//	printf("id is: %d, start is: %d, end is %d \n", id, start, end);
	//	fflush(stdout);

	//	int i;
	//	fseek(file, start + INIT_LINE_SIZE_BYTES, SEEK_SET);
	//	for (i = start; i < end; i++)
	//	{
	//		readPoint(file, &(pointsss[i]));
	//		pointsss[i].id = i;
	//	}
	//}
	//int emptyPointsCounter = 0;

	//for (i = 0; i < numPoints; i++)
	//{
	////	printPoint(&(pointsss[i]));
	//	if (pointsss[i].location.axis[0] == 0 && pointsss[i].location.axis[1] == 0 && pointsss[i].location.axis[2] == 0)
	//	{
	//		emptyPointsCounter++;
	//	}
	//}
	////check if the same points have been read
	//for (i = 0; i < numPoints; i++)
	//{
	//	if (!isVectorEqual(&(points[i].location), &(pointsss[i].location)))
	//	{
	//		printf("\nnot the same read\n");
	//	}
	//}
	//printf("same points\n");
}

static void readPoint(FILE *file, point_t *p)
{
	vector_t *l = &(p->location);
	vector_t *s = &(p->speed);

	fscanf(file, "%lf %lf %lf", &(l->axis[0]), &(l->axis[1]), &(l->axis[2]));
	fscanf(file, "%lf %lf %lf", &(s->axis[0]), &(s->axis[1]), &(s->axis[2]));
}

void file_write_kMeans(cluster_t *clusters, int k, double time, double q)
{
	FILE *output;
	output = fopen(FILE_NAME_WRITE, "w");
	if (!output)
	{
		printf("error on opening input\n");
	}

	fprintf(output, "%s %lf %s %lf \n", FILE_NAME_FORMAT_T, time, FILE_NAME_FORMAT_Q, q);
	fprintf(output, "%s \n", FILE_NAME_FORMAT_CLUSTER);
	
	vector_t *currentCenter;
	int i, j;
	for (i = 0; i < k; i++)
	{
		currentCenter = &(clusters[i].center);
		for (j = 0; j < NUM_OF_ELEMENTS_IN_VECTOR; j++)
		{
			fprintf(output, "%lf ", currentCenter->axis[j]);
		}
		fprintf(output, "\n");
	}
}
