#define _CRT_SECURE_NO_WARNINGS

#include "rwFile.h"
#include "cluster.h"
#include "point.h"
#include "vector.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void readPoints(FILE *file, point_t *points, int numPoints);
static void readPoint(FILE *file, point_t *p);

static void write_vector(FILE *openFile, vector_t *vec);
static void writeClusterCenters(FILE *openFile, cluster_t *clusters, int clustersCount);

void file_read_file(cluster_t  **emptyClusters, int *k, point_t **emptyPoints, int *pointsCount,
			  double *T, double *dt, int *LIMIT, double *QM, char* fileName)
{
	// The first line of the file contains   N    K    T   dT   LIMIT   QM.
	// Next lines are Initial Positions and Velocities of the points(xi, yi, zi, vxi, vyi, vzi)
	FILE *inputFile;
	inputFile = fopen(fileName, "r");
	if (!inputFile)
	{
		printf("file open failed\n");
		fflush(stdout);
	}
	fscanf(inputFile, "%d %d %lf %lf %d %lf", pointsCount, k, T, dt, LIMIT, QM);

	*emptyPoints = (point_t*)calloc(*pointsCount, sizeof(point_t));
	readPoints(inputFile, *emptyPoints, *pointsCount);

	*emptyClusters = (cluster_t*)calloc(*k, sizeof(cluster_t));
	cluster_t *clusters = *emptyClusters;
	int i;
	for(i = 0; i < *k; i++)
	{
		cluster_init(&(clusters[i]), &((*emptyPoints)[i].location));
	}

	fclose(inputFile);
}

static void readPoints(FILE *file, point_t *points, int numPoints)
{
	int i;
	for (i = 0; i < numPoints; i++)
	{
		readPoint(file, &(points[i]));
		points[i].id = i;
	}
}

static void readPoint(FILE *file, point_t *p)
{
	vector_t *l = &(p->location);
	vector_t *s = &(p->speed);

	fscanf(file, "%lf %lf %lf", &(l->axis[0]), &(l->axis[1]), &(l->axis[2]));
	fscanf(file, "%lf %lf %lf", &(s->axis[0]), &(s->axis[1]), &(s->axis[2]));
}

void file_write_kMeans_success(char *fileName, cluster_t *clusters, int k, double time, double q)
{
	FILE *output;
	output = fopen(fileName, "w");
	if (!output)
	{
		printf("error on opening input\n");
	}

	fprintf(output, "%s %lf %s %lf \n", FILE_NAME_FORMAT_T_SUCCESS, time, FILE_NAME_FORMAT_Q_SUCCESS, q);
	fprintf(output, "%s \n", FILE_NAME_FORMAT_CLUSTER);
	
	writeClusterCenters(output, clusters, k);
	
	fclose(output);
}

void file_write_kMeans_fail(char *fileName, cluster_t *clusters, int k, double time, double q)
{
	FILE *output;
	output = fopen(fileName, "w");
	if (!output)
	{
		printf("error on opening input\n");
	}

	fprintf(output, "%s %lf %s %lf \n", FILE_NAME_FORMAT_T_FAIL, time, FILE_NAME_FOMRAT_Q_FAIL, q);
	fprintf(output, "%s \n", FILE_NAME_FORMAT_CLUSTER);

	writeClusterCenters(output, clusters, k);
	
	fclose(output);
}

static void writeClusterCenters(FILE *openFile, cluster_t *clusters, int clustersCount)
{
	vector_t *currentCenter;
	int i;
	for (i = 0; i < clustersCount; i++)
	{
		currentCenter = &(clusters[i].center);
		write_vector(openFile, currentCenter);
	}
}

static void write_vector(FILE *openFile, vector_t *vec)
{
	int j;
	for (j = 0; j < NUM_OF_ELEMENTS_IN_VECTOR; j++)
	{
		fprintf(openFile, "%lf ", vec->axis[j]);
	}
	fprintf(openFile, "\n");
}