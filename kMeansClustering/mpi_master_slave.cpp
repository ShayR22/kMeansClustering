#include "cuda.h"
#include "mpi_master_slave.h"
#include "master.h"
#include "slave.h"
#include "rwFile.h"
#include "cluster.h"
#include "point.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

MPI_Datatype mpi_cluster_type;
MPI_Datatype mpi_clusterHelper_type;
MPI_Datatype mpi_point_type;

static int k;
static int pointsCount;
static int id;
static int numProcs;
static int LIMIT;
static double T;
static double dt;
static double QM;

//init methods
static void process_init(int argc, char *argv[]);
static void create_all_dataTypes();
static void allocate_across_everyone(cluster_t **clusters, clusterHelper_t **helpers, point_t **points);
static void allocate_clusterHelpers(clusterHelper_t **helpers);
static void distribute_points(point_t *points, int pointCount, point_t **pointsWorkLoad, int *numOfHelpersToWorkOn);

//calculation methods
static void find_clusters(cluster_t *clusters, clusterHelper_t *helpers, point_t* points,
						  point_t *startFromPoint, int numOfPointsToWorkOn, double *qResult, double *timeResult);
static void kMeans(cluster_t *clusters, clusterHelper_t *helpers, point_t* allPoints, point_t *pointsWorkLoad, int numOfPointsToWorkOn);
static void distribute_clusterHelpers(clusterHelper_t *helpers, int *numOfHelpersToWorkOn);
static void gather_diameters(cluster_t *clusters, clusterHelper_t *helpers, int numOfHelpersToWorkOn, double *diameterBuffer);
static int is_keep_iterating(double time ,double q);

//ouput methods
static void write_output_file(cluster_t *clusters, double time, double q, char *fileName);
static void print_to_console_results(cluster_t *clusters, double q, double timeTook);


void run(int argc, char *argv[])
{
	process_init(argc, argv);

	double totalStart, totalEnd;
	totalStart = MPI_Wtime();

	create_all_dataTypes();

	//structs for finding clusters
	cluster_t *clusters = NULL;
	clusterHelper_t *helpers = NULL;
	point_t *points = NULL;

	//master read data from file
	if (id == MASTER)
	{
		file_read_file(&clusters, &k, &points, &pointsCount, &T, &dt, &LIMIT, &QM, argv[1]);
	}
	allocate_across_everyone(&clusters, &helpers, &points);

	//each process/machine will work from a point in the points array for a certain range 
	point_t *startFromPoint;
	int numOfPointsToWorkOn;
	distribute_points(points, pointsCount, &startFromPoint, &numOfPointsToWorkOn);

	double time, q;
	find_clusters(clusters, helpers, points, startFromPoint, numOfPointsToWorkOn, &q, &time);
	totalEnd = MPI_Wtime();
	
	if (id == MASTER)
	{
		int i;
		for (i = 0; i < k; i++)
		{
			printf("cluster number %d has %d points and diameter of %lf \n", i, helpers[i].numOfPoints, clusters[i].diameter);
			fflush(stdout);
		}

		write_output_file(clusters, time, q, argv[2]);
		double timeTook = totalEnd - totalStart;
		print_to_console_results(clusters, q, timeTook);
	}
	//free
	free(clusters);
	free(points);
	cluster_free_clusterHelpers_array(helpers, k);

	MPI_Finalize();
}

//Initialize MPI process assigning id and number of process in the system
static void process_init(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	omp_set_nested(true);
	omp_set_dynamic(0);
}

//Initialize all the MPI data types required for the project
static void create_all_dataTypes()
{
	mpi_cluster_type = cluster_create_cluster_mpi_struct();
	mpi_point_type = point_create_mpi_struct();
}

static void allocate_across_everyone(cluster_t **clusters, clusterHelper_t **helpers, point_t **points)
{
	//broadCast number of clusters
	MPI_Bcast(&k, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	//broadCast number of points in the system
	MPI_Bcast(&pointsCount, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	//broadCast LIMIT
	MPI_Bcast(&LIMIT, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	//broadCast dt
	MPI_Bcast(&dt, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	
	if (id != MASTER)
	{
		*clusters = (cluster_t*)calloc(k, sizeof(cluster_t));
		*points = (point_t*)calloc(pointsCount, sizeof(point_t));
	}
	//Master send all the points to all other machines
	MPI_Bcast(*points, pointsCount, mpi_point_type, MASTER, MPI_COMM_WORLD);
	allocate_clusterHelpers(helpers);
}

static void allocate_clusterHelpers(clusterHelper_t **helpers)
{
	/*
	NOTE: all helpers will allocate thier IDs array at the size of pointsCount,
	 although the workload distribution is via points at the kMeans. after that there
	 is a distribution of clusterHelpers for parallel computations of the clusters
	 diameters, therefore at the worst case one cluster will have all the points in it.
	 thus the allocation would be at the size of pointsCount.
	*/

	//master need to allocate more space to recevie results and fold them up
	if (id == MASTER)
	{
		*helpers = (clusterHelper_t*)calloc(k*numProcs, sizeof(clusterHelper_t));
		clusterHelper_allocate(*helpers, k*numProcs, pointsCount);
	}
	else
	{
		*helpers = (clusterHelper_t*)calloc(k, sizeof(clusterHelper_t));
		clusterHelper_allocate(*helpers, k, pointsCount);
	}
}

static void distribute_points(point_t *points, int pointsCount, point_t **pointsWorkLoad, int *numOfHelpersToWorkOn)
{
	if (id == MASTER)
		master_distribute_points(&mpi_point_type, points, pointsCount, numProcs, pointsWorkLoad, numOfHelpersToWorkOn);
	else
		slaveRecvPoints(&mpi_point_type, points, pointsWorkLoad, numOfHelpersToWorkOn);
}

static void find_clusters(cluster_t *clusters, clusterHelper_t *helpers, point_t* points, 
						point_t *startFromPoint, int numOfPointsToWorkOn, double *q, double *time)
{
	double startIterationTime, endIterationTime;
	double *diameterBuffer = (double*)calloc(k, sizeof(double));
	int numOfHelpersToWorkOn;
	*time = 0;
	int outerIterationForQM;
	do
	{
		startIterationTime = MPI_Wtime();
		kMeans(clusters, helpers, points, startFromPoint, numOfPointsToWorkOn);


		distribute_clusterHelpers(helpers, &numOfHelpersToWorkOn);

		cluster_calculate_helpers_diameter(helpers, numOfHelpersToWorkOn, points, pointsCount);

		gather_diameters(clusters, helpers, numOfHelpersToWorkOn, diameterBuffer);


		if (id == MASTER)
		{
			*q = cluster_calculate_quality_measure(clusters, k);
			int iterationNumber = (int)(*time / dt);
			printf("master: iteration %d, q = %lf, timeFrame = %lf\n", iterationNumber, *q, *time);
			fflush(stdout);
			outerIterationForQM = is_keep_iterating(*time, *q);
			if (outerIterationForQM)
			{
				*time += dt;
			}
		}
		MPI_Bcast(&outerIterationForQM, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

		if (outerIterationForQM)
		{
			point_advance_points_by_dt(points, pointsCount, dt);
		}
		endIterationTime = MPI_Wtime();
		if (id == MASTER)
		{
			printf("entire frame iteration took %lf seconds\n\n", (endIterationTime - startIterationTime));
			fflush(stdout);
		}

	} while (outerIterationForQM);
	free(diameterBuffer);
}

static void kMeans(cluster_t *clusters, clusterHelper_t *helpers, point_t* allPoints, point_t *pointsWorkLoad, int numOfPointsToWorkOn)
{

	double start, end;
	start = MPI_Wtime();
	int numOfIteration = 1;
	int innerIterationKmeansAlgorithm = 0;
	double resetTotal =0, addClusterToPointTotal=0, addPointToClusterTotal=0, calcHelperTotal=0, foldAndCenterTotal=0;
	do
	{
		// master broadcast the clusters (maybe no need in init)
		MPI_Bcast(clusters, k, mpi_cluster_type, MASTER, MPI_COMM_WORLD);
		
		cluster_reset_clusterHelpers(helpers, k); 
		cluster_add_to_each_point_its_nearest_cluster(clusters, k, pointsWorkLoad, numOfPointsToWorkOn);
		cluster_add_to_each_clusterHelper_its_points(helpers, k, pointsWorkLoad, numOfPointsToWorkOn);
		cluster_calc_clusterHelper_points_location_sum(helpers, k, allPoints); 

		if (id == MASTER)
		{
			master_recv_clusterHelpers(helpers, k, numProcs);
			master_fold_clusterHelpers(helpers, k, numProcs);
			innerIterationKmeansAlgorithm = cluster_is_center_changes(clusters, helpers, k);	
		}
		else
		{
			slaveSendHelpers(helpers, k);
		}
		
		MPI_Bcast(&innerIterationKmeansAlgorithm, 1, MPI_INT, MASTER, MPI_COMM_WORLD); // broadcast

		numOfIteration++;
	} while (innerIterationKmeansAlgorithm && numOfIteration <= LIMIT);
		  
	end = MPI_Wtime();
	if (id == MASTER)
	{
		printf("time took kMeans is: %lf with %d iterations\n", (end - start), (numOfIteration - 1));
		fflush(stdout);
	}
}

static void distribute_clusterHelpers(clusterHelper_t *helpers, int *numOfHelpersToWorkOn)
{
	if (id == MASTER)
		master_send_clusterHelpers(helpers, k, numProcs, numOfHelpersToWorkOn);
	else
		slave_recv_clusterHelpers(helpers, numOfHelpersToWorkOn);
}

static void gather_diameters(cluster_t *clusters, clusterHelper_t *helpers, int numOfHelpersToWorkOn, double *diameterBuffer)
{
	if (id == MASTER)
	{
		master_copy_diameters_from_its_own_helpers(clusters, helpers, numOfHelpersToWorkOn);
		master_recv_clusters_diameters(clusters, k, numProcs, numOfHelpersToWorkOn, diameterBuffer);
	}
	else
	{
		slave_send_clusterhelpers_diameters(helpers, numOfHelpersToWorkOn, diameterBuffer);
	}
}

static int is_keep_iterating(double time, double q)
{
	if (q < QM || time + dt > T)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

static void write_output_file(cluster_t *clusters, double time, double q, char *fileName)
{
	if (time > T)
	{
		printf("Failed To Find Clusters with Quality Measure good enough\n");
		fflush(stdout);
		file_write_kMeans_fail(fileName, clusters, k, time, q);
	}
	else
	{
		file_write_kMeans_success(fileName, clusters, k, time, q);
	}
}

static void print_to_console_results(cluster_t *clusters, double q, double timeTook)
{
	printf("q = %lf required %lf\n", q, QM);
	printf("Total Time took is %lf\n", timeTook);
	printf("\nFinal Result:\n");
	fflush(stdout);
	cluster_print_clusters(clusters, k);
}
