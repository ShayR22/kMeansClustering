#include "cudaHeader.h"
#include "mpi_master_slave_header.h"
#include "master_header.h"
#include "slave_header.h"
#include "rwFileHeader.h"
#include "cluster_header.h"
#include "point_header.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


MPI_Datatype mpi_cluster_type;
MPI_Datatype mpi_clusterHelper_type;
MPI_Datatype mpi_point_type;



static void process_init(int argc, char *argv[], int *id, int *numProcs);



static void create_all_dataTypes();


static void system_init(int id, int numProcs, cluster_t **clusters, clusterHelper_t **helpers, int *k,
						point_t **points, int *pointsCount, int *LIMIT, double *T, double *dt, double *QM );

static void allocate_clusterHelpers(int id, int numProcs, clusterHelper_t **helpers, int helpersCount, int numPoints);


static void kMeans(int id, int numProcs, cluster_t *clusters, clusterHelper_t *helpers, int k,
						   point_t* allPoints, point_t *pointsWorkLoad, int pointCount, int LIMIT);

static void distribute_points(int id, int numProcs, point_t *points, int pointCount, point_t **pointsWorkLoad, int *numOfPoints);

static void distribute_clusterHelpers(int id, clusterHelper_t *helpers, int helpersCount, int numProcs, int *numOfHelpersToWorkOn);

void run(int argc, char *argv[])
{
	int id, numProcs;
	process_init(argc, argv, &id, &numProcs);
	create_all_dataTypes();

	//TODO remove this
	double totalStart, totalEnd;
	totalStart = MPI_Wtime();
	

	//all data read from file
	cluster_t *clusters = NULL;
	point_t *points = NULL;
	int k, pointsCount, LIMIT;
	double T, dt, QM;

	clusterHelper_t *helpers = NULL;

	system_init(id, numProcs, &clusters, &helpers ,&k, &points, &pointsCount, &LIMIT, &T, &dt, &QM);

	point_t *startFromPoint;
	int numOfPointsToWorkOn;
	distribute_points(id, numProcs, points, pointsCount, &startFromPoint, &numOfPointsToWorkOn);

	int numOfHelpersToWorkOn;
	double time = 0;
	double q;

	double *diameterBuffer = (double*)calloc(k, sizeof(double));

	int testt = 0;
	int outerIterationForQM = 1;
	do
	{
		kMeans(id, numProcs, clusters, helpers, k, points, startFromPoint, numOfPointsToWorkOn, LIMIT);
		distribute_clusterHelpers(id, helpers, k, numProcs, &numOfHelpersToWorkOn);
		cluster_calculate_helpers_diameter(helpers, numOfHelpersToWorkOn, points, pointsCount);
	
		if (id == MASTER)
		{
			master_copy_diameters_from_its_own_helpers(clusters, helpers, numOfHelpersToWorkOn);
			master_recv_clusters_diameters(clusters, k, numProcs, numOfHelpersToWorkOn, diameterBuffer);
		}
		else
		{
			slave_send_clusterhelpers_diameters(helpers, numOfHelpersToWorkOn, diameterBuffer);
		}


		//TODO maybe add this inside if to the up if(id == MASTER)
		if (id == MASTER)
		{
			q = cluster_calculate_quality_measure(clusters, k);
			printf("master: q = %lf\n", q);
			fflush(stdout);
			if (q < QM || time > T)
			{
				outerIterationForQM = 0;
			}
			else
			{
				time += dt;
				outerIterationForQM = 1;
			}
		}

		MPI_Bcast(&outerIterationForQM, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

		if (outerIterationForQM)
		{
			cuda_increment_points(points, pointsCount, dt);
		}

	
	} while (outerIterationForQM);


	totalEnd = MPI_Wtime();
	
	if (id == MASTER)
	{
		if (time > T)
		{
			printf("Failed To Find Clusters with Quality Measure good enough\n");
			fflush(stdout);
		}

		printf("q = %lf required %lf\n", q, QM);
		fflush(stdout);

		printf("Total Time took is %lf\n", totalEnd - totalStart);
		fflush(stdout);

		printf("\nFinal Result:\n");
		fflush(stdout);
		cluster_print_clusters(clusters, k);
		file_write_kMeans(clusters, k, (time - dt), q);
	}

	MPI_Finalize();
}

//Initialize MPI process assigning id and number of process in the system
static void process_init(int argc, char *argv[], int *id, int *numProcs)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,id);
	MPI_Comm_size(MPI_COMM_WORLD, numProcs);
}

//Initialize all the MPI data types required for the project
static void create_all_dataTypes()
{
	mpi_cluster_type = cluster_create_cluster_mpi_struct();
	mpi_clusterHelper_type = cluster_create_clusterHelper_mpi_struct();
	mpi_point_type = point_create_mpi_struct();
}

/*
Initialize the system by:
	- Master reads from file the required data.
	- Master will broadcast the number of: clusters, points, kMeans Limit.
	- Every process will allocate the data required.
	- Master will broadcast all the points it read from the file.
	- Every process will Allocate clusterHelpers based on the data that was distributed.
*/
static void system_init(int id, int numProcs, cluster_t **clusters, clusterHelper_t **helpers, int *k, 
						point_t **points, int *pointsCount, int *LIMIT, double *T, double *dt, double *QM)
{
	//master read data from file
	if (id == MASTER)
	{
		//TODO maybe make this method return int and by that 
		file_read_file(clusters, k, points, pointsCount, T, dt, LIMIT, QM, FILE_NAME_READ);
	}
	
	//broadCast number of clusters
	MPI_Bcast(k, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	//broadCast number of points in the system
	MPI_Bcast(pointsCount, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	//broadCast LIMIT
	MPI_Bcast(LIMIT, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	//broadCast dt
	MPI_Bcast(dt, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	//TODO maybe broadCast more items from the read


	if (id != MASTER)
	{
		*clusters = (cluster_t*)calloc(*k, sizeof(cluster_t));
		*points = (point_t*)calloc(*pointsCount, sizeof(point_t));
	}
	
	//Master send all the points to all other machines
	MPI_Bcast(*points, *pointsCount, mpi_point_type, MASTER, MPI_COMM_WORLD);
	
	//TODO maybe pointsCount isnt right to allocate
	allocate_clusterHelpers(id, numProcs, helpers, *k, *pointsCount);
}


/*
 - Allocate clusterHelpers to the number of clusters there are.
 - Allocate for each helper an index arary as the number of points.
*/
static void allocate_clusterHelpers(int id, int numProcs, clusterHelper_t **helpers, int helpersCount, int numPoints)
{
	//master need to allocate more space to recevie results and fold them up
	if (id == MASTER)
	{
		*helpers = (clusterHelper_t*)calloc(helpersCount*numProcs, sizeof(clusterHelper_t));
		clusterHelper_allocate(*helpers, helpersCount*numProcs, numPoints);
	}
	else
	{
		*helpers = (clusterHelper_t*)calloc(helpersCount, sizeof(clusterHelper_t));
		clusterHelper_allocate(*helpers, helpersCount, numPoints);
	}

}


static void distribute_points(int id, int numProcs, point_t *points, int pointsCount, point_t **pointsWorkLoad, int *numOfPoints)
{
	if (id == MASTER)
		master_distribute_points(&mpi_point_type, points, pointsCount, numProcs, pointsWorkLoad, numOfPoints);
	else
		slaveRecvPoints(&mpi_point_type, points, pointsWorkLoad, numOfPoints);
}

static void kMeans(int id, int numProcs, cluster_t *clusters, clusterHelper_t *helpers, int k,
							  point_t* allPoints, point_t *pointsWorkLoad, int pointCount, int LIMIT)
{
	int numOfIteration = 1;
	double start, end;
	start = MPI_Wtime();
	int innerIterationKmeansAlgorithm = 0;

	double resetStart, resetEnd, addClusterToPointStart, addClusterToPointEnd;
	double addPointToClusterStart, addPointToClusterEnd;
	double calcHelperStart, calcHelperEnd;
	double foldAndCenterStart, foldAndCenterEnd;

	double resetTotal =0, addClusterToPointTotal=0, addPointToClusterTotal=0, calcHelperTotal=0, foldAndCenterTotal=0;

	do
	{
		// master broadcast the clusters (maybe no need in init)
		MPI_Bcast(clusters, k, mpi_cluster_type, MASTER, MPI_COMM_WORLD);

		resetStart = MPI_Wtime();
		//reset all the clusterHelpers
		cluster_reset_clusterHelpers(helpers, k); 
		resetEnd = MPI_Wtime();
		resetTotal += (resetEnd - resetStart);


		addClusterToPointStart = MPI_Wtime();
		//to each point add its nearest cluster
		cluster_add_to_each_point_its_nearest_cluster(clusters, k, pointsWorkLoad, pointCount);
		addClusterToPointEnd = MPI_Wtime();
		addClusterToPointTotal += (addClusterToPointEnd - addClusterToPointStart);

		
		addPointToClusterStart = MPI_Wtime();
		//go over helpers array and fill each of the helpers array with the correct ids
		cluster_add_to_each_clusterHelper_its_points(helpers, k, pointsWorkLoad, pointCount);
		addPointToClusterEnd = MPI_Wtime();
		addPointToClusterTotal += (addPointToClusterEnd - addPointToClusterStart);


		calcHelperStart = MPI_Wtime();
		//calculate the sum of the of vectors of the clusters
		cluster_calc_clusterHelper_points_location_sum(clusters, helpers, k, allPoints); 
		calcHelperEnd = MPI_Wtime();
		calcHelperTotal += (calcHelperEnd - calcHelperStart);

		foldAndCenterStart = MPI_Wtime();
		if (id == MASTER)
		{
			master_recv_clusterHelpers(helpers, k, numProcs, pointCount, &mpi_clusterHelper_type);
			master_fold_clusterHelpers(helpers, k, numProcs);
			innerIterationKmeansAlgorithm = cluster_is_center_changes(clusters, helpers, k);
		}
		else
		{
			slaveSendHelpers(&mpi_clusterHelper_type, helpers, k); 
		}
		foldAndCenterEnd = MPI_Wtime();
		foldAndCenterTotal += (foldAndCenterEnd - foldAndCenterStart);

;
		MPI_Bcast(&innerIterationKmeansAlgorithm, 1, MPI_INT, MASTER, MPI_COMM_WORLD); // broadcast

		numOfIteration++;
	} while (innerIterationKmeansAlgorithm && numOfIteration <= LIMIT);
	
	double totalTime = resetTotal + addPointToClusterTotal + addClusterToPointTotal + calcHelperTotal + foldAndCenterTotal;
	resetTotal /= totalTime;
	addPointToClusterTotal /= totalTime;
	addClusterToPointTotal /= totalTime;
	calcHelperTotal /= totalTime;
	foldAndCenterTotal /= totalTime;
	

	end = MPI_Wtime();
	if (id == MASTER)
	{
		//printf("\nresetHelpers = %lf%c, addClusterToPoint = %lf%c\n", resetTotal, '%', addClusterToPointTotal, '%');
		//printf("addPointToCluster = %lf%c, calcHelpers  = %lf%c, foldAndCenters = %lf%c\n\n", addPointToClusterTotal, '%', calcHelperTotal, '%', foldAndCenterTotal, '%');

		//printf("\nresetHelpers = %lf seconds, addClusterToPoint = %lf seconds\n", resetTotal*totalTime, addClusterToPointTotal*totalTime);
		//printf("addPointToCluster = %lf seconds, calcHelpers  = %lf seconds, foldAndCenters = %lf seconds\n\n", addPointToClusterTotal*totalTime, calcHelperTotal*totalTime, foldAndCenterTotal*totalTime);
	
		
		printf("time took iteration is: %lf with %d iterations\n\n", end - start, numOfIteration - 1);
		fflush(stdout);
	}
}

static void distribute_clusterHelpers(int id ,clusterHelper_t *helpers, int helpersCount, int numProcs, int *numOfHelpersToWorkOn)
{
	if (id == MASTER)
		master_send_clusterHelpers(helpers, helpersCount, numProcs, numOfHelpersToWorkOn, &mpi_clusterHelper_type);
	else
		slave_recv_clusterHelpers(&mpi_clusterHelper_type, helpers, numOfHelpersToWorkOn);

}































