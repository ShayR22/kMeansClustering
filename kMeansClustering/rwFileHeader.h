
#ifndef READ_FILE_HEADER
#include "cluster_header.h"
#include "point_header.h"

//TODO fix file names
#define FILE_NAME_READ "C:/Users/Shuy/Desktop/mevuzar project/mpiProjectFinal/mpiProjectFinal/input.txt"
#define FILE_NAME_WRITE "C:/Users/Shuy/Desktop/mevuzar project/mpiProjectFinal/mpiProjectFinal/output.txt"

#define FILE_NAME_FORMAT_T "First Occurrence t = "
#define FILE_NAME_FORMAT_Q "with q = "
#define FILE_NAME_FORMAT_CLUSTER "Centers of the clusters:"

void file_read_file(cluster_t  **emptyClusters, int *k, point_t **emptyPoints, int *pointsCount,
		   	  double *T, double *dt, int *LIMIT, double *QM, char* fileName);


void file_write_kMeans(cluster_t *clusters, int k, double time, double q);


#endif // !READ_FILE_HEADER
