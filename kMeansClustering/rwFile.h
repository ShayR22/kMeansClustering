
#ifndef READ_FILE_HEADER
#include "cluster.h"
#include "point.h"

#define FILE_NAME_FORMAT_T_SUCCESS "First Occurrence t = "
#define FILE_NAME_FORMAT_T_FAIL "Stoped at t = "

#define FILE_NAME_FORMAT_Q_SUCCESS "with q = "
#define FILE_NAME_FOMRAT_Q_FAIL "with last q = "

#define FILE_NAME_FORMAT_CLUSTER "Centers of the clusters:"

/*
read specific input file which have in this order: 
number of points, number of clusters, time reference, dt (time frame skip), limit (represent the limit of inner
iteration for stabilizing a cluster system), quality meassure.

Parameters:
- emptyClusters: pointer for pointer of clusters array to be initalize in method base on file data.
- k: pointer for int which reprsent the length of emptyClusters array.
- emptyPoints: pointer for pointer of points array to be initalize in method base on file data.
- pointsCount: pointer for int which reprsent the length of emptyPoints array.
- T: pointer for double which will represent the size of the time refernce to search the required quality measure.
- dt: pointer for double which will represent the time skip between different time frames.
- QM: pointer for double which will represent the desirable quality measure. 
- fileName: full path to the input file location.

Returns: void
*/
void file_read_file(cluster_t  **emptyClusters, int *k, point_t **emptyPoints, int *pointsCount,
		   	  double *T, double *dt, int *LIMIT, double *QM, char* fileName);


/*
write clusters locaiton, quality measure that had been found and the time in which the quality measure
has been found.

Parameters:
- fileName: full path to the output file location.
- clusters: pointer of clusters array.
- k: length of clusters array.
- t: time reference in which quality measure has been found or not found
- q: quality measure that have been achieved

Returns: void
*/
void file_write_kMeans_success(char *fileName, cluster_t *clusters, int k, double time, double q);

/*
write clusters locaiton, last quality meassure that was achieved before time run out

Parameters:
- fileName: full path to the output file location.
- clusters: pointer of clusters array.
- k: length of clusters array.
- t: time reference in which quality measure has been found or not found
- q: quality measure that have been achieved

Returns: void
*/
void file_write_kMeans_fail(char *fileName, cluster_t *clusters, int k, double time, double q);

#endif // !READ_FILE_HEADER