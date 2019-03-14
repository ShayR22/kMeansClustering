
#ifndef OPENMP_GENERAL_HEADER
#define OPENMP_GENERAL_HEADER

void openMP_calc_pragma_range(int count, int *threadWorking, int *range, int *id, int *start, int *end);

void openMP_calc_new_deviation_based_on_time(double cudaStart, double cudaEnd, double openMpStart, double openMpEnd, double *cudaDevisionPoints);

#endif // !OPENMP_GENERAL_HEADER

