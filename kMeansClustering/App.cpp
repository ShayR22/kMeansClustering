#include "mpi_master_slave.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
	if (argc != 3)
	{
		printf("this program require 3 arguments: \n");
		printf(" - first: the exe\n");
		printf(" - second: path to where input file is located\n");
		printf(" - third: full path + name of the file to write the ouput \n");
		fflush(stdout);
	}
	else
	{
		run(argc, argv);
	}
}
