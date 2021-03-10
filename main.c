// Libraries
#include <stdio.h>	// Standard Input and Output Library
#include <stdlib.h> // Library for defining functions to perform general functions (malloc())
#include <math.h>	// Math library
#include "mpi.h"	// MPI library

// Prototypes
void worksplit(int *mystart, int *myend, int proc, int nproc, int start, int end);
void checkr(int r, char *txt);

// This structure contains all the data to access a distributed 2d array
typedef struct Maps
{
	int sx, ex, sy, ey;
	int hs;
	//.... maybe we will need more data
} MAP;

// Macro to access field u
// define U(i,j) *( u + ... (map->sx)...)
// define V(i,j)
void createMap(int NPX, int NPY,				   // number of processors in each direction
			   int gsx, int gex, int gsy, int gey, // GLOBAL limits
			   int hs,							   // Halo size
			   MAP *map)
{
	int proc; // From 0 to P-1
	int r;	  // for error checking

	r = MPI_Comm_rank(MPI_COMM_WORLD, &proc);
	checkr(r, "rank");

	// r = MPI_Comm_size(MPI_COMM_WORLD,&nproc); checkr(r,"size");
	// each processor creates its map, filling in its values for sx,ex,sy,ey using worksplit
	// a worksplit for each direction has to be done
	// map->sx = ... ;
	// map->sy = ... ;

	worksplit(&map->sx, &map->ex, proc, NPX, gsx, gex);
	printf("So, as I'm processor %d, I start with x%d and end with x%d\n", proc, map->sx, map->ex);
	worksplit(&map->sy, &map->ey, proc, NPY, gsy, gey);
	printf("So, as I'm processor %d, I start with y%d and end with y%d\n", proc, map->sy, map->ey);
}

void printMap(MAP *map)
{ // prints my memory map
}
double *allocField(MAP *map)
{
	// allocs memory field
	// calculate the size of the field
}
void printField(double *u, MAP *map)
{
	// each processor prints its part of the field
}
// u=v+w
void addField(double *u, double *v, double *w, MAP *map)
{
	// each processor sweeps its part of the field
	// int i,j;
	// for (j=map->sy; j<=map->ey; j++)
	// for (i=map->sx; i<=map->ex; i++)
	// U(i,j)=V(i,j)+W(i,j);
}

int main(int argc, char **argv)
{
	// Declare variables
	int NPX = 1;
	int NPY = 4;
	int gsx = 1;
	int gex = 10;
	int gsy = 1;
	int gey = 1;
	int hs = 2;

	int r; // for error checking

	// double *u, *v, *w;
	MAP map_;
	MAP *map = &map_;

	r = MPI_Init(&argc, &argv);
	checkr(r, "init");

	createMap(NPX, NPY,			  // number of processors in each direction
			  gsx, gex, gsy, gey, // GLOBAL limits
			  hs,				  // Halo size
			  map);

	// createMap(1, 12, 1, 12, 2, map);
	// u = allocField(map);
	// v = allocField(map);
	// w = allocField(map);
	// printField(u, map);
	// free(u);
	// exit(0);
	MPI_Finalize();

	exit(0);
}

// Worksplit function

void worksplit(int *mystart, int *myend, int proc, int nproc, int start, int end)
{
	// Number of tasks
	int ntask = end - start + 1;
	// printf("Number of operations %d\n", ntask);

	// Number of tasks per processor
	int tpProcessor = ntask / nproc;
	// printf("tpProcessor %d\n", tpProcessor);

	// Tasks left
	int remainder_tpProcessor = ntask % nproc;
	// printf("Remainder of tpProcessor %d\n", remainder_tpProcessor);

	// Algorithm for worksplit
	if (remainder_tpProcessor == 0)
	{
		*mystart = proc * tpProcessor + start;
		*myend = *mystart + tpProcessor - 1;
	}
	else
	{
		if (proc < remainder_tpProcessor)
		{
			*mystart = proc * (tpProcessor + 1) + start;
			*myend = *mystart + tpProcessor;
		}
		else
		{
			*mystart = remainder_tpProcessor * (tpProcessor + 1) + (proc - remainder_tpProcessor) * tpProcessor + start;
			*myend = *mystart + tpProcessor - 1;
		}
	}
}

// we use this function to check the return value of every MPI call
void checkr(int r, char *txt)
{
	if (r != MPI_SUCCESS)
	{
		fprintf(stderr, "Error: %s\n", txt);
		exit(-1);
	}
}
