// Libraries
#include <stdio.h>	// Standard Input and Output Library
#include <stdlib.h> // Library for defining functions to perform general functions (malloc())
#include <math.h>	// Math library
#include "mpi.h"	// MPI library

// #define U(x,y) *(u + (x-sx+hs)+(y-sy+hs)*(ex-sx+1+2*hs));
#define U(x, y) *(u + (x - map->sx + map->hs) + (y - map->sy + map->hs) * (map->ex - map->sx + 1 + 2 * map->hs))

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

	// each processor creates its map, filling in its values for sx,ex,sy,ey using worksplit
	// a worksplit for each direction has to be done
	// map->sx = ... ;
	// map->sy = ... ;
	map->hs = hs;

	worksplit(&map->sx, &map->ex, proc, NPX, gsx, gex);
	printf("So, as I'm processor %d, I start with x%d and end with x%d\n", proc, map->sx, map->ex);

	if (proc < NPY)
	{
		worksplit(&map->sy, &map->ey, proc, NPY, gsy, gey);
		printf("So, as I'm processor %d, I start with y%d and end with y%d\n", proc, map->sy, map->ey);
	}
}

void printMap(MAP *map)
{ // prints my memory map
}
double *allocField(MAP *map)
{
	// allocs memory field
	// calculate the size of the field
	double* aux;
	printf("debug: %d %d %d %d %d\n",map->sx,map->ex,map->sy,map->ey,map->hs);
  aux = (double*) malloc(sizeof(double)*(map->ex-map->sx+1+2*map->hs)*(map->ey-map->sy+1+2*map->hs));

	//if memory cannot be allocat
	if(aux == NULL)
  {
    printf("Error! Memory not allocated.\n");
    exit(-1);
  }
  return aux;
}
void printField(double *u, MAP *map)
{
	// each processor prints its part of the field
	for (int j = map->sy; j <= (map->ey); j++)
	{
		for (int i = map->sx; i <= (map->ex); i++)
		{
			U(i, j) = (i - map->sx + map->hs) + (j - map->sy + map->hs) * (map->ex - map->sx + 1 + 2 * map->hs);
			printf("%lf ", U(i, j));
		}
		printf("\n");
	}
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
	int NPY = 1;
	int gsx = 1;
	int gex = 4;
	int gsy = 1;
	int gey = 4;
	int hs = 2;

  int r;     // for error checking
	double u_,v_,w_;
	double* u = &u_;
	// double* v = &v_;
	// double* w = &w_;

	// double *u, *v, *w;
	MAP map_;
	MAP *map = &map_;

	r = MPI_Init(&argc, &argv);
	checkr(r, "init");

	createMap(NPX, NPY,			  // number of processors in each direction
			  gsx, gex, gsy, gey, // GLOBAL limits
			  hs,				          // Halo size
			  map);
	// createMap(1, 12, 1, 12, 2, map);

	u = allocField(map); // this is a pointer!
	printf("memory allocated!\n");
	// v = allocField(map);
	// w = allocField(map);
  printField(u, map);
	printf("field printed!\n");
	// printField(u, map);

	free(u);

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
