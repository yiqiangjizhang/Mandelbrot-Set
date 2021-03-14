/* Mandelbrot set

Date: 27/02/2021
Author/s: Group 1
Subject: High Performance Computing for Aerospace Engineering
Professor: Manel Soria & Arnau Miro


// Problem statement
-------------------------------------------------------------------------
Plot Mandelbrot set using MPI
-------------------------------------------------------------------------

*/

// Libraries
#include <stdio.h>	// Standard Input and Output Library
#include <stdlib.h> // Library for defining functions to perform general functions (malloc())
#include <math.h>	// Math library
#include "mpi.h"	// MPI library

// Macro to access field u, v and w and GLOB
#define U(x, y) *(u + (x - map->sx + map->hs) + (y - map->sy + map->hs) * (map->ex - map->sx + 1 + 2 * map->hs))
#define V(x, y) *(v + (x - map->sx + map->hs) + (y - map->sy + map->hs) * (map->ex - map->sx + 1 + 2 * map->hs))
#define W(x, y) *(w + (x - map->sx + map->hs) + (y - map->sy + map->hs) * (map->ex - map->sx + 1 + 2 * map->hs))
#define GLOB(x, y) *(glob + (x - map->sx) + (y - map->sy) * (map->ex - map->sx + 1))

// This structure contains all the data to access a distributed 2d array
typedef struct Maps
{
	int sx, ex, sy, ey;
	int hs;
	//.... maybe we will need more data
} MAP;

// Prototypes
void worksplit(int *mystart, int *myend, int proc, int nproc, int start, int end);		// Worksplit action
void checkr(int r, char *txt);															// Check state
int proc();																				// Rank of the actual processor
int nproc();																			// Size of total number processors
void createMap(int NPX, int NPY, int gsx, int gex, int gsy, int gey, int hs, MAP *map); // Each processor create its own local map
void createglobalMap(int NPX, int NPY, int gsx, int gex, int gsy, int gey);				// Global map
double *allocGlobalMap(int gsx, int gex, int gsy, int gey);								// Allocate memory for the global field
void printMap(MAP *map);																// Print actual processor rank, mystart and myend
double *allocField(MAP *map);															// Allocate memory for local map
void fillField(double *u, MAP *m);														// Fill local map with values
void printField(double *u, MAP *map);													// Print local map values
void addField(double *u, double *v, double *w, MAP *map);								// Add map values
void fillGlobalField(double *u, MAP *map, double *glob);								// Fill global field with each processor's values

// Main
int main(int argc, char **argv)
{
	// Declare variables
	int NPX = 2; // Number of processors in X axis
	int NPY = 1; // Number of processors in Y axis
	int gsx = 1; // Global X start index
	int gex = 4; // Global X end index
	int gsy = 1; // Global Y start index
	int gey = 1; // Global Y end index
	int hs = 2;	 // Halo size

	// MPI error return value
	int r;

	double *u, *v, *w, *glob;
	MAP map_;
	MAP *map = &map_;

	// Initate MPI
	r = MPI_Init(&argc, &argv);
	checkr(r, "init");

	// Create each processors map
	createMap(NPX, NPY,			  // number of processors in each direction
			  gsx, gex, gsy, gey, // GLOBAL limits
			  hs,				  // Halo size
			  map);

	// Print actual processor rank, mystart and myend
	printMap(map);

	// Allocate memory for local map
	u = allocField(map); // Local vector (pointer)

	// v = allocField(map);
	// w = allocField(map);

	// Allocate memory for global map
	glob = allocGlobalMap(gsx, gex, gsy, gey); // Global vector (pointer)

	printf("memory allocated!\n");

	// Fill map with values
	fillField(u, map);
	printField(u, map);

	printf("field printed!\n");

	// fillField(v, map);
	// printField(v, map);
	// fillField(w, map);
	// printField(w, map);

	// Fill global field with each processor's values
	fillGlobalField(u, map, glob);

	printf("global map field printed!\n");

	// addField(u, v, w, map);

	// printField(u, map);
	// printf("field printed!\n");

	free(u);
	// free(v);
	// free(w);

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

// Action to check the return value of every MPI call
void checkr(int r, char *txt)
{
	if (r != MPI_SUCCESS)
	{
		fprintf(stderr, "Error: %s\n", txt);
		exit(-1);
	}
}

// Rank of the actual processor
int proc()
{
	int a, b;
	a = MPI_Comm_rank(MPI_COMM_WORLD, &b);
	checkr(a, "proc");
	return (b);
}

// Size of total number of processors
int nproc()
{
	int a, b;
	a = MPI_Comm_size(MPI_COMM_WORLD, &b);
	checkr(a, "nproc");
	return (b);
}

// each processor creates its map, filling in its values for sx,ex,sy,ey using worksplit
// a worksplit for each direction has to be done
// map->sx = ... ;
// map->sy = ... ;
void createMap(int NPX, int NPY,				   // number of processors in each direction
			   int gsx, int gex, int gsy, int gey, // GLOBAL limits
			   int hs,							   // Halo size
			   MAP *map)
{

	map->hs = hs;
	map->sy = gsy;
	map->ey = gey;
	worksplit(&(map->sx), &(map->ex), proc(), nproc(), gsx, gex);
}

// Allocate memory for the global field
double *allocGlobalMap(int gsx, int gex, int gsy, int gey) // GLOBAL limits
{
	// allocs memory field
	// calculate the size of the field
	double *aux;
	// printf("debug: %d %d %d %d %d\n", map->sx, map->ex, map->sy, map->ey, map->hs);
	aux = (double *)malloc(sizeof(double) * (gex - gsx + 1) * (gey - gsy + 1));

	//if memory cannot be allocated
	if (aux == NULL)
	{
		printf("Error! Memory not allocated.\n");
		exit(-1);
	}
	return aux;
}

// Print actual processor rank, mystart and myend
void printMap(MAP *map)
{
	printf("So, as I'm processor %d, I start with x%d and end with x%d\n", proc(), map->sx, map->ex);
}

// Allocate memory for local map
double *allocField(MAP *map)
{
	// allocs memory field
	// calculate the size of the field
	double *aux;
	// printf("debug: %d %d %d %d %d\n", map->sx, map->ex, map->sy, map->ey, map->hs);
	aux = (double *)malloc(sizeof(double) * (map->ex - map->sx + 1 + 2 * map->hs) * (map->ey - map->sy + 1 + 2 * map->hs));

	//if memory cannot be allocat
	if (aux == NULL)
	{
		printf("Error! Memory not allocated.\n");
		exit(-1);
	}
	return aux;
}

// Fill local map with values
void fillField(double *u, MAP *map)
{
	// each processor fill its part of the field
	for (int j = map->sy; j <= (map->ey); j++)
	{
		for (int i = map->sx; i <= (map->ex); i++)
		{
			U(i, j) = (i + j / 10.0);
		}
	}
}

// Fill global field with each processor's values
void fillGlobalField(double *u, MAP *map, double *glob)
{
	// each processor fill its part of the global field
	int r;				// for error checking
	int t1 = 0, t2 = 1; // tag, allows to identify messages
	MPI_Status st;		// check tag

	if (proc() == 0)
	{
		printf("global values\n");
		for (int j = map->sy; j <= (map->ey); j++)
		{
			for (int i = map->sx; i <= (map->ex); i++)
			{
				GLOB(i, j) = U(i, j);
				printf("%lf ", GLOB(i, j));
			}
			printf("\n");
		}

		for (int i = 1; i <= (nproc() - 1); i++)
		{
			double vr1[(map->ex - map->sx + 1) * (map->ey - map->sy + 1)];
			int vr2[4];
			r = MPI_Recv(vr1, (map->ex - map->sx + 1) * (map->ey - map->sy + 1), MPI_DOUBLE, i, t1, MPI_COMM_WORLD, &st);
			checkr(r, "receive1");
			// r = MPI_Recv(vr2, 4, MPI_INT, i, t2, MPI_COMM_WORLD, &st);
			// checkr(r, "receive2");

			printf("he rebut %lf %lf \n", *(vr1 + 1), *(vr1 + 2));

			// for (int j2 = vr2[2]; j2 <= vr2[3]; j2++)
			// {
			// 	for (int i2 = vr2[0]; i2 <= vr2[1]; i2++)
			// 	{
			// 		GLOB(i, j) = U(i, j);
			// 		printf("%lf ", GLOB(i, j));
			// 	}
			// 	printf("\n");
			// }
		}
	}
	else
	{
		printf("no soc 0");
		// double *vs1; // local values array
		// vs1 = (double *)malloc(sizeof(double) * (map->ex - map->sx + 1) * (map->ey - map->sy + 1));

		// //if memory cannot be allocated
		// if (vs1 == NULL)
		// {
		// 	printf("Error! Memory not allocated.\n");
		// 	exit(-1);
		// }

		// printf("envio aquests valors\n");

		// for (int j = map->sy; j <= (map->ey); j++)
		// {
		// 	for (int i = map->sx; i <= (map->ex); i++)
		// 	{
		// 		*(vs1 + i) = U(i, j);
		// 	}
		// }
		// printf("%lf %lf \n", *(vs1 ), *(vs1 + 1));

		// int vs2[4]; // local index array
		// vs2[0] = map->sx;
		// vs2[1] = map->ex;
		// vs2[2] = map->sy;
		// vs2[3] = map->ey;

		r = MPI_Ssend(&u, (map->ex - map->sx + 1) * (map->ey - map->sy + 1), MPI_DOUBLE, 0 /*destination*/, t1, MPI_COMM_WORLD);
		checkr(r, "send1");
		// r = MPI_Ssend(vs2, 4, MPI_INT, 0 /*destination*/, t2, MPI_COMM_WORLD);
		// checkr(r, "send2");
		printf("soc 1 i he acabat. \n");
	}
}

// Print local map values
void printField(double *u, MAP *map)
{
	// each processor prints its part of the field
	for (int j = map->sy; j <= (map->ey); j++)
	{
		for (int i = map->sx; i <= (map->ex); i++)
		{
			printf("%lf ", U(i, j));
		}
		printf("\n");
	}
}

// u=v+w
void addField(double *u, double *v, double *w, MAP *map)
{
	// each processor sweeps its part of the field
	for (int j = map->sy; j <= (map->ey); j++)
	{
		for (int i = map->sx; i <= (map->ex); i++)
		{
			U(i, j) = V(i, j) + W(i, j);
		}
		printf("\n");
	}
}
