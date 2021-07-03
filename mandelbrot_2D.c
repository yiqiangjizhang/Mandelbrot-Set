/* Mandelbrot set
Date: 17/03/2021
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

// Macros
#define U(x, y) *(u + (x - map->sx + map->hs) + (y - map->sy + map->hs) * (map->ex - map->sx + 1 + 2 * map->hs))
#define GLOB(x, y) *(glob + (x - map->gsx + map->hs) + (y - map->gsy + map->hs) * (map->gex - map->gsx + 1 + 2 * map->hs))
// #define GLOB(x,y) *(glob + x + (y-1) * (map->gex-map->gsx+1))
#define SR(x, y) *(sr2 + (x - sr1[0] + map->hs) + (y - sr1[2] + map->hs) * (sr1[1] - sr1[0] + 1 + 2 * map->hs))

// This structure contains all the data to access a distributed 2d array
typedef struct Maps
{
	int sx, ex, sy, ey;		// Sizes
	int gsx, gex, gsy, gey; //
	int hs;					// Halo size
} MAP;

// MPI function prototypes
void checkr(int r, char *txt);															// Check state
int proc();																				// Rank of the actual processor
int nproc();																			// Size of total number processors
void worksplit(int *mystart, int *myend, int proc, int nproc, int start, int end);		// Divide tasks between processors
void createMap(int NPX, int NPY, int gsx, int gex, int gsy, int gey, int hs, MAP *map); // Each processor create its own local map
void printMap(MAP *map);																// Print actual processor rank, mystart and myend
double *allocField(MAP *map);															// Allocate memory for local map
double *allocGlobalField(MAP *map);														// Allocate memory for the global field
void fillField(double *u, MAP *m);														// Fill local map with values
void printField(double *u, MAP *map);													// Print local map values
void fillGlobalField(double *u, double *glob, MAP *map);								// Fill global map with values
void printGlobalField(double *glob, MAP *map);											// Print global field values

// Mandelbrot function prototypes
double Mandelbrot(int i, int j, int gex, int gey);

int main(int argc, char **argv)
{
	// Map variables
	int NPX = 2;  // Number of processors in X axis
	int NPY = 2;  // Number of processors in Y axis
	int gsx = 1;  // Global X start index
	int gex = 10; // Global X end index
	int gsy = 1;  // Global Y start index
	int gey = 10; // Global Y end index
	int hs = 0;	  // Halo size
	MAP map_;
	MAP *map = &map_;

	// MPI variables
	int r; // Error checking

	// Memory variables
	double *u, *glob;

	// Start MPI
	r = MPI_Init(&argc, &argv);
	checkr(r, "Initiate");

	// Create local map
	createMap(NPX, NPY, gsx, gex, gsy, gey, hs, map);

	// Print local processor rank, mystart and myend
	// printMap(map);

	// Allocate memory for local map
	u = allocField(map);
	// printf("\nLocal memory allocated successfully! \n\n");

	// Fill map with values
	fillField(u, map);
	// printf("Local map filled successfully! \n\n");

	// Print local map
	// printField(u, map);

	// Allocate memory for global map
	glob = allocGlobalField(map);
	// printf("\nGlobal memory allocated successfully! \n\n");

	// Fill global field with each processor's values
	fillGlobalField(u, glob, map);

	if (proc() == 0)
	{
		// Print global map
		printGlobalField(glob, map);
		// printf("\nGlobal map printed successfully! \n\n");
	}

	// Free allocated memory
	free(u);
	free(glob);

	// End MPI
	MPI_Finalize();

	// End main
	exit(0);
}
// EOF
/* -------------------------------------------------------------------------- */

// Use this function to check the return value of every MPI call

void checkr(int r, char *txt)
{
	if (r != MPI_SUCCESS)
	{
		fprintf(stderr, "Error: %s\n", txt);
		exit(-1);
	}
}

// Indicates the rank of the processor

int proc()
{
	int r, rank;
	r = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	checkr(r, "proc");
	return (rank);
}

// Indicates the number of processors being used

int nproc()
{
	int r, size;
	r = MPI_Comm_size(MPI_COMM_WORLD, &size);
	checkr(r, "nproc");
	return (size);
}

// Divide tasks between processors

void worksplit(int *mystart, int *myend, int proc, int nproc, int start, int end)
{
	// Number of tasks
	int ntask = end - start + 1;
	// Number of tasks per processor
	int interval = ntask / nproc;
	// Tasks left
	int remainder = ntask % nproc;

	if (ntask < nproc)
	{
		printf("Less tasks than processors\n");
		exit(-1);
	}
	if (remainder != 0)
	{
		if (proc < remainder)
		{
			*mystart = start + proc * (interval + 1);
			*myend = *mystart + interval;
		}
		else
		{
			*mystart = start + remainder * (interval + 1) + (proc - remainder) * interval;
			*myend = *mystart + interval - 1;
		}
	}
	else
	{
		*mystart = start + proc * interval;
		*myend = *mystart + (interval - 1);
	}
}

// Each processor creates its map of local values

void createMap(int NPX, int NPY, int gsx, int gex, int gsy, int gey, int hs, MAP *map)
{
	map->gsx = gsx;
	map->gex = gex;
	map->gsy = gsy;
	map->gey = gey;
	map->hs = hs;
	// Worksplit in x direction
	worksplit(&(map->sx), &(map->ex), proc() % NPX, NPX, map->gsx, map->gex);
	// Worksplit in y direction
	worksplit(&(map->sy), &(map->ey), proc() / NPX, NPY, map->gsy, map->gey);
}

// Displays the memory map on the screen

void printMap(MAP *map)
{
	printf("So, as I'm processor %d, I start with x%d and end with x%d\n", proc(), map->sx, map->ex);
	printf("So, as I'm processor %d, I start with y%d and end with y%d\n", proc(), map->sy, map->ey);
}

// Memory allocation of each field (local limits)

double *allocField(MAP *map)
{

	double *myField;
	// printf("debug: %d %d %d %d %d\n", map->sx, map->ex, map->sy, map->ey, map->hs);
	myField = (double *)malloc(sizeof(double) *
							   (map->ex - map->sx + 1 + 2 * map->hs) *
							   (map->ey - map->sy + 1 + 2 * map->hs));

	//if memory cannot be allocated

	if (myField == NULL)
	{
		printf("Error! Memory not allocated.\n");
		exit(-1);
	}
	return myField;
}

// Fill local map with values

void fillField(double *u, MAP *map)
{
	// each processor fill its part of the field
	for (int j = map->sy; j <= (map->ey); j++)
	{
		for (int i = map->sx; i <= (map->ex); i++)
		{
			// U(i, j) = (i + j / 10.0);
			U(i, j) = Mandelbrot(i, j, map->gex, map->gey);
		}
	}
}

double Mandelbrot(int i, int j, int gex, int gey) double x_max = 1;
double x_min = -2;
double y_max = 1.5;
double y_min = -1.5;

double dx = (x_max - x_min) / gex;
double dy = (y_max - y_min) / gey;

double counter = 0;
double iterations = 50;
double Re = 0, Im = 0;
double Re_aux, Im_aux;
double c_real = x_min + i * dx;
double c_imag = y_min + j * dy;
int diverge = 0;
double modulus, value;

while ((counter < iterations) && (!diverge))
{
	Re_aux = Re * Re - Im * Im + c_real;
	Im_aux = 2 * Re * Im + c_imag;

	modulus = Re_aux * Re_aux + Im_aux * Im_aux;

	if (modulus > 4)
	{
		diverge = 1;
	}

	Re = Re_aux;
	Im = Im_aux;

	counter++;
}

value = counter / iterations;

return value;
}

// Print local map values

void printField(double *u, MAP *map)
{
	for (int j = map->sy; j <= (map->ey); j++)
	{
		for (int i = map->sx; i <= (map->ex); i++)
		{
			printf("%lf ", U(i, j));
		}
		printf("\n");
	}
}

// Allocate memory for the global field

double *allocGlobalField(MAP *map)
{
	// allocs memory field
	// calculate the size of the field
	double *myField;
	// printf("debug: %d %d %d %d %d\n", map->sx, map->ex, map->sy, map->ey, map->hs);
	myField = (double *)malloc(sizeof(double) * (map->gex - map->gsx + 1) * (map->gey - map->gsy + 1));

	//if memory cannot be allocated
	if (myField == NULL)
	{
		printf("Error! Memory not allocated.\n");
		exit(-1);
	}
	return myField;
}

// Fill global field with each processor's values
// (each processor fill its part of the global field)

void fillGlobalField(double *u, double *glob, MAP *map)
{
	int r;				// for error checking
	int t1 = 0, t2 = 1; // tags, allow to identify messages
	MPI_Status st;		// check tag

	if (proc() == 0)
	{
		for (int j = map->sy; j <= (map->ey); j++)
		{
			for (int i = map->sx; i <= (map->ex); i++)
			{
				GLOB(i, j) = U(i, j);
				printf("%lf ", U(i, j));
			}
			printf("\n");
		}

		for (int i = 1; i <= (nproc() - 1); i++)
		{
			// Create a receive vector (sx, ex, sy, ey)
			int sr1[4];
			r = MPI_Recv(sr1, 4, MPI_INT, i, t2, MPI_COMM_WORLD, &st);
			checkr(r, "receive");
			// allocs memory of the receiving vector
			double *sr2;
			sr2 = (double *)malloc(sizeof(double) * (sr1[1] - sr1[0] + 1) * (sr1[3] - sr1[2] + 1));

			//if memory cannot be allocated
			if (sr2 == NULL)
			{
				printf("Error! Memory not allocated.\n");
				exit(-1);
			}

			r = MPI_Recv(sr2, (sr1[1] - sr1[0] + 1) * (sr1[3] - sr1[2] + 1), MPI_DOUBLE, i, t1, MPI_COMM_WORLD, &st);
			{
				GLOB(i2, j2) = SR(i2, j2);
				GLOB(i2, j2) = SR(i2, j2);
				// printf("%lf ",SR(i2,j2));
			}
			// printf("\n");
		}

		free(sr2);
	}
}

// proc != 0

else
{
	int vs[4];
	vs[0] = map->sx;
	vs[1] = map->ex;
	vs[2] = map->sy;
	vs[3] = map->ey;

	r = MPI_Send(u, (map->ex - map->sx + 1) * (map->ey - map->sy + 1), MPI_DOUBLE, 0 /*destination*/, t1, MPI_COMM_WORLD);
	checkr(r, "send1");
	r = MPI_Send(vs, 4, MPI_INT, 0 /*destination*/, t2, MPI_COMM_WORLD);
	checkr(r, "send2");
}
}

// Print global field values

void printGlobalField(double *glob, MAP *map)
{
	for (int j = map->gey; j >= map->gsy; j--)
	{
		for (int i = map->gsx; i <= map->gex; i++)
		{
			printf("%lf ", GLOB(i, j));
		}
		printf("\n");
	}
}
