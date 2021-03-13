/* Mandelbrot set

<<<<<<< HEAD:main.c
Date: 15/03/2021
=======
Date: 27/02/2021
>>>>>>> 96cf66f87d7383e3ac9e75d1341c2718f1377163:main_add.c
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
#include <math.h>	  // Math library
#include "mpi.h"	  // MPI library

<<<<<<< HEAD:main.c
// Macro to access field u and GLOB
=======
// Macro to access field u, v and w and GLOB
>>>>>>> 96cf66f87d7383e3ac9e75d1341c2718f1377163:main_add.c
#define U(x, y) *(u + (x - map->sx + map->hs) + (y - map->sy + map->hs) * (map->ex - map->sx + 1 + 2 * map->hs))
#define GLOB(x, y) *(glob + (x - map->sx) + (y - map->sy) * (map->ex - map->sx + 1))

// This structure contains all the data to access a distributed 2d array
typedef struct Maps {
    int sx, ex, sy, ey;
    int hs;
} MAP;

<<<<<<< HEAD:main.c
// Function prototypes
void checkr(int r, char *txt);                                                            // Check state
int proc();                                                                               // Rank of the actual processor
int nproc();                                                                              // Size of total number processors
void worksplit(int *mystart, int *myend, int proc, int nproc, int start, int end);        // Worksplit action
void createMap(int NPX, int NPY, int gsx, int gex, int gsy, int gey, int hs, MAP *map);   // Each processor create its own local map
void printMap(MAP *map);                                                                  // Print actual processor rank, mystart and myend
double *allocField(MAP *map);                                                             // Allocate memory for local map
double *allocGlobalField(int gsx, int gex, int gsy, int gey);							                // Allocate memory for the global field                                                         // Allocate memory for local map
void fillField(double *u, MAP *m);														                            // Fill local map with values
void printField(double *u, MAP *map);													                            // Print local map values

// Main function
int main(int argc, char **argv) {

    // Declare variables
    int NPX = 2; // Number of processors in X axis
  	int NPY = 1; // Number of processors in Y axis
  	int gsx = 1; // Global X start index
  	int gex = 4; // Global X end index
  	int gsy = 1; // Global Y start index
  	int gey = 1; // Global Y end index
  	int hs = 2;	 // Halo size

    // MPI variables
    int r; // Error checking

    // Data variables
    double *u, *glob;
    MAP map_;
    MAP *map = &map_;

    // Start MPI
    r = MPI_Init(&argc, &argv);
      checkr(r,"Initiate");

    // Create local map
    createMap(NPX, NPY,			  // number of processors in each direction
			  gsx, gex, gsy, gey,   // GLOBAL limits
			  hs,				            // Halo size
			  map);

    // Print local processor rank, mystart and myend
    printMap(map);

    // Allocate memory for local map
    u = allocField(map); // Local vector (pointer)
=======
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
>>>>>>> 96cf66f87d7383e3ac9e75d1341c2718f1377163:main_add.c

    // Allocate memory for global map
	  glob = allocGlobalField(gsx, gex, gsy, gey); // Global vector (pointer)

<<<<<<< HEAD:main.c
    printf("\nMemory allocated successfully! \n\n");

    // Fill map with values
    fillField(u, map);

    printf("Local map filled successfully! \n\n");
=======
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
>>>>>>> 96cf66f87d7383e3ac9e75d1341c2718f1377163:main_add.c

    // Print local map
    printField(u, map);

<<<<<<< HEAD:main.c
    printf("\nLocal map printed successfully! \n\n");

    // Fill global field with each processor's values
	  //fillGlobalField(u, map, glob);

    printf("GLobal map filled successfully! \n\n");
=======
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
>>>>>>> 96cf66f87d7383e3ac9e75d1341c2718f1377163:main_add.c

    // Free allocated memory
    free(u);
    free(glob);

    // End MPI
    MPI_Finalize();

    // End main
    exit(0);
}

<<<<<<< HEAD:main.c
// EOF
/* -------------------------------------------------------------------------- */

// we use this function to check the return value of every MPI call

void checkr(int r, char *txt) {
	if (r != MPI_SUCCESS) {
=======
// Action to check the return value of every MPI call
void checkr(int r, char *txt)
{
	if (r != MPI_SUCCESS)
	{
>>>>>>> 96cf66f87d7383e3ac9e75d1341c2718f1377163:main_add.c
		fprintf(stderr, "Error: %s\n", txt);
		exit(-1);
	}
}

<<<<<<< HEAD:main.c
// Indicates the rank of the processor

int proc(){
	int r, rank;
	r = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	checkr(r, "proc");
	return (rank);
}

// Indicates the number of processors being used

int nproc() {
	int r, size;
	r = MPI_Comm_size(MPI_COMM_WORLD, &size);
	checkr(r, "nproc");
	return (size);
}

// Divide tasks between processors

void worksplit(int *mystart, int *myend, int proc, int nproc, int start, int end) {
    // Number of tasks
    int ntask = end-start+1;
    // Number of tasks per processor
    int interval = ntask/nproc;
    // Tasks left
    int remainder = ntask%nproc;

    if (ntask < nproc) {
        printf("Less tasks than processors\n");
        exit(-1);
    }

    if (remainder != 0) {
        if(proc<remainder) {
          *mystart = start + proc * (interval +1);
          *myend = *mystart + interval;
        }
    else {
        *mystart = start + remainder*(interval+1)+(proc-remainder)*interval;
        *myend = *mystart + interval -1;
    }
  }
  else {
      *mystart = start + proc * interval;
      *myend = *mystart + (interval -1);
  }
=======
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
>>>>>>> 96cf66f87d7383e3ac9e75d1341c2718f1377163:main_add.c
}

// each processor creates its map, filling in its values for sx,ex,sy,ey using worksplit

void createMap(int NPX, int NPY,				     // number of processors in each direction
			   int gsx, int gex, int gsy, int gey, // GLOBAL limits
			   int hs,							               // Halo size
			   MAP *map) {

     map->hs = hs;
	   map->sy = gsy;
	   map->ey = gey;
	   worksplit(&(map->sx), &(map->ex), proc(), nproc(), gsx, gex);
}

<<<<<<< HEAD:main.c
// Displays the memory map on the screen

void printMap(MAP *map) {
	printf("So, as I'm processor %d, I start with x%d and end with x%d\n", proc(), map->sx, map->ex);
}

// Memory allocation of each field (local limits)

double *allocField(MAP *map) {

=======
// Allocate memory for the global field
double *allocGlobalMap(int gsx, int gex, int gsy, int gey) // GLOBAL limits
{
	// allocs memory field
	// calculate the size of the field
>>>>>>> 96cf66f87d7383e3ac9e75d1341c2718f1377163:main_add.c
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

<<<<<<< HEAD:main.c
// Allocate memory for the global field (global limits)

double *allocGlobalField(int gsx, int gex, int gsy, int gey) {
=======
// Print actual processor rank, mystart and myend
void printMap(MAP *map)
{
	printf("So, as I'm processor %d, I start with x%d and end with x%d\n", proc(), map->sx, map->ex);
}

// Allocate memory for local map
double *allocField(MAP *map)
{
>>>>>>> 96cf66f87d7383e3ac9e75d1341c2718f1377163:main_add.c
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

// Fill local map with values
<<<<<<< HEAD:main.c

void fillField(double *u, MAP *map) {
=======
void fillField(double *u, MAP *map)
{
>>>>>>> 96cf66f87d7383e3ac9e75d1341c2718f1377163:main_add.c
	// each processor fill its part of the field
	for (int j = map->sy; j <= (map->ey); j++)
	{
		for (int i = map->sx; i <= (map->ex); i++)
		{
			U(i, j) = (i + j / 10.0);
		}
	}
}

<<<<<<< HEAD:main.c
// Print local map values

void printField(double *u, MAP *map) {
=======
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
			r = MPI_Recv((u + map->hs), (map->ex - map->sx + 1) * (map->ey - map->sy + 1), MPI_DOUBLE, i, t1, MPI_COMM_WORLD, &st);
			checkr(r, "receive1");
			r = MPI_Recv(&vr2[0], 4, MPI_INT, i, t2, MPI_COMM_WORLD, &st);
			checkr(r, "receive2");

			printf("he rebut %lf %lf \n", *(u + map->hs), *(u + map->hs + 1));

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
		printf("no soc 0\n");
		int vs[4];
		vs[0] = map->sx;
		vs[1] = map->ex;
		vs[2] = map->sy;
		vs[3] = map->ey;
		printf("he enviat %lf %lf \n", U(3, 1), U(4, 1));
		r = MPI_Send((u + map->hs), (map->ex - map->sx + 1) * (map->ey - map->sy + 1), MPI_DOUBLE, 0 /*destination*/, t1, MPI_COMM_WORLD);
		checkr(r, "send1");
		r = MPI_Send(&vs[0], 4, MPI_INT, 0 /*destination*/, t2, MPI_COMM_WORLD);
		checkr(r, "send2");
		printf("t'has equivocat noi. 0 kelvin. \n");
	}
}

// Print local map values
void printField(double *u, MAP *map)
{
	// each processor prints its part of the field
>>>>>>> 96cf66f87d7383e3ac9e75d1341c2718f1377163:main_add.c
	for (int j = map->sy; j <= (map->ey); j++)
	{
		for (int i = map->sx; i <= (map->ex); i++)
		{
			printf("%lf ", U(i, j));
		}
		printf("\n");
	}
}
<<<<<<< HEAD:main.c
=======

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
>>>>>>> 96cf66f87d7383e3ac9e75d1341c2718f1377163:main_add.c
