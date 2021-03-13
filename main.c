/* Mandelbrot set

Date: 15/03/2021
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

// Macro to access field u and GLOB
#define U(x, y) *(u + (x - map->sx + map->hs) + (y - map->sy + map->hs) * (map->ex - map->sx + 1 + 2 * map->hs))
#define GLOB(x, y) *(glob + (x - map->sx) + (y - map->sy) * (map->ex - map->sx + 1))

// This structure contains all the data to access a distributed 2d array
typedef struct Maps {
    int sx, ex, sy, ey;
    int hs;
} MAP;

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

    // Allocate memory for global map
	  glob = allocGlobalField(gsx, gex, gsy, gey); // Global vector (pointer)

    printf("\nMemory allocated successfully! \n\n");

    // Fill map with values
    fillField(u, map);

    printf("Local map filled successfully! \n\n");

    // Print local map
    printField(u, map);

    printf("\nLocal map printed successfully! \n\n");

    // Fill global field with each processor's values
	  //fillGlobalField(u, map, glob);

    printf("GLobal map filled successfully! \n\n");

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

// we use this function to check the return value of every MPI call

void checkr(int r, char *txt) {
	if (r != MPI_SUCCESS) {
		fprintf(stderr, "Error: %s\n", txt);
		exit(-1);
	}
}

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

// Displays the memory map on the screen

void printMap(MAP *map) {
	printf("So, as I'm processor %d, I start with x%d and end with x%d\n", proc(), map->sx, map->ex);
}

// Memory allocation of each field (local limits)

double *allocField(MAP *map) {

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

// Allocate memory for the global field (global limits)

double *allocGlobalField(int gsx, int gex, int gsy, int gey) {
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

void fillField(double *u, MAP *map) {
	// each processor fill its part of the field
	for (int j = map->sy; j <= (map->ey); j++)
	{
		for (int i = map->sx; i <= (map->ex); i++)
		{
			U(i, j) = (i + j / 10.0);
		}
	}
}

// Print local map values

void printField(double *u, MAP *map) {
	for (int j = map->sy; j <= (map->ey); j++)
	{
		for (int i = map->sx; i <= (map->ex); i++)
		{
			printf("%lf ", U(i, j));
		}
		printf("\n");
	}
}
