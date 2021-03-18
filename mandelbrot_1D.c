/* Mandelbrot set
Date: 18/03/2021
Author/s: Group 1
Subject: High Performance Computing for Aerospace Engineering
Professor: Manel Soria & Arnau Miro

// Problem statement
-------------------------------------------------------------------------
Plot Mandelbrot set using MPI using 1D array
-------------------------------------------------------------------------
*/

// Libraries
#include <stdio.h>	// Standard Input and Output Library
#include <stdlib.h> // Library for defining functions to perform general functions (malloc())
#include <math.h>	// Math library
#include "mpi.h"	// MPI library

// This structure contains all the data to access a distributed 2d array
typedef struct Maps
{
	// Start and end of X-axis and Y-axis
	int sx, ex, sy, ey;
	// Halo size
	int hs;
} MAP;

// Function prototypes
void checkr(int r, char *txt);																  // Check state
int proc();																					  // Rank of the actual processor
int nproc();																				  // Size of total number processors
void worksplit(int *mystart, int *myend, int proc, int nproc, int start, int end);			  // Worksplit action
void createMap(int gsx, int gex, int hs, MAP *map);											  // Each processor create its own local map
void printMap(MAP *map);																	  // Print actual processor rank, mystart and myend
double *allocField(MAP *map, int y_div);													  // Allocate memory for local map
int fillField(double *u, MAP *map, double x_div, double y_div, int proc);					  // Fill local map with values
double distance(double x, double y);														  // Squared distance between points (x,y) and origin (0,0)
void compute(double x, double y, double c_real, double c_imag, double *ans_x, double *ans_y); // Return the n-iteration of the Mandelbrot expression
int mandelbrot(int iter, double c_real, double c_imag);			  // Compute the set of numbers that are inside Mandelbrot's set
void printField(double *u, MAP *map, int counter);

// Main function
int main(int argc, char **argv)
{

	// Declare variables
	int gsx = 1;		 		 // Global X start index
	int gex = 10000;		 	 // Global X end index
	int hs = 0;			 		 // Halo size
	double x_div = gex - gsx + 1;
	double y_div = 10000; // Number of divisions in Y-axis


	// MPI variables
	int r; // Error checking

	// Data variables
	double *u; // Array of local map values
	MAP map_;		  // Map structure
	MAP *map = &map_; // Map pointer

	// Start MPI
	r = MPI_Init(&argc, &argv);
	checkr(r, "Initiate");

	// Create local map
	createMap(gsx, gex, hs, map);

	// printMap(map);

	// Allocate memory for local map
	u = allocField(map, y_div); // Local vector (pointer)

	// Count the amount of numbers inside Mandelbrot's set for each processor
	int counter;

	// Fill map with values and return the amount of numbers inside Mandelbrot's set for each processor
	counter = fillField(u, map, x_div, y_div, proc());

	// Print each processor data
	if (counter != 0) {printField(u, map, counter);}

	// printf("\nMandelbrot printed! \n\n");

	// Free allocated memory
	free(u);

	// End MPI
	MPI_Finalize();

	// End main
	exit(0);
}

// End of File
/* -------------------------------------------------------------------------- */

// Check the return value of every MPI call
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

// Each processor creates its map, filling in its values for sx,ex,sy,ey using worksplit
void createMap(int gsx, int gex, int hs, MAP *map)
{
	map->hs = hs;
	worksplit(&(map->sx), &(map->ex), proc(), nproc(), gsx, gex);
}

// Displays the memory map on the screen
void printMap(MAP *map)
{
	printf("So, as I'm processor %d, I start with x%d and end with x%d\n", proc(), map->sx, map->ex);
}

// Memory allocation of each field (local limits)
double *allocField(MAP *map, int y_div)
{
	double *aux;
	// printf("debug: %d %d %d %d %d\n", map->sx, map->ex, map->sy, map->ey, map->hs);
	aux = (double *)malloc(sizeof(double) * (2 * (map->ex - map->sx + 1) + 2 * map->hs) * (y_div + 2 * map->hs));

	//if memory cannot be allocated
	if (aux == NULL)
	{
		printf("Error! Memory not allocated.\n");
		exit(-1);
	}
	return aux;
}

// Fill local map with values
int fillField(double *u, MAP *map, double x_div, double y_div, int proc)
{
	double y_min = -1, y_max = 1;	// Limits of the area of study for the real part
	double x_min = -2, x_max = 0.5; // Limits of the area of study for the imaginary part
	int iter = 250;					// maximum iterations
	double interval_x = (x_max - x_min) / x_div;
	double interval_y = (y_max - y_min) / y_div;
	int convergence;
	double c_real = -2 + ((x_max - x_min)/x_div * (map->ex - map-> sx + 1)) * proc;
	double c_imag = y_min;
	int counter = 0;

	for (int i = map->sx; i <= map->ex; i++)
	{
		for (int j = 1; j < y_div; j++)
		{
			// printf("c_real = %lf c_imag = %lf \n", c_real, c_imag);
			convergence = mandelbrot(iter, c_real, c_imag);

			// If the point does not belong to the Mandelbrot set
			if (convergence == 1)
			{
				// printf("%lf %lf \n", px, py);
				*(u + counter) = c_real;
				*(u + counter + 1) = c_imag;
				counter = counter + 2;
			}
			c_imag = c_imag + interval_y;
		}
		c_imag = y_min;
		c_real = c_real + interval_x;
	}
	return counter;
}

// Squared distance between points (x,y) and origin (0,0)
double distance(double x, double y)
{
	/* Computes the square of the distance to the origin.
		Receives: doubles x and y
		Sends: x^2+y^2
	*/
	return (x * x + y * y);
}

// Return the n-iteration of the Mandelbrot
void compute(double x, double y, double c_real, double c_imag, double *ans_x, double *ans_y)
{
	// Calculates the n-iteration
	*ans_x = x * x - y * y + c_real; // Real part from the Mandelbrot expression
	*ans_y = 2 * x * y + c_imag;	 // Imaginary part from the Mandelbrot expression
}

// Compute the set of numbers that are inside Mandelbrot's set
int mandelbrot(int iter, double c_real, double c_imag)
{
	int counter = 0;
	double x = 0, y = 0;
	double ans_x;
	double ans_y;
	int loc_convergence = 1;

	while ((counter < iter) && (loc_convergence != 0))
	{
		compute(x, y, c_real, c_imag, &ans_x, &ans_y);

		if (distance(x, y) > 4) // The complex number is out of the Mandelbrot set
		{
			loc_convergence = 0;
		}

		x = ans_x; // Updates the real part for the next iteration
		y = ans_y; // Updates the imaginary part for the next iteration

		counter++;
	}

	return loc_convergence;

}

void printField(double *u, MAP *map, int counter)
{
		// MPI variables
		int r; // Error checking
		int p;
		for (p = 0; p < nproc(); p++)
		{
				if (proc() == p)
				{
						for (int i = 0; i < counter; i=i+2)
						{
								printf("%lf %lf \n",*(u+i),*(u+i+1));
						}
				}
		// Waits until all the processors arrive to the barrier
		r = MPI_Barrier(MPI_COMM_WORLD); checkr(r,"Barrier");
		}
}
