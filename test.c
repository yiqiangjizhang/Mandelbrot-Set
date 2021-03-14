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
#include <stdio.h>  // Standard Input and Output Library
#include <stdlib.h> // Library for defining functions to perform general functions (malloc())
#include <math.h>   // Math library
#include "mpi.h"    // MPI library

// Prototypes
void checkr(int r, char *txt);
int quants();
int whoami();
double *linspace(double start, double end, int n);

// Main
int main(int argc, char **argv)
{
    // MPI variables
    int r; // Error checking
    double *u;
    u = linspace(1, 10, 10);
    for (int i = 0; i < 10; i++)
    {
        printf("%f\n", *(u + i));
    }
    return 0;
}

// Check
void checkr(int r, char *txt)
{
    if (r != MPI_SUCCESS)
    {
        fprintf(stderr, "Error: %s\n", txt);
        exit(-1);
    }
}

// Size of the processors
int quants()
{
    int a, b;
    a = MPI_Comm_size(MPI_COMM_WORLD, &b);
    checkr(a, "quants");
    return (b);
}

// Rank of the processors
int whoami()
{
    int a, b;
    a = MPI_Comm_rank(MPI_COMM_WORLD, &b);
    checkr(a, "whoami");
    return (b);
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

// Create vector
double *linspace(double start, double end, int n)
{
    // Step size
    double stepsize;
    double *u;

    u = (double *)malloc(sizeof(double) * n);

    // Step size
    stepsize = (end - start) / (n - 1);

    for (int i = 0; i <= n; i++)
    {
        u[i] = start + stepsize * i;
    }
    return u;
}
