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
void worksplit(int *mystart, int *myend, int proc, int nproc, int start, int end);

// Main
int main(int argc, char **argv)
{
    // MPI variables
    int r; // Error checking

    // X axis
    double *u;
    // Y axis
    double *v;

    // Variables
    int start = 0;
    int end = 10;
    int mystart;
    int myend;
    int nproc = 2;

    // Number of points
    int n = 10;
    double sx = -2;
    double ex = 2;
    double sy = -1.5;
    double ey = 1.5;

    // Initate MPI
    r = MPI_Init(&argc, &argv);
    checkr(r, "init");

    u = (double *)malloc(sizeof(double) * n);
    v = (double *)malloc(sizeof(double) * n);

    u = linspace(sx, ex, n);
    v = linspace(sy, ey, n);

    for (int i = 0; i < 10; i++)
    {
        printf("%f\n", *(u + i));
    }

    for (int i = 0; i < 10; i++)
    {
        printf("%f\n", *(v + i));
    }

    // if (whoami() == 0)
    // {
    // worksplit(&mystart, &myend, whoami(), quants(), start, end);

    // // for (int i = 0; i < n; i++)
    // // {
    // // printf("%f\n", *(u + i));
    // // }
    // }
    // else
    // {
    // // worksplit(&mystart, &myend, whoami(), quants(), start - 1, end);
    // // for (int i = 0; i < n; i++)
    // // {
    // // printf("%f\n", *(u + i));
    // // }
    // }

    // // Loop through all the processors
    // for (int proc = 0; (proc <= nproc - 1) && (proc < end - start + 1); proc++)
    // {
    // worksplit(&mystart, &myend, proc, nproc, start, end);
    // }
    // for (int proc = end - start + 1; proc <= nproc - 1; proc++)
    // {
    // printf("So, as I'm processor %d, all tasks have been asigned and I can idle \n", proc);
    // }

    free(u);
    free(v);
    MPI_Finalize();
    exit(0);
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
    // Print off a hello world message
    printf("So, as I'm processor %d, I start with %d and end with %d\n", proc, *mystart, *myend);
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
