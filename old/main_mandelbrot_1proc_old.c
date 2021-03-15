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

// Macro to access field u and GLOB
#define U(x, y) *(u + (x - map->sx + map->hs) + (y - map->sy + map->hs) * (map->ex - map->sx + 1 + 2 * map->hs))
#define GLOB(x) *(glob + (x - gsx))

// This structure contains all the data to access a distributed 2d array
typedef struct Maps
{
    int sx, ex, sy, ey;
    int hs;
} MAP;

// Function prototypes
void checkr(int r, char *txt);                                                          // Check state
int proc();                                                                             // Rank of the actual processor
int nproc();                                                                            // Size of total number processors
void worksplit(int *mystart, int *myend, int proc, int nproc, int start, int end);      // Worksplit action
void createMap(int NPX, int NPY, int gsx, int gex, int gsy, int gey, int hs, MAP *map); // Each processor create its own local map
void printMap(MAP *map);                                                                // Print actual processor rank, mystart and myend
double *allocField(MAP *map);                                                           // Allocate memory for local map
double *allocGlobalField(int gsx, int gex, int gsy, int gey);                           // Allocate memory for the global field
void fillField(double *u, MAP *m);                                                      // Fill local map with values
void printField(double *u, MAP *map);                                                   // Print local map values
void fillGlobalField(double *u, MAP *map, double *glob, int gsx);                       // Fill global map with values
void printGlobalField(double *glob, int gsx, int gex, int gsy, int gey);                // Print global field values
double distance(double x, double y);

void compute(double x, double y, double c_real, double c_imag, double *ans_x, double *ans_y);
void mandelbrot(double *px, double *py, int iter, double c_real, double c_imag);

// Main function
int main(int argc, char **argv)
{

    // Declare variables
    int NPX = 10; // Number of processors in X axis
    int NPY = 1;  // Number of processors in Y axis
    int gsx = 1;  // Global X start index
    int gex = 50; // Global X end index
    int gsy = 1;  // Global Y start index
    int gey = 1;  // Global Y end index
    int hs = 0;   // Halo size

    // MPI variables
    int r; // Error checking

    // Data variables
    double *u, *glob;
    MAP map_;
    MAP *map = &map_;

    // Complex number
    double r_min, r_max;
    double y_min = -1, y_max = 1, y_div = 2000;
    double x_min = -2, x_max = 1, x_div = 3000;

    int iter = 50; // maximum iterations
    double px, py, c_real = x_min, c_imag = y_min;

    // Start MPI
    r = MPI_Init(&argc, &argv);
    checkr(r, "Initiate");

    // Create local map
    createMap(NPX, NPY,           // number of processors in each direction
              gsx, gex, gsy, gey, // GLOBAL limits
              hs,                 // Halo size
              map);

    // Allocate memory for local map
    // u = allocField(map); // Local vector (pointer)

    // Allocate memory for global map
    // glob = allocGlobalField(gsx, gex, gsy, gey); // Global vector (pointer)

    // // Fill map with values
    // fillField(u, map);

    // printf("Local map filled successfully! \n\n");

    double interval_x = (x_max - x_min) / x_div;
    double interval_y = (y_max - y_min) / y_div;

    for (int i = 1; i <= x_div; i++)
    {
        for (int j = 1; j <= y_div; j++)
        {
            // printf("c_real = %lf c_imag = %lf \n", c_real, c_imag);
            mandelbrot(&px, &py, iter, c_real, c_imag);

            if ((px != -2) && (py != -2))
            {
                //printf("px = %lf and py = % lf \n", px, py);
                printf("%lf %lf \n", px, py);
            }
            c_imag = c_imag + interval_y;
            // printf("c_imag = %lf \n", c_imag);
        }
        c_imag = y_min;
        c_real = c_real + interval_x;
        // printf("c_real = %lf \n", *p_c_real);
    }

    // printf("\nMandelbrot printed! \n\n");

    // Fill global field with each processor's values
    // fillGlobalField(u, map, glob, gsx);

    // printf("Global map filled successfully! \n\n");

    // Print global map

    // if (proc() == 0)
    // {
    // 	printGlobalField(glob, gsx, gex, gsy, gey);
    // 	printf("Global map printed successfully! \n\n");
    // }

    // Free allocated memory
    // free(u);
    // free(glob);

    // End MPI
    MPI_Finalize();

    // End main
    exit(0);
}

// EOF
/* -------------------------------------------------------------------------- */

// we use this function to check the return value of every MPI call

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

// each processor creates its map, filling in its values for sx,ex,sy,ey using worksplit

void createMap(int NPX, int NPY,                   // number of processors in each direction
               int gsx, int gex, int gsy, int gey, // GLOBAL limits
               int hs,                             // Halo size
               MAP *map)
{

    map->hs = hs;
    map->sy = gsy;
    map->ey = gey;
    worksplit(&(map->sx), &(map->ex), proc(), nproc(), gsx, gex);
}

// Displays the memory map on the screen

void printMap(MAP *map)
{
    printf("So, as I'm processor %d, I start with x%d and end with x%d\n", proc(), map->sx, map->ex);
}

// Memory allocation of each field (local limits)

double *allocField(MAP *map)
{

    double *aux;
    // printf("debug: %d %d %d %d %d\n", map->sx, map->ex, map->sy, map->ey, map->hs);
    aux = (double *)malloc(sizeof(double) * (map->ex - map->sx + 1 + 2 * map->hs) * (map->ey - map->sy + 1 + 2 * map->hs));

    //if memory cannot be allocated
    if (aux == NULL)
    {
        printf("Error! Memory not allocated.\n");
        exit(-1);
    }
    return aux;
}

// Allocate memory for the global field (global limits)

double *allocGlobalField(int gsx, int gex, int gsy, int gey)
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

// Fill global field with each processor's values

void fillGlobalField(double *u, MAP *map, double *glob, int gsx)
{
    // each processor fill its part of the global field
    int r;              // for error checking
    int t1 = 0, t2 = 1; // tag, allows to identify messages
    MPI_Status st;      // check tag

    if (proc() == 0)
    {
        printf("Global values\n");
        for (int j = map->sy; j <= (map->ey); j++)
        {
            for (int i = map->sx; i <= (map->ex); i++)
            {
                GLOB(i) = U(i, j);
                //printf("%lf ", GLOB(i, j, map->sx, map->ey, map->sy));
            }
            printf("\n");
        }

        for (int i = 1; i <= (nproc() - 1); i++)
        {
            // Create a receive vector (sx, ex, sy, ey)
            int vr1[4];
            r = MPI_Recv(vr1, 4, MPI_INT, i, t2, MPI_COMM_WORLD, &st);
            checkr(r, "receive");

            // allocs memory of the receiving vector
            double *vr2;
            vr2 = (double *)malloc(sizeof(double) * (vr1[1] - vr1[0] + 1) * (vr1[3] - vr1[2] + 1));

            //if memory cannot be allocated
            if (vr2 == NULL)
            {
                printf("Error! Memory not allocated.\n");
                exit(-1);
            }

            r = MPI_Recv(vr2, (vr1[1] - vr1[0] + 1) * (vr1[3] - vr1[2] + 1), MPI_DOUBLE, i, t1, MPI_COMM_WORLD, &st);
            checkr(r, "receive2");

            // printf("he rebut %lf %lf \n", *(vr2), *(vr2 + 1));

            for (int j2 = vr1[2]; j2 <= vr1[3]; j2++)
            {
                for (int i2 = vr1[0]; i2 <= vr1[1]; i2++)
                {
                    GLOB(i2) = *(vr2 + i2 - vr1[0]);
                    // printf("%lf ", GLOB(i2, j2, vr1[0], vr1[1], vr1[2]));
                }
                //printf("\n");
            }

            free(vr2);
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
        // printf("he enviat %lf %lf \n", U(3, 1), U(4, 1));
        // printf("he enviat %lf %lf \n", *(u), *(u + 1));

        r = MPI_Send(u, (map->ex - map->sx + 1) * (map->ey - map->sy + 1), MPI_DOUBLE, 0 /*destination*/, t1, MPI_COMM_WORLD);
        checkr(r, "send1");
        r = MPI_Send(vs, 4, MPI_INT, 0 /*destination*/, t2, MPI_COMM_WORLD);
        checkr(r, "send2");
    }
}

// Print global field values

void printGlobalField(double *glob, int gsx, int gex, int gsy, int gey)
{
    for (int j = gsy; j <= gey; j++)
    {
        for (int i = gsx; i <= gex; i++)
        {
            printf("%lf ", GLOB(i));
        }
        printf("\n");
    }
}

//

double distance(double x, double y)
{
    /* Computes the square of the distance to the origin.
		Receives: doubles x and y
		Sends: x^2+y^2
	*/
    return (x * x + y * y);
}

void compute(double x, double y, double c_real, double c_imag, double *ans_x, double *ans_y)
{
    // Calculates the n-iteration
    *ans_x = x * x - y * y + c_real;
    *ans_y = 2 * x * y + c_imag;
}

void mandelbrot(double *px, double *py, int iter, double c_real, double c_imag)
{
    int counter = 0;
    double x = 0, y = 0;
    int imlittle = 1;
    double ans_x;
    double ans_y;
    while ((counter < iter) && (imlittle == 1))
    {
        compute(x, y, c_real, c_imag, &ans_x, &ans_y);

        if (distance(x, y) > 4)
        {
            imlittle = 0;
            *px = -2;
            *py = -2;
        }

        x = ans_x;
        y = ans_y;

        counter++;
    }

    if (imlittle == 1)
    {
        *px = c_real;
        *py = c_imag;
    }
}