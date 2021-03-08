#include <stdio.h>
#include <stdlib.h>
//#include "mpi.h"
// This structure contains all the data to access a distributed 2d array
typedef struct Maps {
int sx,ex,sy,ey;
int hs;
//.... maybe we will need more data
} MAP;
// Macro to access field u
// define U(i,j) *( u + ... (map->sx)...)
// define V(i,j)
void createMap(int NPX, int NPY, // number of processors in each direction
int gsx,int gex,int gsy,int gey, // GLOBAL limits
int hs, // Halo size
MAP *map) {
// each processor creates its map, filling in its values for sx,ex,sy,ey using worksplit
// a worksplit for each direction has to be done
// map->sx = ... ;
// map->sy = ... ;
}

void printMap(MAP *map) { // prints my memory map
}
double *allocField(MAP* map) {
// allocs memory field
// calculate the size of the field
}
void printField(double *u,MAP* map) {
// each processor prints its part of the field
}
// u=v+w
void addField(double *u, double *v, double *w,MAP* map) {
// each processor sweeps its part of the field
// int i,j;
// for (j=map->sy; j<=map->ey; j++)
// for (i=map->sx; i<=map->ex; i++)
// U(i,j)=V(i,j)+W(i,j);
}
int main() {
double *u,*v,*w;
MAP map_;
MAP* map=&map_;
createMap(1,12,1,12,2,map);
u=allocField(map);
v=allocField(map);
w=allocField(map);
printField(u,map);
free(u);
}
