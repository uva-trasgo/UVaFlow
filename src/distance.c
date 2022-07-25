#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "distance.h"

double distance_PQ ( mesh_t *mesh, double *P, double *Q )
{
   double distance = (Q[0]-P[0]) * (Q[0]-P[0]) + (Q[1]-P[1]) * (Q[1]-P[1]);
   if (mesh->nDim == 3)
      distance += (Q[2]-P[2]) * (Q[2]-P[2]);
   return sqrt(distance);
}

