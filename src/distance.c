#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "distance.h"

double distance_PQ_2D ( int nDim, double Px, double Py, double *Q )
{
   return sqrt( (Q[0]-Px) * (Q[0]-Px) + (Q[1]-Py) * (Q[1]-Py) );
}

double distance_PQ_3D ( int nDim, double Px, double Py, double Pz, double *Q )
{
   return sqrt( (Q[0]-Px) * (Q[0]-Px) + (Q[1]-Py) * (Q[1]-Py) + (Q[2]-Pz) * (Q[2]-Pz) );
}

