#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "rk4.h"
#include "interpolation.h"

void runge_kutta_4 ( double *Pcoords, double t0, double tend, mesh_t *mesh, double *result, int nsteps) //, int it )
{
   int i;
   struct timeval start;
   struct timeval end;
   double time;
   double h;

   double k1[mesh->nDim];
   double k2[mesh->nDim];
   double k3[mesh->nDim];
   double k4[mesh->nDim];
   double x0[mesh->nDim];
   double newP[mesh->nDim];

   h = (tend - t0) / nsteps;

   x0[0] = Pcoords[0];
   x0[1] = Pcoords[1];
   if ( mesh->nDim == 3 ) x0[2] = Pcoords[2]; 

   result[0] = Pcoords[0];
   result[1] = Pcoords[1];
   if ( mesh->nDim == 3 ) result[2] = Pcoords[2];

   for ( i = 0; i < nsteps; i++ )
   {
	/* K1 */
   	linear_interpolation_approach2 ( t0, x0, mesh, k1 );
   	k1[0] = k1[0] * h;
   	k1[1] = k1[1] * h;
	if ( mesh->nDim == 3 ) k1[2] = k1[2] * h;
 
	/* K2 */
   	newP[0] = x0[0] + 1/2 * k1[0];
   	newP[1] = x0[1] + 1/2 * k1[1];
        if ( mesh->nDim == 3 ) newP[2] = x0[2] + 1/2 * k1[2];

   	linear_interpolation_approach2 ( t0+h/2.0, newP, mesh, k2 );

   	k2[0] = k2[0] * h;
   	k2[1] = k2[1] * h;
	if ( mesh->nDim == 3 ) k2[2] = k2[2] * h;

	/* K3 */  
   	newP[0] = x0[0] + 1/2 * k2[0];
   	newP[1] = x0[1] + 1/2 * k2[1];
	if ( mesh->nDim == 3 ) newP[2] = x0[2] + 1/2 * k2[2];

   	linear_interpolation_approach2 ( t0+h/2.0, newP, mesh, k3 );

   	k3[0] = k3[0] * h;
  	k3[1] = k3[1] * h;
	if ( mesh->nDim == 3 ) k3[2] = k3[2] * h;

	/* K4 */
   	newP[0] = x0[0] + k3[0];
   	newP[1] = x0[1] + k3[1];
	if ( mesh->nDim == 3 ) newP[2] = x0[2] + k3[2];

   	linear_interpolation_approach2 ( t0+h, newP, mesh, k4 );

   	k4[0] = k4[0] * h;
        k4[1] = k4[1] * h;
	if ( mesh->nDim == 3 ) k4[2] = k4[2] * h;

	/* Result */
   	result[0] += 1/6.0 * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
   	result[1] += 1/6.0 * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
	if ( mesh->nDim == 3 ) result[2] += 1/6.0 * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);

	x0[0] += h;
	x0[1] += h;
	if ( mesh->nDim == 3 ) x0[2] += h;
   }
	//printf("Point %f, %f, result %f, %f\n", Pcoords[0], Pcoords[1], result[0], result[1]);
/*
   free(k1);
   free(k2);
   free(k3);
   free(k4);
   free(newP);
   free(x0);
*/
}
