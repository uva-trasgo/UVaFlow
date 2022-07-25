#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "rk4.h"
#include "interpolation.h"
/*
void runge_kutta_45 ( double *Pcoords, double t0, double tend, mesh_t *mesh, double *result, int nsteps, int it )
{ 
    int i;   
    double k1[mesh->nDim];
    double k2[mesh->nDim];
    double k3[mesh->nDim];
    double k4[mesh->nDim];
    double k5[mesh->nDim];
    double k6[mesh->nDim];
    double h = (tend - t0) / nsteps;
    double x0[mesh->nDim];
    double newP[mesh->nDim];

    x0[0] = Pcoords[0];
    x0[1] = Pcoords[1];
    if ( mesh->nDim == 3 ) x0[2] = Pcoords[2];

    for ( i = 0; i < nsteps; i++ )
    {
        // K1 
        linear_interpolation_approach2 ( t0, x0, mesh, k1, it );

        k1[0] = k1[0] * h;
        k1[1] = k1[1] * h;
        if ( mesh->nDim == 3 ) k1[2] = k1[2] * h;

        // K2 
        newP[0] = x0[0] + 1/4 * k1[0];
        newP[1] = x0[1] + 1/4 * k1[1];
        if ( mesh->nDim == 3 ) newP[2] = x0[2] + 1/4 * k1[2];

        linear_interpolation_approach2 ( t0 + 1/4 * h, newP, mesh, k2, it );
        k2[0] = k2[0] * h;
        k2[1] = k2[1] * h;
        if ( mesh->nDim == 3 ) k2[2] = k2[2] * h;

        // K3 
        newP[0] = x0[0] + 3/32 * k1[0] + 9/32 * k2[0];
        newP[1] = x0[1] + 3/32 * k1[1] + 9/32 * k2[1];
        if ( mesh->nDim == 3 ) newP[2] = x0[2] + 3/32 * k1[2] + 9/32 * k2[2];

        linear_interpolation_approach2 ( t0 + 3/8 * h, newP, mesh, k3, it );

        k3[0] = k3[0] * h;
        k3[1] = k3[1] * h;
        if ( mesh->nDim == 3 ) k3[2] = k3[2] * h;

        // K4 
        newP[0] = x0[0] + 1932/2197 * k1[0] - 7200/2197 * k2[0] + 7296/2197 * k3[0];
        newP[1] = x0[1] + 1932/2197 * k1[1] - 7200/2197 * k2[1] + 7296/2197 * k3[1];
        if ( mesh->nDim == 3 ) newP[2] = x0[2] + 1932/2197 * k1[2] - 7200/2197 * k2[2] + 7296/2197 * k3[2];

        linear_interpolation_approach2 ( t0 + 12/13 * h, newP, mesh, k4, it );

        k4[0] = k4[0] * h;
        k4[1] = k4[1] * h;
        if ( mesh->nDim == 3 ) k4[2] = k4[2] * h;

        // K5 
        newP[0] = x0[0] + 439/216 * k1[0] - 8 * k2[0] + 3680/513 * k3[0] - 845/4104 * k4[0];
        newP[1] = x0[1] + 439/216 * k1[1] - 8 * k2[1] + 3680/513 * k3[1] - 845/4104 * k4[1];
        if ( mesh->nDim == 3 ) newP[2] = x0[2] + 439/216 * k1[2] - 8 * k2[2] + 3680/513 * k3[2] - 845/4104 * k4[2];

        linear_interpolation_approach2 ( t0 + h, newP, mesh, k5, it );

        k5[0] = k5[0] * h;
        k5[1] = k5[1] * h;
        if ( mesh->nDim == 3 ) k5[2] = k5[2] * h;

        // K6 
        newP[0] = x0[0] - 8/27 * k1[0] + 2 * k2[0] - 3544/2565 * k3[0] + 1859/4104 * k4[0] - 11/40 * k5[0];
        newP[1] = x0[1] - 8/27 * k1[1] + 2 * k2[1] - 3544/2565 * k3[1] + 1859/4104 * k4[1] - 11/40 * k5[1];
        if ( mesh->nDim == 3 ) newP[2] = x0[2] - 8/27 * k1[2] + 2 * k2[2] - 3544/2565 * k3[2] + 1859/4104 * k4[2] - 11/40 * k5[2];

        linear_interpolation_approach2 ( t0 + 1/2 * h, newP, mesh, k6, it );

        k6[0] = k6[0] * h;
        k6[1] = k6[1] * h;
        if ( mesh->nDim == 3 ) k6[2] = k6[2] * h;
    }

    // Rk 4th order
    //result[0] = Pcoords[0] + 25/216 * k1[0] + 1408/2565 * k3[0] + 2197/4104 * k4[0] - 1/5 * k5[0];
    //result[1] = Pcoords[1] + 25/216 * k1[1] + 1408/2565 * k3[1] + 2197/4104 * k4[1] - 1/5 * k5[1];
    //if ( mesh->nDim == 3 ) result[2] = Pcoords[2] + 25/216 * k1[2] + 1408/2565 * k3[2] + 2197/4104 * k4[2] - 1/5 * k5[2];
    
    // Rk 5th order
    result[0] = Pcoords[0] + 16/135 * k1[0] + 6656/12825 * k3[0] + 28561/56430 * k4[0] - 9/50 * k5[0] + 2/55 * k6[0];
    result[1] = Pcoords[1] + 16/135 * k1[1] + 6656/12825 * k3[1] + 28561/56430 * k4[1] - 9/50 * k5[1] + 2/55 * k6[1];
    if ( mesh->nDim == 3 ) result[2] = Pcoords[2] + 16/135 * k1[2] + 6656/12825 * k3[2] + 28561/56430 * k4[2] - 9/50 * k5[2] + 2/55 * k6[2];
}
*/

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
