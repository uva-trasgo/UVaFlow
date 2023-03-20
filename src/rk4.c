#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "rk4.h"
#include "interpolation.h"
#include "kdtree.h"

void runge_kutta_4 ( double Pcoords_x, double Pcoords_y, double Pcoords_z, double t0, double tend, double *result, int nsteps, int nDim, int nPoints, int nTimes, double *times, int nVertsPerFace, int nFaces, int *faces, double *coords_x, double *coords_y, double *coords_z, double *velocities, struct kdtree *kd, int *nFacesPerPoint, int *facesPerPoint) //, int it )
{
   int i;
   struct timeval start;
   struct timeval end;
   double time;
   double h;

   double k1[nDim];
   double k2[nDim];
   double k3[nDim];
   double k4[nDim];
   double x0[nDim];
   double newP[nDim];


   h = (tend - t0) / nsteps;

   x0[0] = Pcoords_x; //[0];
   x0[1] = Pcoords_y; //[1];
   if ( nDim == 3 ) x0[2] = Pcoords_z; //[2]; 

   result[0] = Pcoords_x;//[0];
   result[1] = Pcoords_y;//[1];
   if ( nDim == 3 ) result[2] = Pcoords_z;//[2];

   for ( i = 0; i < nsteps; i++ )
   {
	/* K1 */
	if ( nDim == 2 ) 
	   	linear_interpolation_approach2_2D ( t0, x0, times, k1, nDim, nPoints, nTimes, nVertsPerFace, nFaces, faces, coords_x, coords_y, velocities, kd, nFacesPerPoint, facesPerPoint );
	else
		linear_interpolation_approach2_3D ( t0, x0, times, k1, nDim, nPoints, nTimes, nVertsPerFace, nFaces, faces, coords_x, coords_y, coords_z, velocities, kd, nFacesPerPoint, facesPerPoint );

   	k1[0] = k1[0] * h;
   	k1[1] = k1[1] * h;
	if ( nDim == 3 ) k1[2] = k1[2] * h;
 
	/* K2 */
   	newP[0] = x0[0] + 1/2 * k1[0];
   	newP[1] = x0[1] + 1/2 * k1[1];
        if ( nDim == 3 )
	{
		newP[2] = x0[2] + 1/2 * k1[2];
		linear_interpolation_approach2_3D ( t0+h/2.0, newP, times, k2, nDim, nPoints, nTimes, nVertsPerFace, nFaces, faces, coords_x, coords_y, coords_z, velocities, kd, nFacesPerPoint, facesPerPoint );
	}
	else
		linear_interpolation_approach2_2D ( t0+h/2.0, newP, times, k2, nDim, nPoints, nTimes, nVertsPerFace, nFaces, faces, coords_x, coords_y, velocities, kd, nFacesPerPoint, facesPerPoint );


   	k2[0] = k2[0] * h;
   	k2[1] = k2[1] * h;
	if ( nDim == 3 ) k2[2] = k2[2] * h;

	/* K3 */  
   	newP[0] = x0[0] + 1/2 * k2[0];
   	newP[1] = x0[1] + 1/2 * k2[1];
	if ( nDim == 3 )
	{
		newP[2] = x0[2] + 1/2 * k2[2];
		linear_interpolation_approach2_3D ( t0+h/2.0, newP, times, k3, nDim, nPoints, nTimes, nVertsPerFace, nFaces, faces, coords_x, coords_y, coords_z, velocities, kd, nFacesPerPoint, facesPerPoint );
	}
	else
		linear_interpolation_approach2_2D ( t0+h/2.0, newP, times, k3, nDim, nPoints, nTimes, nVertsPerFace, nFaces, faces, coords_x, coords_y, velocities, kd, nFacesPerPoint, facesPerPoint );


   	k3[0] = k3[0] * h;
  	k3[1] = k3[1] * h;
	if ( nDim == 3 ) k3[2] = k3[2] * h;

	/* K4 */
   	newP[0] = x0[0] + k3[0];
   	newP[1] = x0[1] + k3[1];
	if ( nDim == 3 )
	{
		newP[2] = x0[2] + k3[2];
		linear_interpolation_approach2_3D ( t0+h, newP, times, k4, nDim, nPoints, nTimes, nVertsPerFace, nFaces, faces, coords_x, coords_y, coords_z, velocities, kd, nFacesPerPoint, facesPerPoint );
	}
	else
		linear_interpolation_approach2_2D ( t0+h, newP, times, k4, nDim, nPoints, nTimes, nVertsPerFace, nFaces, faces, coords_x, coords_y, velocities, kd, nFacesPerPoint, facesPerPoint );


   	k4[0] = k4[0] * h;
        k4[1] = k4[1] * h;
	if ( nDim == 3 ) k4[2] = k4[2] * h;

	/* Result */
   	result[0] += 1/6.0 * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
   	result[1] += 1/6.0 * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
	if ( nDim == 3 ) result[2] += 1/6.0 * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);

	x0[0] = result[0];
	x0[1] = result[1];
	if ( nDim == 3 ) x0[2] = result[2];
	t0 += h;
   }
}
