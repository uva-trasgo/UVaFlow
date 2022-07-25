#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "interpolation.h"
#include "location.h"
#include "distance.h"

int reset_coordinates ( mesh_t *mesh, double *Pcoords )
{
   double distance, min_distance;
   int ip, ipmin;
   min_distance = distance_PQ ( mesh, mesh->points[0].coordinates, Pcoords );
   ipmin = 0;
   for ( ip = 1; ip < mesh->nPoints; ip++ )
   {
      distance = distance_PQ ( mesh, mesh->points[ip].coordinates, Pcoords );
      if ( distance < min_distance )
      {
         ipmin = ip;
         min_distance = distance;
      }
   }
   return ipmin;
}

void interpolate_triangle (   double *V1, double *velV1, 
                              double *V2, double *velV2, 
                              double *V3, double *velV3,
                              double *P, mesh_t *mesh,
                              double *interpolated_vel )
{
   double denom, wv1, wv2, wv3;

   if (mesh->nDim == 3)
   {
      fprintf(stderr, "nDim == 3 not supported in interpolate yet\n");
      exit(-1);
   }

   denom = (V2[1]-V3[1]) * (V1[0]-V3[0]) + (V3[0]-V2[0]) * (V1[1]-V3[1]);
   wv1   = (V2[1]-V3[1]) * (P[0]-V3[0])  + (V3[0]-V2[0]) * (P[1]-V3[1]);
   wv1   = wv1 / denom;
   wv2   = (V3[1]-V1[1]) * (P[0]-V3[0])  + (V1[0]-V3[0]) * (P[1]-V3[1]);
   wv2   = wv2 / denom;
   wv3   = 1 - wv1 - wv2;

   interpolated_vel[0] = wv1 * velV1[0] + wv2 * velV2[0] + wv3 * velV3[0];
   interpolated_vel[1] = wv1 * velV1[1] + wv2 * velV2[1] + wv3 * velV3[1];
}

void interpolate_3D_tetrahedral ( double *V1, double *velV1,
                      double *V2, double *velV2,
                      double *V3, double *velV3,
                      double *V4, double *velV4,
                      double *P, mesh_t *mesh,
                      double *interpolated_vel )
{
	int d;
	double num, denom;
	double lambda[mesh->nDim + 1];

	/* 1. Compute lambda values (barycentric coordinates) */
	
	/* 1.1 Compute lambda3 (located at lambda[2]) */
	num = 	(V4[0]-P[0]) * ( V1[1]*(V4[2]-V2[2]) + V2[1]*(V1[2]-V4[2]) + V4[1]*(V2[2]-V1[2]) ) +
			(V4[1]-P[1]) * ( V1[0]*(V2[2]-V4[2]) + V2[0]*(V4[2]-V1[2]) + V4[0]*(V1[2]-V2[2]) ) +
			(V4[2]-P[2]) * ( V1[0]*(V4[1]-V2[1]) + V2[0]*(V1[1]-V4[1]) + V4[0]*(V4[1]-V1[1]) );
	denom = (V3[0]-V4[0]) * ( V1[1]*(V2[2]-V4[2]) + V2[1]*(V4[2]-V1[2]) + V4[1]*(V1[2]-V2[2]) ) +
			(V2[0]-V4[0]) * ( V1[1]*(V4[2]-V3[2]) + V3[1]*(V1[2]-V4[2]) + V4[1]*(V3[2]-V1[2]) ) +
			(V1[0]-V4[0]) * ( V2[1]*(V3[2]-V4[2]) + V3[1]*(V4[2]-V2[2]) + V4[1]*(V2[2]-V3[2]) );
	lambda[2] = num/denom;

	/* 1.2 Compute lambda2 (located at lambda[1]) */
	num =	-( V1[1] - V4[1] ) * ( P[0] - lambda[2] * V3[0] ) +
			V1[0] * ( -lambda[2] * V3[1] + lambda[2] * V4[1] + P[1] + V4[1] ) + 
			V4[0] * ( -lambda[2] * V1[1] + lambda[2] * V3[1] - P[1] - V1[1] );
	denom = V1[0]*(V2[1]-V4[1]) + V2[0]*(V4[1]-V1[1]) + V4[0]*(V1[1]-V2[1]);
	lambda[1] = num / denom;

	/* 1.3 Compute lambda1 (located at lambda[0]) */
	lambda[0] = ( P[0] - (V2[0]-V4[0]) * lambda[1] - (V3[0]-V4[0]) * lambda[2] - V4[0] ) / (V1[0]-V4[0]);

	/* 1.4 Compute lambda4 (located at lambda[3]) */
	lambda[3] = 1 - lambda[0] - lambda[1] - lambda[2];

	/* 2. Interpolate veocity */
	for ( d = 0; d < mesh->nDim; d++ )
	{
		interpolated_vel[d] += 	lambda[0] * velV1[d] +  
								lambda[1] * velV2[d] +
								lambda[2] * velV3[d] +
								lambda[3] * velV4[d];
	}
}

void linear_interpolation_approach2 ( double t, double *Pcoords, mesh_t *mesh, double *interpolation)
{
   int ip, iface, itime, itprev, itpost, tsearch = 1;
   double t0, t1;
   double v0[mesh->nDim];
   double v1[mesh->nDim];
   double interp1[mesh->nDim];
   double interp2[mesh->nDim];
   double interp3[mesh->nDim];

   ip = find_point_in_mesh ( Pcoords, mesh );
   if ( ip == -1 ) // P is not a mesh point
   {
      iface = find_point_simplex_in_mesh ( Pcoords, mesh );
      if (iface > -1) // There is a mesh face or volume containing given point
      {
         // Find closest time values in mesh->times array and interpolate
	 for ( itime = 0; itime < mesh->nTimes && tsearch; itime++ )
	 {
	    if ( mesh->times[itime] == t ) // t exists in the mesh times elements
            {
		if ( mesh->nDim == 2 )
		{
                	interpolate_triangle (  mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].coordinates,  &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].velocity[mesh->nDim*itime]),
                                        mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].velocity[mesh->nDim*itime]),
                                        mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].velocity[mesh->nDim*itime]),
                                        Pcoords, mesh, interpolation );
		} 
		else // mesh->nDim == 3
		{
			interpolate_3D_tetrahedral ( mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].coordinates,  &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].velocity[mesh->nDim*itime]),
                                        mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].velocity[mesh->nDim*itime]),
                                        mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].velocity[mesh->nDim*itime]),
					mesh->points[mesh->faces[iface*mesh->nVertsPerFace+3]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+3]].velocity[mesh->nDim*itime]),
                                        Pcoords, mesh, interpolation );
		}
		tsearch = 0;
            }
            else
            {
               if ( mesh->times[itime] > t ) // There exist 2 times in provided mesh data surrounding given t
               {
                 itprev = itime-1;
                 itpost = itime;
                 t0  = mesh->times[itprev];
                 t1  = mesh->times[itpost];

		 if ( mesh->nDim == 2 )
                 {
			 interpolate_triangle (  mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].coordinates,   &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].velocity[mesh->nDim*itprev]),
                                 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].velocity[mesh->nDim*itprev]),
                                 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].velocity[mesh->nDim*itprev]),
                                 Pcoords, mesh, interp1 );

		         interpolate_triangle (  mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].coordinates,   &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].velocity[mesh->nDim*itpost]),
                                 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].velocity[mesh->nDim*itpost]),
                                 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].velocity[mesh->nDim*itpost]),
                                 Pcoords, mesh, interp2 );
		 }
		 else // mesh->nDim == 3
		 {
			interpolate_3D_tetrahedral (  mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].coordinates,   &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].velocity[mesh->nDim*itprev]),
                                 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].velocity[mesh->nDim*itprev]),
                                 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].velocity[mesh->nDim*itprev]),
				 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+3]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+3]].velocity[mesh->nDim*itprev]),
                                 Pcoords, mesh, interp1 );

                         interpolate_3D_tetrahedral (  mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].coordinates,   &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].velocity[mesh->nDim*itpost]),
                                 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].velocity[mesh->nDim*itpost]),
                                 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].velocity[mesh->nDim*itpost]),
				 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+3]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+3]].velocity[mesh->nDim*itpost]),
                                 Pcoords, mesh, interp2 );

		 }

	         interpolation[0] = (interp1[0] * (t1-t) + interp2[0] * (t-t0))/(t1-t0);
        	 interpolation[1] = (interp1[1] * (t1-t) + interp2[1] * (t-t0))/(t1-t0);
		 if ( mesh->nDim == 3 ) interpolation[2] = (interp1[2] * (t1-t) + interp2[2] * (t-t0))/(t1-t0);

                 tsearch = 0;
               }
            }
         }
         if (tsearch) // There is not enough data to interpolate
         {
            fprintf( stderr, "Error: There is not enough information to perform linear_interpolation for t=%f\n", t );
            exit(-1);
         }
      }
      else // The point coords are outside the given mesh data
      {
         // Pcoords will be reset to the closest mesh point coords
         ip = reset_coordinates ( mesh, Pcoords );
      }
   }
   if ( ip > -1 ) // Either the point is a mesh point or has been reseted to one of the mesh points
   {
      // Find closest time values in mesh->times array and interpolate
      for ( itime = 0; itime < mesh->nTimes && tsearch; itime++ )
      {
         if ( mesh->times[itime] == t )
         {
            interpolation[0] = mesh->points[ip].velocity[itime*mesh->nDim];
            interpolation[1] = mesh->points[ip].velocity[itime*mesh->nDim+1];
	    if ( mesh->nDim == 3 ) interpolation[2] = mesh->points[ip].velocity[itime*mesh->nDim+2];
            tsearch = 0;
         }
         else
         {
            if ( mesh->times[itime] > t ) // There exist 2 times in provided mesh data surrounding given t
            {
               itprev = itime-1;
               itpost = itime;
               t0  = mesh->times[itprev];
               t1  = mesh->times[itpost];

               v0[0] = mesh->points[ip].velocity[itprev*mesh->nDim];
               v0[1] = mesh->points[ip].velocity[itprev*mesh->nDim+1];
	       if ( mesh->nDim == 3 ) v0[2] = mesh->points[ip].velocity[itprev*mesh->nDim+2];

               v1[0] = mesh->points[ip].velocity[itpost*mesh->nDim];
               v1[1] = mesh->points[ip].velocity[itpost*mesh->nDim+1];
	       if ( mesh->nDim == 3 ) v1[2] = mesh->points[ip].velocity[itpost*mesh->nDim+2];

               interpolation[0] = (v0[0] * (t1-t) + v1[0] * (t-t0))/(t1-t0);
               interpolation[1] = (v0[1] * (t1-t) + v1[1] * (t-t0))/(t1-t0);
               if ( mesh->nDim == 3 ) interpolation[2] = (v0[2] * (t1-t) + v1[2] * (t-t0))/(t1-t0);

               tsearch = 0;
            }
         }
      }
      if (tsearch) // There is not enough data to interpolate
      {
         fprintf( stderr, "Error: There is not enough information to perform linear_interpolation for t=%f\n", t );
         exit(-1);
      }
   }
}

void linear_interpolation_approach1 ( double t, double *Pcoords, mesh_t *mesh, double *interpolation)
{
   int ip, iface, itime, itprev, itpost, tsearch = 1;
   double t0, t1;
   double v0[mesh->nDim];
   double v1[mesh->nDim];
   double interp1[mesh->nDim];
   double interp2[mesh->nDim];
   double interp3[mesh->nDim];
   double interp4[mesh->nDim];

   ip = find_point_in_mesh ( Pcoords, mesh );
   if ( ip == -1 ) // P is not a mesh point
   {
      iface = find_point_simplex_in_mesh ( Pcoords, mesh );
      if (iface > -1) // There is a mesh face containing given point
      {
         linear_interpolation_approach1 ( t, mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].coordinates,   mesh, interp1 );
         linear_interpolation_approach1 ( t, mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].coordinates, mesh, interp2 );
         linear_interpolation_approach1 ( t, mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].coordinates, mesh, interp3 );
	 if ( mesh->nDim == 3 ) linear_interpolation_approach1 ( t, mesh->points[mesh->faces[iface*mesh->nVertsPerFace+3]].coordinates, mesh, interp4 );

	 if ( mesh->nDim == 2)
 	 {
         	interpolate_triangle (  mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].coordinates,   interp1,
                                 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].coordinates, interp2, 
                                 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].coordinates, interp3, 
                                 Pcoords, mesh, interpolation );
	 }
	 else // mesh->nDim == 3
 	 {
		interpolate_3D_tetrahedral (  mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].coordinates,   interp1, 
                                 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].coordinates, interp2, 
                                 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].coordinates, interp3, 
				 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+3]].coordinates, interp4,
                                 Pcoords, mesh, interpolation );
	 }
      }
      else // The point coords are outside the given mesh data
      {
         // Pcoords will be reset to the closest mesh point coords
         ip = reset_coordinates ( mesh, Pcoords );
      }
   }

   if ( ip > -1 ) // Either the point is a mesh point or has been reseted to one of the mesh points
   {
      // Find closest time values in mesh->times array and interpolate
      for ( itime = 0; itime < mesh->nTimes && tsearch; itime++ )
      {
         if ( mesh->times[itime] == t ) // Exact velocity value is already stored
         {
            interpolation[0] = mesh->points[ip].velocity[itime*mesh->nDim];
            interpolation[1] = mesh->points[ip].velocity[itime*mesh->nDim+1];
	    if ( mesh->nDim == 3 ) interpolation[2] = mesh->points[ip].velocity[itime*mesh->nDim+2];
            tsearch = 0;
         }
         else
         {
            if ( mesh->times[itime] > t ) // There exist 2 times in provided mesh data surrounding given t
            {
               itprev = itime-1;
               itpost = itime;
               t0  = mesh->times[itprev];
               t1  = mesh->times[itpost];

               v0[0] = mesh->points[ip].velocity[itprev*mesh->nDim];
               v0[1] = mesh->points[ip].velocity[itprev*mesh->nDim+1];
	       if ( mesh->nDim == 3 ) v0[2] = mesh->points[ip].velocity[itprev*mesh->nDim+2];

               v1[0] = mesh->points[ip].velocity[itpost*mesh->nDim];
               v1[1] = mesh->points[ip].velocity[itpost*mesh->nDim+1];
	       if ( mesh->nDim == 3 ) v1[2] = mesh->points[ip].velocity[itpost*mesh->nDim+2];

               interpolation[0] = (v0[0] * (t1-t) + v1[0] * (t-t0))/(t1-t0);
               interpolation[1] = (v0[1] * (t1-t) + v1[1] * (t-t0))/(t1-t0);
	       if ( mesh->nDim == 3 ) interpolation[2] = (v0[2] * (t1-t) + v1[2] * (t-t0))/(t1-t0);

               tsearch = 0;
            }
         }
      }
      if (tsearch) // There is not enough data to interpolate
      {
         fprintf( stderr, "Error: There is not enough information to perform linear_interpolation for t=%f\n", t );
         exit(-1);
      }
   }
}
