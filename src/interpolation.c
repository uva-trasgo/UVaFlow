#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "interpolation.h"
#include "location.h"
#include "distance.h"

int reset_coordinates ( int nDim, int nPoints, double *coords, double *Pcoords )
{
   double distance, min_distance;
   int ip, ipmin;
   min_distance = distance_PQ ( nDim, &coords[0], Pcoords );
   ipmin = 0;
   for ( ip = 1; ip < nPoints; ip++ )
   {
      distance = distance_PQ ( nDim, &coords[ip*nDim], Pcoords );
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
                              double *P,
                              double *interpolated_vel,
			      int nDim )
{
   double denom, wv1, wv2, wv3;

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
                      double *P,
                      double *interpolated_vel,
		      int nDim )
{
	int d;
	double num, denom;
	double lambda[nDim + 1];

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
	for ( d = 0; d < nDim; d++ )
	{
		interpolated_vel[d] += 	lambda[0] * velV1[d] +  
								lambda[1] * velV2[d] +
								lambda[2] * velV3[d] +
								lambda[3] * velV4[d];
	}
}

void linear_interpolation_approach2 ( double t, double *Pcoords, double *times, double *interpolation, int nDim, int nPoints, int nTimes, int nVertsPerFace, int nFaces, int *faces, double *coords, double*velocities)
{
   int ip, iface, itime, itprev, itpost, tsearch = 1;
   int iv1, iv2, iv3, iv4;
   double t0, t1;
   double v0[nDim];
   double v1[nDim];
   double interp1[nDim];
   double interp2[nDim];
   double interp3[nDim];

   ip = find_point_in_mesh ( Pcoords, nDim, nPoints, coords );
   //printf("Point in mesh\n");
   if ( ip == -1 ) // P is not a mesh point
   {
      iface = find_point_simplex_in_mesh ( Pcoords, nDim, nFaces, nVertsPerFace, coords, faces);
      if (iface > -1) // There is a mesh face or volume containing given point
      {
         // Find closest time values in mesh->times array and interpolate
	 for ( itime = 0; itime < nTimes && tsearch; itime++ )
	 {
	    if ( times[itime] == t ) // t exists in the mesh times elements
            {
		if ( nDim == 2 )
		{
			iv1 = faces[iface*nVertsPerFace];
			iv2 = faces[iface*nVertsPerFace+1];
			iv3 = faces[iface*nVertsPerFace+2];
			interpolate_triangle ( &coords[iv1 * nDim], &velocities[itime * nPoints * nDim + iv1 * nDim],
					       &coords[iv2 * nDim], &velocities[itime * nPoints * nDim + iv2 * nDim],
					       &coords[iv3 * nDim], &velocities[itime * nPoints * nDim + iv3 * nDim],
                                               Pcoords, interpolation, nDim );
			/*
                	interpolate_triangle (  mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].coordinates,  &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].velocity[mesh->nDim*itime]),
                                        mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].velocity[mesh->nDim*itime]),
                                        mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].velocity[mesh->nDim*itime]),
                                        Pcoords, mesh, interpolation );
			*/
		} 
		else // mesh->nDim == 3
		{
			iv1 = faces[iface*nVertsPerFace];
                        iv2 = faces[iface*nVertsPerFace+1];
                        iv3 = faces[iface*nVertsPerFace+2];
			iv4 = faces[iface*nVertsPerFace+3];
                        interpolate_3D_tetrahedral ( &coords[iv1 * nDim], &velocities[itime * nPoints * nDim + iv1 * nDim],
                                               &coords[iv2 * nDim], &velocities[itime * nPoints * nDim + iv2 * nDim],
                                               &coords[iv3 * nDim], &velocities[itime * nPoints * nDim + iv3 * nDim],
					       &coords[iv4 * nDim], &velocities[itime * nPoints * nDim + iv4 * nDim],
                                               Pcoords, interpolation, nDim );
			/*
			interpolate_3D_tetrahedral ( mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].coordinates,  &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].velocity[mesh->nDim*itime]),
                                        mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].velocity[mesh->nDim*itime]),
                                        mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].velocity[mesh->nDim*itime]),
					mesh->points[mesh->faces[iface*mesh->nVertsPerFace+3]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+3]].velocity[mesh->nDim*itime]),
                                        Pcoords, mesh, interpolation );
			*/
		}
		tsearch = 0;
            }
            else
            {
               if ( times[itime] > t ) // There exist 2 times in provided mesh data surrounding given t
               {
                 itprev = itime-1;
                 itpost = itime;
                 t0  = times[itprev];
                 t1  = times[itpost];

		 if ( nDim == 2 )
                 {
			iv1 = faces[iface*nVertsPerFace];
                        iv2 = faces[iface*nVertsPerFace+1];
                        iv3 = faces[iface*nVertsPerFace+2];
                        interpolate_triangle ( &coords[iv1 * nDim], &velocities[itprev * nPoints * nDim + iv1 * nDim],
                                               &coords[iv2 * nDim], &velocities[itprev * nPoints * nDim + iv2 * nDim],
                                               &coords[iv3 * nDim], &velocities[itprev * nPoints * nDim + iv3 * nDim],
                                               Pcoords, interp1, nDim );
			interpolate_triangle ( &coords[iv1 * nDim], &velocities[itpost * nPoints * nDim + iv1 * nDim],
                                               &coords[iv2 * nDim], &velocities[itpost * nPoints * nDim + iv2 * nDim],
                                               &coords[iv3 * nDim], &velocities[itpost * nPoints * nDim + iv3 * nDim],
                                               Pcoords, interp2, nDim );
			/*
			 interpolate_triangle (  mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].coordinates,   &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].velocity[mesh->nDim*itprev]),
                                 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].velocity[mesh->nDim*itprev]),
                                 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].velocity[mesh->nDim*itprev]),
                                 Pcoords, mesh, interp1 );

		         interpolate_triangle (  mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].coordinates,   &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace]].velocity[mesh->nDim*itpost]),
                                 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+1]].velocity[mesh->nDim*itpost]),
                                 mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].coordinates, &(mesh->points[mesh->faces[iface*mesh->nVertsPerFace+2]].velocity[mesh->nDim*itpost]),
                                 Pcoords, mesh, interp2 );
			*/
		 }
		 else // mesh->nDim == 3
		 {
			iv1 = faces[iface*nVertsPerFace];
                        iv2 = faces[iface*nVertsPerFace+1];
                        iv3 = faces[iface*nVertsPerFace+2];
                        iv4 = faces[iface*nVertsPerFace+3];
                        interpolate_3D_tetrahedral ( &coords[iv1 * nDim], &velocities[itprev * nPoints * nDim + iv1 * nDim],
                                               &coords[iv2 * nDim], &velocities[itprev * nPoints * nDim + iv2 * nDim],
                                               &coords[iv3 * nDim], &velocities[itprev * nPoints * nDim + iv3 * nDim],
                                               &coords[iv4 * nDim], &velocities[itprev * nPoints * nDim + iv4 * nDim],
                                               Pcoords, interp1, nDim );
                        interpolate_3D_tetrahedral ( &coords[iv1 * nDim], &velocities[itpost * nPoints * nDim + iv1 * nDim],
                                               &coords[iv2 * nDim], &velocities[itpost * nPoints * nDim + iv2 * nDim],
                                               &coords[iv3 * nDim], &velocities[itpost * nPoints * nDim + iv3 * nDim],
                                               &coords[iv4 * nDim], &velocities[itpost * nPoints * nDim + iv4 * nDim],
                                               Pcoords, interp2, nDim );
			/*
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
			*/
		 }

	         interpolation[0] = (interp1[0] * (t1-t) + interp2[0] * (t-t0))/(t1-t0);
        	 interpolation[1] = (interp1[1] * (t1-t) + interp2[1] * (t-t0))/(t1-t0);
		 if ( nDim == 3 ) interpolation[2] = (interp1[2] * (t1-t) + interp2[2] * (t-t0))/(t1-t0);

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
         ip = reset_coordinates ( nDim, nPoints, coords, Pcoords );
      }
   }
   if ( ip > -1 ) // Either the point is a mesh point or has been reseted to one of the mesh points
   {
      // Find closest time values in mesh->times array and interpolate
      for ( itime = 0; itime < nTimes && tsearch; itime++ )
      {
         if ( times[itime] == t )
         {
		//printf("t equal\n");
            interpolation[0] = velocities[itime*nPoints*nDim + ip*nDim];//mesh->points[ip].velocity[itime*mesh->nDim];
            interpolation[1] = velocities[itime*nPoints*nDim + ip*nDim + 1];//mesh->points[ip].velocity[itime*mesh->nDim+1];
	    //printf("interp OK\n");
	    if ( nDim == 3 ) interpolation[2] = velocities[itime*nPoints*nDim + ip*nDim + 2];//mesh->points[ip].velocity[itime*mesh->nDim+2];
            tsearch = 0;
         }
         else
         {
            if ( times[itime] > t ) // There exist 2 times in provided mesh data surrounding given t
            {
		//printf("times surr\n");
               itprev = itime-1;
               itpost = itime;
               t0  = times[itprev];
               t1  = times[itpost];

               v0[0] = velocities[itprev*nPoints*nDim + ip*nDim];//   mesh->points[ip].velocity[itprev*mesh->nDim];
               v0[1] = velocities[itprev*nPoints*nDim + ip*nDim + 1];// mesh->points[ip].velocity[itprev*mesh->nDim+1];
	       if ( nDim == 3 ) v0[2] = velocities[itprev*nPoints*nDim + ip*nDim + 2]; //mesh->points[ip].velocity[itprev*mesh->nDim+2];

               v1[0] = velocities[itpost*nPoints*nDim + ip*nDim]; //mesh->points[ip].velocity[itpost*mesh->nDim];
               v1[1] = velocities[itpost*nPoints*nDim + ip*nDim + 1];//mesh->points[ip].velocity[itpost*mesh->nDim+1];
	       if ( nDim == 3 ) v1[2] = velocities[itpost*nPoints*nDim + ip*nDim + 2];//mesh->points[ip].velocity[itpost*mesh->nDim+2];

               interpolation[0] = (v0[0] * (t1-t) + v1[0] * (t-t0))/(t1-t0);
               interpolation[1] = (v0[1] * (t1-t) + v1[1] * (t-t0))/(t1-t0);
               if ( nDim == 3 ) interpolation[2] = (v0[2] * (t1-t) + v1[2] * (t-t0))/(t1-t0);

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
