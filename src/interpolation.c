#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "interpolation.h"
#include "location.h"
#include "distance.h"
#include "search.h"
#include "kdtree.h"

int reset_coordinates ( struct kdtree *kd, double *Pcoords, int nDim )
{
   struct kdres *nearest;
   double pos[nDim];
   int data;
   nearest = kd_nearest (kd, Pcoords);
   data    = kd_res_item(nearest, pos);
   return data;
}

void interpolate_triangle ( 	double *coords_x, double *coords_y, double *velocities,
				int iv1, int iv2, int iv3, int itime,
				int nPoints, double *P, double *interpolated_vel, int nDim )
{
   double denom, wv1, wv2, wv3;

   denom = (coords_y[iv2]-coords_y[iv3]) * (coords_x[iv1]-coords_x[iv3]) + (coords_x[iv3]-coords_x[iv2]) * (coords_y[iv1]-coords_y[iv3]);
   wv1   = (coords_y[iv2]-coords_y[iv3]) * (P[0]-coords_x[iv3])  + (coords_x[iv3]-coords_x[iv2]) * (P[1]-coords_y[iv3]);
   wv1   = wv1 / denom;
   wv2   = (coords_y[iv3]-coords_y[iv1]) * (P[0]-coords_x[iv3])  + (coords_x[iv1]-coords_x[iv3]) * (P[1]-coords_y[iv3]);
   wv2   = wv2 / denom;
   wv3   = 1 - wv1 - wv2;

   interpolated_vel[0] = wv1 * velocities[itime*nPoints*nDim + iv1*nDim]     + wv2 * velocities[itime*nPoints*nDim + iv2*nDim]     + wv3 * velocities[itime*nPoints*nDim + iv3*nDim];
   interpolated_vel[1] = wv1 * velocities[itime*nPoints*nDim + iv1*nDim + 1] + wv2 * velocities[itime*nPoints*nDim + iv2*nDim + 1] + wv3 * velocities[itime*nPoints*nDim + iv3*nDim + 1];
}

void interpolate_3D_tetrahedral ( double *coords_x, double *coords_y, double *coords_z, double *velocities,
                                  int iv1, int iv2, int iv3, int iv4, int itime,
                                  int nPoints, double *P, double *interpolated_vel, int nDim )
{
	int d;
        double num, denom;
        double lambda[nDim + 1];

        double V1[3], V2[3], V3[3], V4[3];
        V1[0] = coords_x[iv1]; V1[1] = coords_y[iv1]; V1[2] = coords_z[iv1];
        V2[0] = coords_x[iv2]; V2[1] = coords_y[iv2]; V2[2] = coords_z[iv2];
        V3[0] = coords_x[iv3]; V3[1] = coords_y[iv3]; V3[2] = coords_z[iv3];
        V4[0] = coords_x[iv4]; V4[1] = coords_y[iv4]; V4[2] = coords_z[iv4];

        double velV1[3], velV2[3], velV3[3], velV4[3];
        velV1[0] = velocities[itime * nPoints * nDim + iv1 * nDim];
        velV1[1] = velocities[itime * nPoints * nDim + iv1 * nDim + 1];
        velV1[2] = velocities[itime * nPoints * nDim + iv1 * nDim + 2];
        velV2[0] = velocities[itime * nPoints * nDim + iv2 * nDim];
        velV2[1] = velocities[itime * nPoints * nDim + iv2 * nDim + 1];
        velV2[2] = velocities[itime * nPoints * nDim + iv2 * nDim + 2];
        velV3[0] = velocities[itime * nPoints * nDim + iv3 * nDim];
        velV3[1] = velocities[itime * nPoints * nDim + iv3 * nDim + 1];
        velV3[2] = velocities[itime * nPoints * nDim + iv3 * nDim + 2];
        velV4[0] = velocities[itime * nPoints * nDim + iv4 * nDim];
        velV4[1] = velocities[itime * nPoints * nDim + iv4 * nDim + 1];
        velV4[2] = velocities[itime * nPoints * nDim + iv4 * nDim + 2];
	 
    lambda[0] = (-V3[1] * V2[2] + V4[1] * V2[2] + V2[1] * V3[2] - V4[1] * V3[2] - V2[1] * V4[2] + V3[1] * V4[2]) * (P[0] - V4[0]) +
               (V3[0] * V2[2] - V4[0] * V2[2] - V2[0] * V3[2] + V4[0] * V3[2] + V2[0] * V4[2] - V3[0] * V4[2]) * (P[1] - V4[1]) +
               (-V3[0] * V2[1] + V4[0] * V2[1] + V2[0] * V3[1] - V4[0] * V3[1] - V2[0] * V4[1] + V3[0] * V4[1]) * (P[2] - V4[2]);

    lambda[1] = (V3[1] * V1[2] - V4[1] * V1[2] - V1[1] * V3[2] + V4[1] * V3[2] + V1[1] * V4[2] - V3[1] * V4[2]) * (P[0] - V4[0]) +
               (-V3[0] * V1[2] + V4[0] * V1[2] + V1[0] * V3[2] - V4[0] * V3[2] - V1[0] * V4[2] + V3[0] * V4[2]) * (P[1] - V4[1]) +
               (V3[0] * V1[1] - V4[0] * V1[1] - V1[0] * V3[1] + V4[0] * V3[1] + V1[0] * V4[1] - V3[0] * V4[1]) * (P[2] - V4[2]);

    lambda[2] = (-V2[1] * V1[2] + V4[1] * V1[2] + V1[1] * V2[2] - V4[1] * V2[2] - V1[1] * V4[2] + V2[1] * V4[2]) * (P[0] - V4[0]) +
               (V2[0] * V1[2] - V4[0] * V1[2] - V1[0] * V2[2] + V4[0] * V2[2] + V1[0] * V4[2] - V2[0] * V4[2]) * (P[1] - V4[1]) +
               (-V2[0] * V1[1] + V4[0] * V1[1] + V1[0] * V2[1] - V4[0] * V2[1] - V1[0] * V4[1] + V2[0] * V4[1]) * (P[2] - V4[2]);

    double det = (V1[0] - V4[0]) * (V2[1] - V4[1]) * (V3[2] - V4[2]) + (V2[0] - V4[0]) * (V3[1] - V4[1]) * (V1[2] - V4[2]) + (V1[1] - V4[1]) * (V2[2] - V4[2]) * (V3[0] - V4[0]) -
                (V3[0] - V4[0]) * (V2[1] - V4[1]) * (V1[2] - V4[2]) - (V2[2] - V4[2]) * (V3[1] - V4[1]) * (V1[0] - V4[0]) - (V1[1] - V4[1]) * (V2[0] - V4[0]) * (V3[2] - V4[2]);

    lambda[0] = lambda[0]/det;
    lambda[1] = lambda[1]/det;
    lambda[2] = lambda[2]/det;

    lambda[3] = 1 - lambda[0] - lambda[1] - lambda[2];

	for ( d = 0; d < nDim; d++ )
        {
                interpolated_vel[d] +=  lambda[0] * velV1[d] +
                                                                lambda[1] * velV2[d] +
                                                                lambda[2] * velV3[d] +
                                                                lambda[3] * velV4[d];
        }

}

void interpolate_3D_tetrahedral_old ( double *coords_x, double *coords_y, double *coords_z, double *velocities,
				  int iv1, int iv2, int iv3, int iv4, int itime,


				  int nPoints, double *P, double *interpolated_vel, int nDim )
{
	int d;
	double num, denom;
	double lambda[nDim + 1];

	double V1[3], V2[3], V3[3], V4[3];
	V1[0] = coords_x[iv1]; V1[1] = coords_y[iv1]; V1[2] = coords_z[iv1];
	V2[0] = coords_x[iv2]; V2[1] = coords_y[iv2]; V2[2] = coords_z[iv2];
	V3[0] = coords_x[iv3]; V3[1] = coords_y[iv3]; V3[2] = coords_z[iv3];
	V4[0] = coords_x[iv4]; V4[1] = coords_y[iv4]; V4[2] = coords_z[iv4];

	double velV1[3], velV2[3], velV3[3], velV4[3];
	velV1[0] = velocities[itime * nPoints * nDim + iv1 * nDim];
	velV1[1] = velocities[itime * nPoints * nDim + iv1 * nDim + 1];
	velV1[2] = velocities[itime * nPoints * nDim + iv1 * nDim + 2];
	velV2[0] = velocities[itime * nPoints * nDim + iv2 * nDim];
	velV2[1] = velocities[itime * nPoints * nDim + iv2 * nDim + 1];
	velV2[2] = velocities[itime * nPoints * nDim + iv2 * nDim + 2];
	velV3[0] = velocities[itime * nPoints * nDim + iv3 * nDim];
	velV3[1] = velocities[itime * nPoints * nDim + iv3 * nDim + 1];
	velV3[2] = velocities[itime * nPoints * nDim + iv3 * nDim + 2];
	velV4[0] = velocities[itime * nPoints * nDim + iv4 * nDim];
	velV4[1] = velocities[itime * nPoints * nDim + iv4 * nDim + 1];
	velV4[2] = velocities[itime * nPoints * nDim + iv4 * nDim + 2];


	/* 1. Compute lambda values (barycentric coordinates) */
	
	/* 1.1 Compute lambda3 (located at lambda[2]) */
	num = ( V2[0] - V4[0] - V2[1] + V4[1] ) * ( V1[1] - V4[1] ) * ( P[2] + (V1[2] - V4[2])/(V1[1] - V4[1]) * (V4[1] - P[1]) - V4[2] - (P[0]-P[1]+V4[1]-V4[0])/(V2[0]-V4[0]-V2[1]+V4[1]) );
	denom = (V3[1] - V4[1] - V3[0] + V4[0]) * (V1[1]-V4[1]) + (V3[2]*V1[1] - V3[2]*V4[1] - V4[2]*V1[1] - V1[2]*V3[1] + V1[2]*V4[1] + V4[2]*V3[1]) * (V2[0]-V4[0] -V2[1]+V4[1]);
	lambda[2] = num/denom;

	/* 1.2 Compute lambda2 (located at lambda[1]) */
	num = P[0] - P[1] + V4[1] - V4[0] + lambda[2] * ( V3[1] - V4[1] - V3[0] + V4[0] );
	denom = V2[0] - V4[0] - V2[1] + V4[1];
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

void linear_interpolation_approach2_2D ( double t, double *Pcoords, double *times, double *interpolation, int nDim, int nPoints, int nTimes, int nVertsPerFace, int nFaces, int *faces, double *coords_x, double *coords_y, double *velocities, struct kdtree *kd, int *nFacesPerPoint, int *facesPerPoint)
{
   int ip, iface, itime, itprev, itpost, tsearch = 1;
   int nFacesNearestPoint, ifp, fsearch, faceNearestPoint;
   int iv1, iv2, iv3;

   double t0, t1;
   double v0[nDim];
   double v1[nDim];
   double interp1[nDim];
   double interp2[nDim];
   double interp3[nDim];

   double pos[2];
   int data;

   struct kdres *nearest;

   tsearch = binarySearch(times, nTimes, t, &itime);
   nearest = kd_nearest (kd, Pcoords);
   data    = kd_res_item(nearest, pos);

   if ( pos[0] != Pcoords[0] || pos[1] != Pcoords[1] )
   {
      ip = -1;
   }
   else
   {
      ip = data;
   }

   if ( ip == -1 ) // P is not a mesh point
   {
      // TODO (future): Consider closest point might not belong to the triangle containing the given point
      iface = -1;
      data > 0 ? nFacesNearestPoint = nFacesPerPoint[data] - nFacesPerPoint[data-1] : nFacesPerPoint[0];
      fsearch = 1;
      for ( ifp = 0; ifp < nFacesNearestPoint && fsearch; ifp++ )
      {
         data > 0 ? faceNearestPoint = facesPerPoint[nFacesPerPoint[data-1]+ifp] : facesPerPoint[ifp];
         iv1 = faces[faceNearestPoint*nVertsPerFace];
         iv2 = faces[faceNearestPoint*nVertsPerFace+1];
         iv3 = faces[faceNearestPoint*nVertsPerFace+2];
	 if ( p_inside_triangle ( iv1, iv2, iv3, coords_x, coords_y, Pcoords, nDim ) )
	 {
	    iface = faceNearestPoint;
            fsearch = 0;
         }
      } 

      if (iface > -1) // There is a mesh face or volume containing the given point
      {
         // Find closest time values in mesh->times array and interpolate
	 if (tsearch) // times[itime] == t // t exists in the mesh times elements
         {
	 	iv1 = faces[iface*nVertsPerFace];
		iv2 = faces[iface*nVertsPerFace+1];
		iv3 = faces[iface*nVertsPerFace+2];
			
		interpolate_triangle (  coords_x, coords_y, velocities, 
					iv1, iv2, iv3, itime, 
					nPoints, Pcoords, interpolation, nDim );
            }
            else
            {
               if ( tsearch == 0 ) // times[itime] > t // There exist 2 times in provided mesh data surrounding given t
               {
	                 itprev = itime-1;
                	 itpost = itime;
        	         t0  = times[itprev];
	                 t1  = times[itpost];

			iv1 = faces[iface*nVertsPerFace];
                       	iv2 = faces[iface*nVertsPerFace+1];
        		iv3 = faces[iface*nVertsPerFace+2];

			interpolate_triangle ( coords_x, coords_y, velocities,
                                        iv1, iv2, iv3, itprev,
                                        nPoints, Pcoords, interp1, nDim );

			interpolate_triangle ( coords_x, coords_y, velocities,
                                        iv1, iv2, iv3, itpost,
                                        nPoints, Pcoords, interp2, nDim );

		         interpolation[0] = (interp1[0] * (t1-t) + interp2[0] * (t-t0))/(t1-t0);
        		 interpolation[1] = (interp1[1] * (t1-t) + interp2[1] * (t-t0))/(t1-t0);
               }
            }
            if (tsearch < 0) // There is not enough data to interpolate
            {
            	fprintf( stderr, "Error: There is not enough information to perform linear_interpolation for t=%f\n", t );
            	exit(-1);
            }
      }
      else // The point coords are outside the given mesh data
      {
         // Pcoords will be reset to the closest mesh point coords
         ip = reset_coordinates ( kd, Pcoords, nDim );
      }
   }
   if ( ip > -1 ) // Either the point is a mesh point or has been reseted to one of the mesh points
   {
      // Find closest time values in mesh->times array and interpolate
         if (tsearch) // times[itime] == t // t exists in the mesh times elements
         {
            interpolation[0] = velocities[itime*nPoints*nDim + ip*nDim];//mesh->points[ip].velocity[itime*mesh->nDim];
            interpolation[1] = velocities[itime*nPoints*nDim + ip*nDim + 1];//mesh->points[ip].velocity[itime*mesh->nDim+1];
         }
         else
         {
            if ( tsearch == 0 ) // times[itime] > t // There exist 2 times in provided mesh data surrounding given t
            {
               itprev = itime-1;
               itpost = itime;
               t0  = times[itprev];
               t1  = times[itpost];

               v0[0] = velocities[itprev*nPoints*nDim + ip*nDim];//   mesh->points[ip].velocity[itprev*mesh->nDim];
               v0[1] = velocities[itprev*nPoints*nDim + ip*nDim + 1];// mesh->points[ip].velocity[itprev*mesh->nDim+1];

               v1[0] = velocities[itpost*nPoints*nDim + ip*nDim]; //mesh->points[ip].velocity[itpost*mesh->nDim];
               v1[1] = velocities[itpost*nPoints*nDim + ip*nDim + 1];//mesh->points[ip].velocity[itpost*mesh->nDim+1];

               interpolation[0] = (v0[0] * (t1-t) + v1[0] * (t-t0))/(t1-t0);
               interpolation[1] = (v0[1] * (t1-t) + v1[1] * (t-t0))/(t1-t0);
            }
         }
      if (tsearch < 0) // There is not enough data to interpolate
      {
         fprintf( stderr, "Error: There is not enough information to perform linear_interpolation for t=%f\n", t );
         exit(-1);
      }
   }
}

void linear_interpolation_approach2_3D ( double t, double *Pcoords, double *times, double *interpolation, int nDim, int nPoints, int nTimes, int nVertsPerFace, int nFaces, int *faces, double *coords_x, double *coords_y, double *coords_z, double*velocities, struct kdtree *kd, int *nFacesPerPoint, int *facesPerPoint)
{
   int ip, iface, itime, itprev, itpost, tsearch = 1;
   int nFacesNearestPoint, ifp, fsearch, faceNearestPoint;
   int iv1, iv2, iv3, iv4;

   double t0, t1;
   double v0[nDim];
   double v1[nDim];
   double interp1[nDim];
   double interp2[nDim];
   double interp3[nDim];

   double pos[3];
   int data;

   struct kdres *nearest;

   tsearch = binarySearch(times, nTimes, t, &itime);

   nearest = kd_nearest (kd, Pcoords);
   data = kd_res_item(nearest, pos);

   if ( pos[0] != Pcoords[0] || pos[1] != Pcoords[1] || pos[2] != Pcoords[2] )
   {  
      ip = -1;
   }
   else
   {
      ip = data;
   }

   if ( ip == -1 ) // P is not a mesh point
   {
      // TODO: Consider closest point might not belong to the triangle containing the given point
      iface = -1;
      data > 0 ? nFacesNearestPoint = nFacesPerPoint[data] - nFacesPerPoint[data-1] : nFacesPerPoint[0];
      fsearch = 1;
      for ( ifp = 0; ifp < nFacesNearestPoint && fsearch; ifp++ )
      {
         data > 0 ? faceNearestPoint = facesPerPoint[nFacesPerPoint[data-1]+ifp] : facesPerPoint[ifp];
         iv1 = faces[faceNearestPoint*nVertsPerFace];
         iv2 = faces[faceNearestPoint*nVertsPerFace+1];
         iv3 = faces[faceNearestPoint*nVertsPerFace+2];
         iv4 = faces[faceNearestPoint*nVertsPerFace+3];
         if ( p_inside_tetrahedron ( iv1, iv2, iv3, iv4, coords_x, coords_y, coords_z, Pcoords, nDim ) )
         {
            iface = faceNearestPoint;
            fsearch = 0;
         }
      }

      if (iface > -1) // There is a mesh face or volume containing given point
      {
         // Find closest time values in mesh->times array and interpolate
            if (tsearch) // times[itime] == t // t exists in the mesh times elements
            {
                        iv1 = faces[iface*nVertsPerFace];
                        iv2 = faces[iface*nVertsPerFace+1];
                        iv3 = faces[iface*nVertsPerFace+2];
                        iv4 = faces[iface*nVertsPerFace+3];
			interpolate_3D_tetrahedral ( coords_x, coords_y, coords_z, velocities,
						     iv1, iv2, iv3, iv4, itime,
						     nPoints, Pcoords, interpolation, nDim );
            }
            else
            {
               if ( tsearch == 0 ) // times[itime] > t // There exist 2 times in provided mesh data surrounding given t
               {
                 itprev = itime-1;
                 itpost = itime;
                 t0  = times[itprev];
                 t1  = times[itpost];

                        iv1 = faces[iface*nVertsPerFace];
                        iv2 = faces[iface*nVertsPerFace+1];
                        iv3 = faces[iface*nVertsPerFace+2];
                        iv4 = faces[iface*nVertsPerFace+3];

			interpolate_3D_tetrahedral ( coords_x, coords_y, coords_z, velocities,
                                                     iv1, iv2, iv3, iv4, itprev, 
                                                     nPoints, Pcoords, interp1, nDim );

			interpolate_3D_tetrahedral ( coords_x, coords_y, coords_z, velocities,
                                                     iv1, iv2, iv3, iv4, itpost,
                                                     nPoints, Pcoords, interp2, nDim );
		
                 interpolation[0] = (interp2[0] * (t1-t) + interp1[0] * (t-t0))/(t1-t0);
                 interpolation[1] = (interp2[1] * (t1-t) + interp1[1] * (t-t0))/(t1-t0);
                 interpolation[2] = (interp2[2] * (t1-t) + interp1[2] * (t-t0))/(t1-t0);

               }
            }
         if (tsearch < 0) // There is not enough data to interpolate
         {
            fprintf( stderr, "Error: There is not enough information to perform linear_interpolation for t=%f\n", t );
            exit(-1);
         }
      }
      else // The point coords are outside the given mesh data
      {
         // Pcoords will be reset to the closest mesh point coords
         ip = reset_coordinates ( kd, Pcoords, nDim );
      }
   }
   if ( ip > -1 ) // Either the point is a mesh point or has been reseted to one of the mesh points
   {
      // Find closest time values in mesh->times array and interpolate
         if (tsearch) // times[itime] == t // t exists in the mesh times elements
         {
            interpolation[0] = velocities[itime*nPoints*nDim + ip*nDim];//mesh->points[ip].velocity[itime*mesh->nDim];
            interpolation[1] = velocities[itime*nPoints*nDim + ip*nDim + 1];//mesh->points[ip].velocity[itime*mesh->nDim+1];
            interpolation[2] = velocities[itime*nPoints*nDim + ip*nDim + 2];//mesh->points[ip].velocity[itime*mesh->nDim+2];
            }
	else
	{
		if ( tsearch == 0 ) // times[itime] > t // There exist 2 times in provided mesh data surrounding given t
            {
               itprev = itime-1;
               itpost = itime;
               t0  = times[itprev];
               t1  = times[itpost];

               v0[0] = velocities[itprev*nPoints*nDim + ip*nDim];//   mesh->points[ip].velocity[itprev*mesh->nDim];
               v0[1] = velocities[itprev*nPoints*nDim + ip*nDim + 1];// mesh->points[ip].velocity[itprev*mesh->nDim+1];
               v0[2] = velocities[itprev*nPoints*nDim + ip*nDim + 2]; //mesh->points[ip].velocity[itprev*mesh->nDim+2];

               v1[0] = velocities[itpost*nPoints*nDim + ip*nDim]; //mesh->points[ip].velocity[itpost*mesh->nDim];
               v1[1] = velocities[itpost*nPoints*nDim + ip*nDim + 1];//mesh->points[ip].velocity[itpost*mesh->nDim+1];
               v1[2] = velocities[itpost*nPoints*nDim + ip*nDim + 2];//mesh->points[ip].velocity[itpost*mesh->nDim+2];

               interpolation[0] = (v0[0] * (t1-t) + v1[0] * (t-t0))/(t1-t0);
               interpolation[1] = (v0[1] * (t1-t) + v1[1] * (t-t0))/(t1-t0);
               interpolation[2] = (v0[2] * (t1-t) + v1[2] * (t-t0))/(t1-t0);

         }
	}
      if (tsearch < 0) // There is not enough data to interpolate
      {
         fprintf( stderr, "Error: There is not enough information to perform linear_interpolation for t=%f\n", t );
         exit(-1);
      }
   }
}
