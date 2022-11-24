#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "location.h"
#include "distance.h"

int find_point_simplex_in_mesh ( double *Pcoords, int nDim, int nFaces, int nVertsPerFace, double *coords, int *faces )
{
   int res = -1;
   int iface, iv1, iv2, iv3, iv4, search_face = 1;
   for ( iface = 0; iface < nFaces && search_face; iface++ )
   {
      iv1 = faces[iface*nVertsPerFace];
      iv2 = faces[iface*nVertsPerFace+1];
      iv3 = faces[iface*nVertsPerFace+2];
      iv4 = faces[iface*nVertsPerFace+3];
      if ( (nDim == 2) && (p_inside_triangle(&coords[iv1*nDim], &coords[iv2*nDim], &coords[iv3*nDim], Pcoords, nDim)) )
      {
         search_face = 0;
	 res = iface;
      }
      else if ( (nDim == 3) && (p_inside_tetrahedron(&coords[iv1*nDim], &coords[iv2*nDim], &coords[iv3*nDim], &coords[iv4*nDim], Pcoords, nDim)) )
      {
         search_face = 0;
         res = iface;
      }
   }
   return res;
}

int find_point_in_mesh ( double *Pcoords, int nDim, int nPoints, double *coords )
{
   int i, pindex = -1;
   for ( i = 0; i < nPoints && pindex < 0; i++ )
   {
      if ( (coords[i*nDim] == Pcoords[0]) && (coords[i*nDim+1] == Pcoords[1]) )
      {
         if ( (nDim == 2) || ( (nDim == 3) && (coords[i*nDim+2] == Pcoords[2]) ) )
         {
            pindex = i;
         }
      }
   }
   return pindex;
}

int p_inside_triangle ( double *V1, double *V2, double *V3, double *P, int nDim )
{
   double v0[nDim];
   double v1[nDim];
   double v2[nDim];
   double dot00, dot01, dot02, dot11, dot12, invDenom, u, v;

   // Compute vectors v0, v1, v2
   v0[0] = V3[0] - V1[0];
   v0[1] = V3[1] - V1[1];
   v1[0] = V2[0] - V1[0];
   v1[1] = V2[1] - V1[1];
   v2[0] = P[0]  - V1[0];
   v2[1] = P[1]  - V1[1];

   // Compute dot products
   dot00 = v0[0]*v0[0] + v0[1]*v0[1];
   dot01 = v0[0]*v1[0] + v0[1]*v1[1];
   dot02 = v0[0]*v2[0] + v0[1]*v2[1];
   dot11 = v1[0]*v1[0] + v1[1]*v1[1];
   dot12 = v1[0]*v2[0] + v1[1]*v2[1];
     
   // Compute barycentric coordinates
   invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
   u = (dot11 * dot02 - dot01 * dot12) * invDenom;
   v = (dot00 * dot12 - dot01 * dot02) * invDenom;

   // Check if point is in triangle
   return (u >= 0) && (v >= 0) && (u + v < 1);
}

int share_side ( double *V1, double *V2, double *V3, double *V4, double *P, int nDim )
{
   int i;
   double cross[nDim];
   double V2V1[nDim];
   double V3V1[nDim];
   double V4V1[nDim];
   double PV1[nDim];
   double dot = 0.0, dotP = 0.0; 

   for ( i = 0; i < nDim; i++ )
   {
      V2V1[i] = V2[i] - V1[i];
      V3V1[i] = V3[i] - V1[i];
      V4V1[i] = V4[i] - V1[i];
      PV1[i]  = P[i]  - V1[i];
   }

   cross[0] = V2V1[1] * V3V1[2] - V2V1[2] * V3V1[1];
   cross[1] = V2V1[2] * V3V1[0] - V2V1[0] * V3V1[2];
   cross[2] = V2V1[0] * V3V1[1] - V2V1[1] * V3V1[0];

   for ( i = 0; i < nDim; i++ )
   {
      dot  += cross[i] * V4V1[i];
      dotP += cross[i] * PV1[i];
   }

   return ( dot == 0 && dotP == 0 ) || 
          ( dot <  0 && dotP <  0 ) || 
          ( dot >  0 && dotP >  0 );

}

int p_inside_tetrahedron ( double *V1, double *V2, double *V3, double *V4, double *P, int nDim )
{
   
   int share_1234P = share_side ( V1, V2, V3, V4, P, nDim );
   int share_2341P = share_side ( V2, V3, V4, V1, P, nDim );
   int share_3412P = share_side ( V3, V4, V1, V2, P, nDim );
   int share_4123P = share_side ( V4, V1, V2, V3, P, nDim );
   
   return ( share_1234P && share_2341P && share_3412P && share_4123P );
}
