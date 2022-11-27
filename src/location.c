#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "location.h"
#include "distance.h"

int find_point_simplex_in_mesh_2D ( double *Pcoords, int nDim, int nFaces, int nVertsPerFace, double *coords_x, double *coords_y, int *faces )
{
   int res = -1;
   int iface, iv1, iv2, iv3, search_face = 1;
   for ( iface = 0; iface < nFaces && search_face; iface++ )
   {
      iv1 = faces[iface*nVertsPerFace];
      iv2 = faces[iface*nVertsPerFace+1];
      iv3 = faces[iface*nVertsPerFace+2];
      if ( p_inside_triangle(iv1, iv2, iv3, coords_x, coords_y, Pcoords, nDim) )
      {
         search_face = 0;
	 res = iface;
      }
   }
   return res;
}

int find_point_simplex_in_mesh_3D ( double *Pcoords, int nDim, int nFaces, int nVertsPerFace, double *coords_x, double *coords_y, double *coords_z, int *faces )
{
   int res = -1;
   int iface, iv1, iv2, iv3, iv4, search_face = 1;
   for ( iface = 0; iface < nFaces && search_face; iface++ )
   {
      iv1 = faces[iface*nVertsPerFace];
      iv2 = faces[iface*nVertsPerFace+1];
      iv3 = faces[iface*nVertsPerFace+2];
      iv4 = faces[iface*nVertsPerFace+3];
      if ( p_inside_tetrahedron(iv1, iv2, iv3, iv4, coords_x, coords_y, coords_z, Pcoords, nDim) )
      {
         search_face = 0;
         res = iface;
      }
   }
   return res;
}

int find_point_in_mesh_2D ( double *Pcoords, int nDim, int nPoints, double *coords_x, double *coords_y )
{
   int i, pindex = -1;
   for ( i = 0; i < nPoints && pindex < 0; i++ )
   {
      if ( (coords_x[i] == Pcoords[0]) && (coords_y[i] == Pcoords[1]) )
      {
            pindex = i;
      }
   }
   return pindex;
}

int find_point_in_mesh_3D ( double *Pcoords, int nDim, int nPoints, double *coords_x, double *coords_y, double *coords_z )
{
   int i, pindex = -1;
   for ( i = 0; i < nPoints && pindex < 0; i++ )
   {
      if ( (coords_x[i] == Pcoords[0]) && (coords_y[i] == Pcoords[1]) && (coords_z[i] == Pcoords[2]) )
      {
            pindex = i;
      }
   }
   return pindex;
}

int p_inside_triangle ( int iv1, int iv2, int iv3, double *coords_x, double *coords_y, double *P, int nDim )
{
   double v0[nDim];
   double v1[nDim];
   double v2[nDim];
   double dot00, dot01, dot02, dot11, dot12, invDenom, u, v;

   // Compute vectors v0, v1, v2
   v0[0] = coords_x[iv3] - coords_x[iv1];
   v0[1] = coords_y[iv3] - coords_y[iv1];
   v1[0] = coords_x[iv2] - coords_x[iv1];
   v1[1] = coords_y[iv2] - coords_y[iv1];
   v2[0] = P[0]          - coords_x[iv1];
   v2[1] = P[1]          - coords_y[iv1];

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

int p_inside_tetrahedron ( int iv1, int iv2, int iv3, int iv4, double *coords_x, double *coords_y, double *coords_z, double *P, int nDim )
{
   double V1[3], V2[3], V3[3], V4[3];
   V1[0] = coords_x[iv1]; V1[1] = coords_y[iv1]; V1[2] = coords_z[iv1];
   V2[0] = coords_x[iv2]; V2[1] = coords_y[iv2]; V2[2] = coords_z[iv2];
   V3[0] = coords_x[iv3]; V3[1] = coords_y[iv3]; V3[2] = coords_z[iv3];
   V4[0] = coords_x[iv4]; V4[1] = coords_y[iv4]; V4[2] = coords_z[iv4];

   int share_1234P = share_side ( V1, V2, V3, V4, P, nDim );
   int share_2341P = share_side ( V2, V3, V4, V1, P, nDim );
   int share_3412P = share_side ( V3, V4, V1, V2, P, nDim );
   int share_4123P = share_side ( V4, V1, V2, V3, P, nDim );
   
   return ( share_1234P && share_2341P && share_3412P && share_4123P );
}
