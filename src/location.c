#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "location.h"
#include "distance.h"

int find_point_simplex_in_mesh ( double *Pcoords, mesh_t *mesh )
{
   int res = -1;
   int iface, iv1, iv2, iv3, iv4, search_face = 1;
   for ( iface = 0; iface < mesh->nFaces && search_face; iface++ )
   {
      iv1 = mesh->faces[iface*mesh->nVertsPerFace];
      iv2 = mesh->faces[iface*mesh->nVertsPerFace+1];
      iv3 = mesh->faces[iface*mesh->nVertsPerFace+2];
      iv4 = mesh->faces[iface*mesh->nVertsPerFace+3];
      if ( (mesh->nDim == 2) && (p_inside_triangle(mesh, mesh->points[iv1].coordinates, mesh->points[iv2].coordinates, mesh->points[iv3].coordinates, Pcoords)) )
      {
         search_face = 0;
	 res = iface;
      }
      else if ( (mesh->nDim == 3) && (p_inside_tetrahedron(mesh, mesh->points[iv1].coordinates, mesh->points[iv2].coordinates, mesh->points[iv3].coordinates, mesh->points[iv4].coordinates, Pcoords)) )
      {
         search_face = 0;
         res = iface;
      }
   }
   return res;
}

int find_point_in_mesh ( double *Pcoords, mesh_t *mesh )
{
   int i, pindex = -1;
   for ( i = 0; i < mesh->nPoints && pindex < 0; i++ )
   {
      if ( (mesh->points[i].coordinates[0] == Pcoords[0]) && (mesh->points[i].coordinates[1] == Pcoords[1]) )
      {
         if ( (mesh->nDim == 2) || ( (mesh->nDim == 3) && (mesh->points[i].coordinates[2] == Pcoords[2]) ) )
         {
            pindex = i;
         }
      }
   }
   return pindex;
}

int p_inside_triangle ( mesh_t *mesh, double *V1, double *V2, double *V3, double *P )
{
   double v0[mesh->nDim];
   double v1[mesh->nDim];
   double v2[mesh->nDim];
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

int share_side ( mesh_t *mesh, double *V1, double *V2, double *V3, double *V4, double *P )
{
   int i;
   double cross[mesh->nDim];
   double V2V1[mesh->nDim];
   double V3V1[mesh->nDim];
   double V4V1[mesh->nDim];
   double PV1[mesh->nDim];
   double dot = 0.0, dotP = 0.0; 

   for ( i = 0; i < mesh->nDim; i++ )
   {
      V2V1[i] = V2[i] - V1[i];
      V3V1[i] = V3[i] - V1[i];
      V4V1[i] = V4[i] - V1[i];
      PV1[i]  = P[i]  - V1[i];
   }

   cross[0] = V2V1[1] * V3V1[2] - V2V1[2] * V3V1[1];
   cross[1] = V2V1[2] * V3V1[0] - V2V1[0] * V3V1[2];
   cross[2] = V2V1[0] * V3V1[1] - V2V1[1] * V3V1[0];

   for ( i = 0; i < mesh->nDim; i++ )
   {
      dot  += cross[i] * V4V1[i];
      dotP += cross[i] * PV1[i];
   }

   return ( dot == 0 && dotP == 0 ) || 
          ( dot <  0 && dotP <  0 ) || 
          ( dot >  0 && dotP >  0 );

}

int p_inside_tetrahedron ( mesh_t *mesh, double *V1, double *V2, double *V3, double *V4, double *P )
{
   
   int share_1234P = share_side ( mesh, V1, V2, V3, V4, P );
   int share_2341P = share_side ( mesh, V2, V3, V4, V1, P );
   int share_3412P = share_side ( mesh, V3, V4, V1, V2, P );
   int share_4123P = share_side ( mesh, V4, V1, V2, V3, P );
   
   return ( share_1234P && share_2341P && share_3412P && share_4123P );
}
