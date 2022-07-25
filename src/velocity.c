#include <stdio.h>
#include <stdlib.h>
#include "velocity.h"

void compute_velocity ( double t, double *coords, mesh_t *mesh, double *velocity )
{
   int i;
   for ( i = 0; i < mesh->nDim; i++ )
      velocity[i] = coords[i] + t;
}

void compute_velocity_vector ( mesh_t *mesh )
{
   int it, ip;
   for ( ip = 0; ip < mesh->nPoints; ip++ )
   {
      mesh->points[ip].velocity = malloc( sizeof(double) * mesh->nDim * mesh->nTimes );
      if ( mesh->points[ip].velocity == NULL )
      {
         fprintf( stderr, "Error: Allocating memory for velocity in compute_velocity_vector function\n" );
         exit(-1);
      }
      for ( it = 0; it < mesh->nTimes; it++ )
         compute_velocity ( mesh->times[it], mesh->points[ip].coordinates, mesh, &(mesh->points[ip].velocity[it*mesh->nDim]) );
   }
}