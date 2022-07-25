#include <stdio.h>
#include <stdlib.h>
#include "read_python_files.h"

void read_coordinates ( char *filename, mesh_t *mesh )
{
   int ip, d, check_EOF;
   //char buffer[255];
   double buffer;
   int buffer_size;
   FILE *file;

   // Open file
   file = fopen( filename, "r" );

   // First element must be nPoints
   //check_EOF = fscanf(file, "%s", buffer);
   check_EOF = fscanf(file, "%d\n", &buffer_size);
   if ( check_EOF == EOF )
   {
      fprintf( stderr, "Error: Unexpected EOF in read_coordinates\n" );
      exit(-1);
   }

   mesh->nPoints = buffer_size;// atoi(buffer);
   mesh->points  = malloc (sizeof(point_t) * buffer_size);//atoi(buffer));
   if ( mesh->points == NULL )
   {
      fprintf( stderr, "Error: Allocating memory for mesh->points in read_coordinates function\n" );
      exit(-1);
   }
   for ( ip = 0; ip < mesh->nPoints; ip++ )
   {
      mesh->points[ip].index = ip;
      mesh->points[ip].coordinates = malloc(sizeof(double)*mesh->nDim);
      if ( mesh->points[ip].coordinates == NULL )
      {
         fprintf( stderr, "Error: Allocating memory for mesh->points[ip].coordinates in read_coordinates function\n" );
         exit(-1);
      }
   }

   // Rest of read elements will be points' coordinates
   for ( ip = 0; ip < mesh->nPoints; ip++ )
   {
      for ( d = 0; d < mesh->nDim; d++ )
      {
         //check_EOF = fscanf(file, "%s", buffer);
         check_EOF = fscanf(file, "%lf\n", &buffer);
         if ( check_EOF == EOF )
         {
            fprintf( stderr, "Error: Unexpected EOF in read_coordinates\n" );
            exit(-1);
         }
         mesh->points[ip].coordinates[d] = buffer; //atof(buffer);
      }
   }

   // Close file
   fclose(file);
}

void read_faces ( char *filename, mesh_t *mesh )
{
   int iface, ielem, check_EOF;
   //char buffer[255];
   int buffer;
   FILE *file;

   // Open file
   file = fopen( filename, "r" );

   // First element must be nFaces
   //check_EOF = fscanf(file, "%s", buffer);
   check_EOF = fscanf(file, "%d\n", &buffer);
   if ( check_EOF == EOF )
   {
      fprintf( stderr, "Error: Unexpected EOF in read_faces\n" );
      exit(-1);
   }
   mesh->nFaces = buffer; //atoi(buffer);
   mesh->faces = malloc (sizeof(int) * mesh->nVertsPerFace * mesh->nFaces);
   if ( mesh->faces == NULL )
   {
      fprintf( stderr, "Error: Allocating memory for mesh->faces in read_faces function\n" );
      exit(-1);
   }

   // Rest of read elements will be faces points' indices
   for ( iface = 0; iface < mesh->nFaces; iface++ )
   {
      for ( ielem = 0; ielem < mesh->nVertsPerFace; ielem++ )
      {
         //check_EOF = fscanf(file, "%s", buffer);
         check_EOF = fscanf(file, "%d\n", &buffer);
         if ( check_EOF == EOF )
         {
            fprintf( stderr, "Error: Unexpected EOF in read_faces\n" );
            exit(-1);
         }
         mesh->faces[iface * mesh->nVertsPerFace + ielem] = buffer; // atoi(buffer);
      }
   }

   // Close file
   fclose(file);
}

void read_time ( char *filename, mesh_t *mesh )
{
   int it, check_EOF;
   //char buffer[255];
   double buffer;
   int buffer_size;
   FILE *file;

   // Open file
   file = fopen( filename, "r" );
   
   // First element must be nTimes
   //check_EOF = fscanf(file, "%s", buffer);
   check_EOF = fscanf(file, "%d\n", &buffer_size);
   if ( check_EOF == EOF )
   {
      fprintf( stderr, "Error: Unexpected EOF in read_time\n" );
      exit(-1);
   }
   mesh->nTimes = buffer_size;// atoi(buffer);
   mesh->times  = malloc (sizeof(double) * mesh->nTimes);
   if ( mesh->times == NULL )
   {
      fprintf( stderr, "Error: Allocating memory for mesh->times in read_time function\n" );
      exit(-1);
   }

   // Rest of read elements will be time data
   for ( it = 0; it < mesh->nTimes; it++ )
   {
      //check_EOF = fscanf(file, "%s", buffer);
      check_EOF = fscanf(file, "%lf\n", &buffer);
      if ( check_EOF == EOF )
      {
         fprintf( stderr, "Error: Unexpected EOF in read_time\n" );
         exit(-1);
      }
      mesh->times[it] = buffer; // atof(buffer);
   }

   // Close file
   fclose(file);
}

void read_velocity ( char *filename, mesh_t *mesh )
{
   int ip, it, d, check_EOF;
   //char buffer[255];
   double buffer;
   FILE *file;

   // Open file
   file = fopen( filename, "r" );

   // Set velocity vectors space
   for ( ip = 0; ip < mesh->nPoints; ip++ )
   {
      mesh->points[ip].velocity = malloc(sizeof(double) * mesh->nDim * mesh->nTimes);
      if ( mesh->points[ip].velocity == NULL )
      {
         fprintf( stderr, "Error: Allocating memory for mesh->points[i].velocity in read_velocity function\n" );
         exit(-1);
      }
   }

   // Read elements will be points' velocity data
   for ( it = 0; it < mesh->nTimes; it++ )
   {
      for ( ip = 0; ip < mesh->nPoints; ip++ )
      {
         for ( d = 0; d < mesh->nDim; d++ )
         {
            //check_EOF = fscanf(file, "%s", buffer);
            check_EOF = fscanf(file, "%lf\n", &buffer);
            if ( check_EOF == EOF )
            {
               fprintf( stderr, "Error: Unexpected EOF in read_velocity\n" );
               exit(-1);
            }            
            mesh->points[ip].velocity[it*mesh->nDim+d] = buffer; // atof(buffer);
         }
      }
   }

   // Close file
   fclose(file);
}

void assign_faces_to_verts ( mesh_t *mesh ) 
{
	int faces[mesh->nFaces];
	int ip, iface, ipf, count;
	for ( ip = 0; ip < mesh->nPoints; ip++ )
	{
		count = 0;
		for ( iface = 0; iface < mesh->nFaces; iface++ )
		{
			for ( ipf = 0; ipf < mesh->nVertsPerFace; ipf++ )
			{
				if ( mesh->faces[iface * mesh->nVertsPerFace + ipf] == ip )
				{
					faces[count] = iface;
					count++;
				}
			}
		}
		mesh->points[ip].nFaces = count;
		mesh->points[ip].faces = malloc ( sizeof(int) * count );
		for ( iface = 0; iface < count; iface++ )
		{
			mesh->points[ip].faces[iface] = faces[iface];
		}
	}
}
