#include <stdio.h>
#include <stdlib.h>
#include "read_python_files.h"

void read_coordinates ( char *filename, int nDim, int nPoints, double *coords_x, double *coords_y, double *coords_z )
{
	int ip, d, check_EOF;
	char buffer[255];
	FILE *file;

	// Open file
	file = fopen( filename, "r" );

	// First element must be nPoints
	check_EOF = fscanf(file, "%s", buffer);
	if ( check_EOF == EOF )
	{
		fprintf( stderr, "Error: Unexpected EOF in read_coordinates\n" );
		exit(-1);
	}

	// Rest of read elements will be points' coordinates
	for ( ip = 0; ip < nPoints; ip++ )
	{
		check_EOF = fscanf(file, "%s", buffer);
                if ( check_EOF == EOF )
                {
                        fprintf( stderr, "Error: Unexpected EOF in read_coordinates\n" );
                        exit(-1);
                }
                coords_x[ip] = atof(buffer);

		check_EOF = fscanf(file, "%s", buffer);
                if ( check_EOF == EOF )
                {
                        fprintf( stderr, "Error: Unexpected EOF in read_coordinates\n" );
                        exit(-1);
                }
                coords_y[ip] = atof(buffer);

		if ( nDim == 3 )
		{
			check_EOF = fscanf(file, "%s", buffer);
        	        if ( check_EOF == EOF )
                	{
                        	fprintf( stderr, "Error: Unexpected EOF in read_coordinates\n" );
	                        exit(-1);
        	        }
                	coords_z[ip] = atof(buffer);
		}
	}

	// Close file
	fclose(file);
}

void read_faces ( char *filename, int nDim, int nVertsPerFace, int nFaces, int *faces )
{
   int iface, ielem, check_EOF;
   char buffer[255];
   FILE *file;

   // Open file
   file = fopen( filename, "r" );

   // First element must be nFaces
   check_EOF = fscanf(file, "%s", buffer);
   if ( check_EOF == EOF )
   {
      fprintf( stderr, "Error: Unexpected EOF in read_faces\n" );
      exit(-1);
   }

   // Rest of read elements will be faces points' indices
   for ( iface = 0; iface < nFaces; iface++ )
   {
      for ( ielem = 0; ielem < nVertsPerFace; ielem++ )
      {
         check_EOF = fscanf(file, "%s", buffer);
         if ( check_EOF == EOF )
         {
            fprintf( stderr, "Error: Unexpected EOF in read_faces\n" );
            exit(-1);
         }
         faces[iface * nVertsPerFace + ielem] = atoi(buffer);
      }
   }

   // Close file
   fclose(file);
}

void read_times ( char *filename, int nTimes, double *times) 
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

   // Rest of read elements will be time data
   for ( it = 0; it < nTimes; it++ )
   {
      //check_EOF = fscanf(file, "%s", buffer);
      check_EOF = fscanf(file, "%lf\n", &buffer);
      if ( check_EOF == EOF )
      {
         fprintf( stderr, "Error: Unexpected EOF in read_time\n" );
         exit(-1);
      }
      times[it] = buffer; // atof(buffer);
   }

   // Close file
   fclose(file);
}

void read_velocities ( char *filename, int nPoints, int nDim, int nTimes, double *velocity )
{
   int ip, it, d, check_EOF;
   //char buffer[255];
   double buffer;
   FILE *file;

   // Open file
   file = fopen( filename, "r" );

   // Read elements will be points' velocity data
   for ( it = 0; it < nTimes; it++ )
   {
      for ( ip = 0; ip < nPoints; ip++ )
      {
         for ( d = 0; d < nDim; d++ )
         {
            //check_EOF = fscanf(file, "%s", buffer);
            check_EOF = fscanf(file, "%lf\n", &buffer);
            if ( check_EOF == EOF )
            {
               fprintf( stderr, "Error: Unexpected EOF in read_velocity\n" );
               exit(-1);
            }            
            velocity[it*nPoints*nDim +ip*nDim + d] = buffer; // atof(buffer);
         }
      }
   }

   // Close file
   fclose(file);
}

/*
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
*/

void create_nFacesPerPoint_vector ( int nDim, int nPoints, int nFaces, int nVertsPerFace, int *faces, int *nFacesPerPoint )
{
	int ip, iface, ipf;
	for ( ip = 0; ip < nPoints; ip++ )
        {
		nFacesPerPoint[ip] = 0;
	}
	for ( iface = 0; iface < nFaces; iface++ )
	{
		for ( ipf = 0; ipf < nVertsPerFace; ipf++ )
		{
			ip = faces[iface * nVertsPerFace + ipf];
			nFacesPerPoint[ip] = nFacesPerPoint[ip] + 1;
		}
	}
	for ( ip = 1; ip < nPoints; ip++ )
        {
                nFacesPerPoint[ip] = nFacesPerPoint[ip] + nFacesPerPoint[ip-1];
        }	
}

void create_facesPerPoint_vector ( int nDim, int nPoints, int nFaces, int nVertsPerFace, int *faces, int *nFacesPerPoint, int *facesPerPoint )
{
	int ip, count, iface, ipf, nFacesP, iFacesP;

        for ( ip = 0; ip < nPoints; ip++ )
        {       
                count   = 0;
		iFacesP = ( ip == 0 ) ? 0 : nFacesPerPoint[ip-1];
		nFacesP = ( ip == 0 ) ? nFacesPerPoint[ip] : nFacesPerPoint[ip] - nFacesPerPoint[ip-1];
                for ( iface = 0; ( iface < nFaces ) && ( count < nFacesP ); iface++ )
                {     
                      for ( ipf = 0; ipf < nVertsPerFace; ipf++ )
                      {       
                              if ( faces[iface * nVertsPerFace + ipf] == ip )
                              {
					facesPerPoint[iFacesP + count] = iface;
					count++;
                              }
                      }
                }
        }
}

