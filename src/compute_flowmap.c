#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "read_python_files.h"
#include "interpolation.h"
#include "rk4.h"
#include "omp.h"
#include "kdtree.h"
#include <assert.h>
#include <sys/resource.h>

int main(int argc, char *argv[])
{
   struct timeval globalstart;
   struct timeval globalend;  
   struct rusage r_usage;
   int ret;

   gettimeofday(&globalstart, NULL);
   // Check usage
   if (argc != 10)
   {
	printf("USAGE: ./executable <nDim> <t_eval> <coords_file> <faces_file> <times_file> <vel_file> <nsteps_rk4> <sched_policy> <print>\n");
	printf("\texecutable:   compute_flowmap\n");
	printf("\tnDim:         dimensions of the space (2D/3D)\n");
	printf("\tt_eval:       t to evaluate.\n");
	printf("\tcoords_file:  file where mesh coordinates are stored.\n");
        printf("\tfaces_file:   file where mesh faces are stored.\n");
        printf("\ttimes_file:   file where original time data are stored.\n");
        printf("\tvel_file:     file where original velocity data is stored.\n");
        printf("\tnsteps_rk4:   number of iterations to perform in the RK4 call.\n");
	printf("\tsched_policy: SEQUENTIAL (1) / OMP_STATIC (2) / OMP_DYNAMIC (3) / OMP_GUIDED (4)\n");
	printf("\tprint to file? (0-NO, 1-YES)\n");
        return 1;
   }

   int d, ip, it, itprev;

   FILE  *fp_w;   

   int    nDim;
   int    nPoints;
   int    nTimes;
   int    nFaces;
   int    nVertsPerFace;
   int    policy;
   int    sched_chunk_size;
   int    nsteps_rk4;

   double t_eval;
   double time;

   int    *faces;
   int    *nFacesPerPoint;
   int    *facesPerPoint;

   double *coords_x;
   double *coords_y;
   double *coords_z = NULL;

   double *times;

   double *velocities;

   struct timeval start;
   struct timeval end;

   int  check_EOF;
   char buffer[255];
   FILE *file;

   struct kdtree *kd;
   int xlim, ylim, zlim;

   /* Mesh dimension obtained (from input arguments) */
   /* Number of vertices per face according to dim   */
   /* Create kdtree struct                           */
   nDim = atoi(argv[1]);
   if ( nDim == 2 )
   {
      kd = kd_create(2);
      nVertsPerFace = 3; // 2D: faces are triangles
      xlim = 2; ylim = 1; // 2D: double gire flow
   }
   else
   {
      if ( nDim == 3)
      {
         kd = kd_create(3);
         nVertsPerFace = 4; // 3D: faces (volumes) are tetrahedrons
         xlim = 1; ylim = 1; zlim = 1; // 3D: ABC flow
      }
      else
      {
         printf("Wrong dimension provided (2 or 3 supported)\n"); 
         return 1;
      }
   }

   /* Time instant to evaluate (from input arguments) */
   t_eval = atof(argv[2]);

   /* Read coordinates, faces, times and velocities information                    */
   /* from the Python-generated files and generate the corresponding coords vector */

   // 1. Get nPoints
   file = fopen( argv[3], "r" );
   check_EOF = fscanf(file, "%s", buffer);
   if ( check_EOF == EOF )
   {
      fprintf( stderr, "Error: Unexpected EOF in read_coordinates\n" );
      exit(-1);
   }
   nPoints = atoi(buffer);
   fclose(file);

   // 2. Get nFaces
   file = fopen( argv[4], "r" );
   check_EOF = fscanf(file, "%s", buffer);
   if ( check_EOF == EOF )
   {
      fprintf( stderr, "Error: Unexpected EOF in read_faces\n" );
      exit(-1);
   }
   nFaces = atoi(buffer);   
   fclose(file);

   // 3. Get nTimes
   file = fopen( argv[5], "r" );
   check_EOF = fscanf(file, "%s", buffer);
   if ( check_EOF == EOF )
   {
      fprintf( stderr, "Error: Unexpected EOF in read_times\n" );
      exit(-1);
   }
   nTimes = atoi(buffer);
   fclose(file);

   // 4. Read mesh coords values from the given file
   printf("Reading mesh points coordinates...                                  ");
   coords_x = (double *) malloc ( sizeof(double) * nPoints );
   coords_y = (double *) malloc ( sizeof(double) * nPoints );
   if ( nDim == 3) coords_z = (double *) malloc ( sizeof(double) * nPoints );
   read_coordinates(argv[3], nDim, nPoints, coords_x, coords_y, coords_z); 
   printf("DONE\n");

   // 4.1 Populate kdtree with coords
   printf("Populating kdtree with mesh coordinates...                          ");
   if (nDim == 2)
	for ( ip = 0; ip < nPoints; ip++ )
      		assert(kd_insert2(kd, coords_x[ip], coords_y[ip], ip) == 0);
   else
        for ( ip = 0; ip < nPoints; ip++ )
                assert(kd_insert3(kd, coords_x[ip], coords_y[ip], coords_z[ip], ip) == 0);
 
   printf("DONE\n");

   // 5. Read mesh faces values from the given file
   printf("Reading mesh faces vertices...                                      "); 
   faces = (int *) malloc ( sizeof(int) * nFaces * nVertsPerFace );
   read_faces(argv[4], nDim, nVertsPerFace, nFaces, faces); 
   printf("DONE\n");

   // 5.1 Create auxiliar faces per point vectors information
   
   printf("Creating auxiliar faces info per point...                           ");
   nFacesPerPoint = malloc ( sizeof(int) * nPoints );
   create_nFacesPerPoint_vector ( nDim, nPoints, nFaces, nVertsPerFace, faces, nFacesPerPoint );
   facesPerPoint  = malloc ( sizeof(int) * nFacesPerPoint[nPoints-1] );
   create_facesPerPoint_vector ( nDim, nPoints, nFaces, nVertsPerFace, faces, nFacesPerPoint, facesPerPoint );
   printf("DONE\n");

   // 6. Read time values from the given file
   printf("Reading known time values...                                        ");
   times = (double *) malloc ( sizeof(double) * nTimes );
   read_times(argv[5], nTimes, times);
   printf("DONE\n");

   // 7. Read velocities from the given file
   printf("Reading known velocity values...                                    ");fflush(stdout);
   velocities = (double *) malloc ( sizeof(double) * nTimes * nPoints * nDim );
   read_velocities(argv[6], nPoints, nDim, nTimes, velocities);
   printf("DONE\n");

   /* Set number of steps of each RK4 function call */
   nsteps_rk4 = atoi(argv[7]);
 
   /* Set scheduling policy and its chunk size */
   policy = atoi(argv[8]);
   if ( policy < 1 || policy > 4 )
   {
      printf("Wrong sched policy. Should be a 1-4 value.\nSEQUENTIAL (1) / OMP_STATIC (2) / OMP_DYNAMIC (3) / OMP_GUIDED (4)\n");
      return 1;
   }

   #pragma omp parallel
   {
	if (omp_get_thread_num() == 1)
	switch ( policy )
	{
		case 1:
			printf("\nComputing the flowmap using the sequential version ");
			break;

		case 2:
			printf("\nComputing the flowmap using the OpenMP-Static version with %d threads ", omp_get_num_threads());
			break;
	
		case 3:
			printf("\nComputing the flowmap using the OpenMP-Dynamic version with %d threads ", omp_get_num_threads());
                        break;

		case 4:
			printf("\nComputing the flowmap using the OpenMP-Guyided version with %d threads ", omp_get_num_threads());
                        break;
	}		
	fflush(stdout);
   }

   /* Solve IVPs (using RK4) */
   double *result = malloc( sizeof(double) * nPoints * nDim );
   #pragma omp parallel for default(none) shared(result, coords_x, coords_y, nDim, nPoints) private(ip)
   for ( ip = 0; ip < nPoints; ip++ )
   {
   	result[ip * nDim]     = coords_x[ip];
        result[ip * nDim + 1] = coords_y[ip];
   }
   if ( nDim == 3 )
   {
	#pragma omp parallel for default(none) shared(result, coords_z, nDim, nPoints) private(ip)
   	for ( ip = 0; ip < nPoints; ip++ )
   	{
        	result[ip * nDim + 2] = coords_z[ip];
   	}
   }

   // Get memory usage information
   getrusage(RUSAGE_SELF,&r_usage);
   printf("(memory usage = %ldKB)...  ",r_usage.ru_maxrss); fflush(stdout);

   if ( policy == 1 )
   {
      gettimeofday(&start, NULL);
      for ( ip = 0; ip < nPoints; ip++ )
      {
	    itprev = 0;
            while ( ( itprev+1 < nTimes ) && ( times[itprev+1] < t_eval ) )
            {
	       itprev++;
               runge_kutta_4 ( result[ ip * nDim ], result[ ip * nDim + 1 ], result[ ip * nDim + 2], 
                         times[itprev-1],
                         times[itprev], 
                         &result[ ip * nDim ], 
                         nsteps_rk4,
			 nDim, nPoints, nTimes, times, 
			 nVertsPerFace, nFaces, faces, 
			 coords_x, coords_y, coords_z, velocities, kd, nFacesPerPoint, facesPerPoint);
            }
            runge_kutta_4 ( result[ ip * nDim ], result[ ip * nDim + 1 ], result[ ip * nDim + 2],
                         times[itprev],
                         t_eval,
                         &result[ ip * nDim ],
                         nsteps_rk4,
			 nDim, nPoints, nTimes, times,
			 nVertsPerFace, nFaces, faces, 
			 coords_x, coords_y, coords_z, velocities, kd, nFacesPerPoint, facesPerPoint);
      }
      gettimeofday(&end, NULL);
   }
   else if ( policy == 2 )
   {
      gettimeofday(&start, NULL);
      #pragma omp parallel for default(none) shared(nFacesPerPoint, facesPerPoint, kd, nPoints, nDim, coords_x, coords_y, coords_z, velocities, times, nTimes, t_eval, nsteps_rk4, result, nVertsPerFace, nFaces, faces) private(ip, it, itprev, r_usage) schedule(static)
      for ( ip = 0; ip < nPoints; ip++ )
      {
            itprev = 0;
            while ( ( itprev+1 < nTimes ) && ( times[itprev+1] < t_eval ) )
            {
               itprev++;
               runge_kutta_4 ( result[ ip * nDim ], result[ ip * nDim + 1 ], result[ ip * nDim + 2],
                         times[itprev-1],
                         times[itprev],      
                         &result[ ip * nDim ],                                      
                         nsteps_rk4,
			 nDim, nPoints, nTimes, times,
			 nVertsPerFace, nFaces, faces, 
			 coords_x, coords_y, coords_z, velocities, kd, nFacesPerPoint, facesPerPoint);
            }
            runge_kutta_4 ( result[ ip * nDim ], result[ ip * nDim + 1 ], result[ ip * nDim + 2],
                         times[itprev],
                         t_eval,
                         &result[ ip * nDim ],
                         nsteps_rk4,
			 nDim, nPoints, nTimes, times,
			 nVertsPerFace, nFaces, faces, 
			 coords_x, coords_y, coords_z, velocities, kd, nFacesPerPoint, facesPerPoint);
      }
      gettimeofday(&end, NULL);
   }
   else if ( policy == 3 )
   {
      gettimeofday(&start, NULL);
      #pragma omp parallel for default(none) shared(nFacesPerPoint, facesPerPoint, kd, nPoints, nDim, coords_x, coords_y, coords_z, velocities, times, nTimes, t_eval, nsteps_rk4, result, nVertsPerFace, nFaces, faces) private(ip, it, itprev) schedule(dynamic)
      for ( ip = 0; ip < nPoints; ip++ )
      {
            itprev = 0;
            while ( ( itprev+1 < nTimes ) && ( times[itprev+1] < t_eval ) )
            {
               itprev++;
               runge_kutta_4 ( result[ ip * nDim ], result[ ip * nDim + 1 ], result[ ip * nDim + 2],
                         times[itprev-1],
                         times[itprev],      
                         &result[ ip * nDim ],                                      
                         nsteps_rk4,
			 nDim, nPoints, nTimes, times,
			 nVertsPerFace, nFaces, faces, 
			 coords_x, coords_y, coords_z, velocities, kd, nFacesPerPoint, facesPerPoint);
            }
            runge_kutta_4 ( result[ ip * nDim ], result[ ip * nDim + 1 ], result[ ip * nDim + 2],
                         times[itprev],
                         t_eval,
                         &result[ ip * nDim ],
                         nsteps_rk4,
			 nDim, nPoints, nTimes, times,
			 nVertsPerFace, nFaces, faces, 
			 coords_x, coords_y, coords_z, velocities, kd, nFacesPerPoint, facesPerPoint);
      }
      gettimeofday(&end, NULL);
   }
   else if ( policy == 4 )
   {
      gettimeofday(&start, NULL);
      #pragma omp parallel for default(none) shared(nFacesPerPoint, facesPerPoint, kd, nPoints, nDim, coords_x, coords_y, coords_z, velocities, times, nTimes, t_eval, nsteps_rk4, result, nVertsPerFace, nFaces, faces) private(ip, it, itprev) schedule(guided)
      for ( ip = 0; ip < nPoints; ip++ )
      {
            itprev = 0;
            while ( ( itprev+1 < nTimes ) && ( times[itprev+1] < t_eval ) )
            {
               itprev++;
               runge_kutta_4 ( result[ ip * nDim ], result[ ip * nDim + 1 ], result[ ip * nDim + 2],
                         times[itprev-1],
                         times[itprev],      
                         &result[ ip * nDim ],                                      
                         nsteps_rk4,
			 nDim, nPoints, nTimes, times,
			 nVertsPerFace, nFaces, faces, 
			 coords_x, coords_y, coords_z, velocities, kd, nFacesPerPoint, facesPerPoint);
            }
            runge_kutta_4 ( result[ ip * nDim ], result[ ip * nDim + 1 ], result[ ip * nDim + 2],
                         times[itprev],
                         t_eval,
                         &result[ ip * nDim ],
                         nsteps_rk4,
			 nDim, nPoints, nTimes, times,
			 nVertsPerFace, nFaces, faces, 
			 coords_x, coords_y, coords_z, velocities, kd, nFacesPerPoint, facesPerPoint);
      }
      gettimeofday(&end, NULL);
   }
   else
   {
      fprintf( stderr, "Error: wrong scheduling policy\n" );
      exit(-1);
   }

   printf("DONE\n");

   /* Show execution time */
   time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1000000.0;
   printf("\nFlowmap computation elapsed time = %f seconds\n", time);

   /* Print res to file */
   if ( atoi(argv[9]) )
   {
      fp_w = fopen("output_for_ftle.csv", "w");
      for ( ip = 0; ip < nPoints; ip++)
      {
	    if ( nDim == 2 )
   	    {
	       fprintf(fp_w, "%1.14lf\n%1.14lf\n", 
                        result[ ip * nDim ], 
                        result[ ip * nDim + 1]);
	    }
	    else
            {
	       fprintf(fp_w, "%1.14lf\n%1.14lf\n%1.14lf\n", 
			result[ ip * nDim ], 
			result[ ip * nDim + 1],
			result[ ip * nDim + 2]);
            }
      }
      fclose(fp_w);
   }

   gettimeofday(&globalend, NULL);
   time = (globalend.tv_sec - globalstart.tv_sec) + (globalend.tv_usec - globalstart.tv_usec)/1000000.0;
   printf("TOTAL elapsed time (including reading) = %f seconds\n", time);

   // Free
   kd_free(kd);
   free (coords_x);
   free (coords_y);
   if ( nDim == 3 ) free (coords_z);
   free (velocities);
   free (times);
   free (result);
   free (faces);
   free (nFacesPerPoint);
   free (facesPerPoint);

   return 0;
}
