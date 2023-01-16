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

#define pi 3.142857

void dgyre_velocity ( double x, double y, double t, double *vel )
{
   double A = 0.1;
   double omega = 2 * pi/10;
   double epsilon = 0.25;
   double a_t = epsilon*sin(omega*t);
   double b_t = 1 - 2*epsilon*sin(omega*t);
   double f = a_t*x*x+b_t*x;
   double dfdx = 2*a_t*x+b_t;
   vel[0] = -A*pi*sin(pi*f)*cos(pi*y);
   vel[1] = pi * A*cos(pi*f)*sin(pi*y)*dfdx;
}

void abc_velocity ( double x, double y, double z, double t, double *vel )
{
   double A = sqrt(3);
   double B = sqrt(2);
   double C = 1;
   double omega = 2*pi/10;
   double epsilon = 0.1;
   double A_t = A + epsilon * cos(omega*t);
   vel[0] = A_t * sin(z) + C * cos(y);
   vel[1] = B * sin(x) + A_t * cos(z);
   vel[2] = C * sin(y) + B * cos(x);
}

int main(int argc, char *argv[])
{
   struct timeval globalstart;
   struct timeval globalend;  
   struct rusage r_usage;
   int ret;

   gettimeofday(&globalstart, NULL);
   // Check usage
   if (argc != 13)
   {
	printf("USAGE: ./executable <nDim> <t_eval> <faces_file> <times_file> <nsteps_rk4> <sched_policy> <print> <nx> <ny> <nz> <nt> <tlim>\n");
	printf("\texecutable:   compute_flowmap\n");
	printf("\tnDim:         dimensions of the space (2D/3D)\n");
	printf("\tt_eval:       t to evaluate.\n");
        printf("\tfaces_file:   file where mesh faces are stored.\n");
        printf("\ttimes_file:   file where original time data are stored.\n");
        printf("\tnsteps_rk4:   number of iterations to perform in the RK4 call.\n");
	printf("\tsched_policy: SEQUENTIAL (1) / OMP_STATIC (2) / OMP_DYNAMIC (3) / OMP_GUIDED (4)\n");
	printf("\tprint to file? (0-NO, 1-YES)\n");
        printf("\tnx:           number of points in X axis\n");
        printf("\tny:           number of points in Y axis\n");
        printf("\tnz:           number of points in Z axis\n");
        printf("\tnt:           number of known time instants\n");
        printf("\ttlim:         max known time instant\n");
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
   double tlim = atof(argv[12]);

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

   /* Number of points in each axis */
   int nx = atoi(argv[8]);
   int ny = atoi(argv[9]);
   int nz = atoi(argv[10]);
   int nt = atoi(argv[11]);

   /* Time instant to evaluate (from input arguments) */
   t_eval = atof(argv[2]);

   /* Read coordinates, faces, times and velocities information                    */
   /* from the Python-generated files and generate the corresponding coords vector */

   // 1. Set nPoints
   nPoints = nx * ny;
   if ( nDim == 3 ) nPoints = nPoints * nz;

   // 2. Get nFaces and faces ids
   file = fopen( argv[3], "r" );
   check_EOF = fscanf(file, "%s", buffer);
   if ( check_EOF == EOF )
   {
      fprintf( stderr, "Error: Unexpected EOF in read_faces\n" );
      exit(-1);
   }
   nFaces = atoi(buffer);   
   fclose(file);

   // 3. Get nTimes
   nTimes = nt;
   double dx = (double) xlim/ (double) nx; double xval;
   double dy = (double) ylim/ (double) ny; double yval;
   double dz;
   double dt;

   // 4. Generate coordinates
   printf("Reading mesh points coordinates...                     ");
   coords_x = (double *) malloc ( sizeof(double) * nPoints );
   coords_y = (double *) malloc ( sizeof(double) * nPoints );
   if ( nDim == 3) coords_z = (double *) malloc ( sizeof(double) * nPoints );

   ip = 0;
   int ix, iy, iz;
   if ( nDim == 2 )
   {
      for (ix = 0; ix < nx; ix++ )
      {
         xval = ( ix == (nx - 1) ) ? xlim : dx * ix;
         #pragma omp parallel for default(none) shared(coords_x, coords_y, ny, ylim, dy) private(ip, iy, yval) firstprivate(ix, xval)
         for (iy = 0; iy < ny; iy++ )
         {
            ip = ix * ny + iy;
            yval = ( iy == (ny - 1) ) ? ylim : dy * iy;
            coords_x[ip] = xval;
            coords_y[ip] = yval;
         }
      }
   }
   else // nDim==3
   {
      dz = (double) zlim/ (double) nz;
      for (ix = 0; ix < nx; ix++ )
      {
         xval = ( ix == (nx - 1) ) ? xlim : dx * ix;
         for (iy = 0; iy < ny; iy++ )
         {
            yval = ( iy == (ny - 1) ) ? ylim : dy * iy;
            #pragma omp parallel for default(none) shared(coords_x, coords_y, coords_z, ny, ylim, dy, nz, zlim, dz) private(ip, iz) firstprivate(ix, iy, xval, yval)
            for (iz = 0; iz < nz; iz++ )
            {
               ip = ix * ny * nz + iy * nz + iz;
               coords_x[ip] = xval;
               coords_y[ip] = yval;
               coords_z[ip] = ( iz == (nz - 1) ) ? zlim : dz * iz;
            }
         }
      }
   }
   printf("DONE\n");

   // 4.1 Populate kdtree with coords
   printf("Populating kdtree with mesh coordinates...             ");
   if (nDim == 2)
	for ( ip = 0; ip < nPoints; ip++ )
      		assert(kd_insert2(kd, coords_x[ip], coords_y[ip], ip) == 0);
   else
        for ( ip = 0; ip < nPoints; ip++ )
                assert(kd_insert3(kd, coords_x[ip], coords_y[ip], coords_z[ip], ip) == 0);
 
   printf("DONE\n");

   // 5. Read mesh faces values from the given file
   printf("Reading mesh faces vertices...                         "); 
   faces = (int *) malloc ( sizeof(int) * nFaces * nVertsPerFace );
   read_faces(argv[3], nDim, nVertsPerFace, nFaces, faces); 
   printf("DONE\n");

   // 5.1 Create auxiliar faces per point vectors information
   printf("Creating auxiliar faces info per point...              ");
   nFacesPerPoint = malloc ( sizeof(int) * nPoints );
   create_nFacesPerPoint_vector ( nDim, nPoints, nFaces, nVertsPerFace, faces, nFacesPerPoint );
   facesPerPoint  = malloc ( sizeof(int) * nFacesPerPoint[nPoints-1] );
   create_facesPerPoint_vector ( nDim, nPoints, nFaces, nVertsPerFace, faces, nFacesPerPoint, facesPerPoint );
   printf("DONE\n");

   // 6. Read time values from the given file
   printf("Reading known time values...                           ");
   times = (double *) malloc ( sizeof(double) * nTimes );
   read_times(argv[4], nTimes, times);
   printf("DONE\n");

   // 7. Read velocities from the given file
   printf("Reading known velocity values...                       ");fflush(stdout);
   velocities = (double *) malloc ( sizeof(double) * nTimes * nPoints * nDim );
   int iv = 0, itim;
   if ( nDim == 2 )
   {
      for ( itim = 0; itim < nTimes; itim++ )
      {
         for ( ip = 0; ip < nPoints; ip++ )
         {
            dgyre_velocity(coords_x[ip], coords_y[ip], times[itim], &velocities[iv]);
            iv = iv + 2;
         }
      }
   }
   else // nDim==3
   {
     for ( itim = 0; itim < nTimes; itim++ )
      {
         for ( ip = 0; ip < nPoints; ip++ )
         {
            abc_velocity(coords_x[ip], coords_y[ip], coords_z[ip], times[itim], &velocities[iv]);
            iv = iv + 3;
         }
      }
   }
   printf("DONE\n");

   /* Set number of steps of each RK4 function call */
   nsteps_rk4 = atoi(argv[5]);
 
   /* Set scheduling policy */
   policy = atoi(argv[6]);
   if ( policy < 1 || policy > 4 )
   {
      printf("Wrong sched policy. Should be a 1-4 value.\nSEQUENTIAL (1) / OMP_STATIC (2) / OMP_DYNAMIC (3) / OMP_GUIDED (4)\n");
      return 1;
   }

   #pragma omp parallel
   {
	if ( omp_get_thread_num() == 0 )
		printf("\nComputing the flowmap with %d threads ", omp_get_num_threads()); fflush(stdout);
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
   printf("(memory usage = %ld KB)... ",r_usage.ru_maxrss); fflush(stdout);

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
   printf("\nFlowmap computation elapsed time = %f seconds\n", nPoints, time);

   /* Print res to file */
   if ( atoi(argv[7]) )
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
   free (result);
   kd_free(kd);

   gettimeofday(&globalend, NULL);
   time = (globalend.tv_sec - globalstart.tv_sec) + (globalend.tv_usec - globalstart.tv_usec)/1000000.0;
   printf("TOTAL elapsed time (including reading) = %f seconds\n", nPoints, time);
   return 0;
}
