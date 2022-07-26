#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "read_python_files.h"
#include "velocity.h"
#include "interpolation.h"
#include "compute_flowmap.h"
#include "rk4.h"
#include "omp.h"

int main(int argc, char *argv[])
{
   // Check usage
   if (argc != 13)
   {
	printf("USAGE: ./executable <nDim> <ntimes_eval> <t0_eval> <tdelta_eval> <coords_file> <faces_file> <times_file> <vel_file> <nsteps_rk4> <sched_policy> <print>\n");
	printf("\texecutable:  compute_flowmap (sequential)\n");
	printf("\tnDim:        dimensions of the space (2D/3D)\n");
	printf("\tntimes_eval: quantity of t values to choose (starting at t0 and adding tdelta each time) which we want to evaluate.\n");
	printf("\tt0_eval:      first t to evaluate.\n");
	printf("\ttdelta_eval:  t increment in each t to evaluate after t0.\n");
	printf("\tcoords_file:  file where mesh coordinates are stored.\n");
        printf("\tfaces_file:   file where mesh faces are stored.\n");
        printf("\ttimes_file:   file where original time data are stored.\n");
        printf("\tvel_file:     file where original velocity data is stored.\n");
        printf("\tnsteps_rk4:   number of iterations to perform in the RK4 call.\n");
	printf("\tsched_policy: SEQUENTIAL (1) / OMP_STATIC (2) / OMP_DYNAMIC (3) / OMP_GUIDED (4)\n");
        printf("\tchunk_size:   size of the chunk for the chosen scheduling policy\n");
	printf("\tprint to file? (0-NO, 1-YES)\n");
	return 1;
   }

   int d, ip, it, npoints, nteval, nsteps_rk4, itprev;

   FILE *fp_w;   

   double    t0, tdelta, time, t_start;

   struct timeval start;
   struct timeval end;


   int policy = atoi(argv[10]);
   int sched_chunk_size = atoi(argv[11]);

   mesh_t mesh;
   mesh.nDim = atoi(argv[1]);
   

   if ( mesh.nDim == 2 ) mesh.nVertsPerFace = 3; // 2D: faces are triangles
   else mesh.nVertsPerFace = 4; // 3D: faces (volumes) are tetrahedrons

   /* Read coordinates, faces, time values and velocity data from Python-generated files */
   printf("Reading mesh points coordinates... ");
   read_coordinates(argv[5], &mesh); 
   printf("DONE\n"); fflush(stdout);
   printf("Reading mesh faces faces...        ");
   read_faces(argv[6], &mesh); 
   printf("DONE\n"); fflush(stdout);
   printf("Reading mesh times...              ");
   read_time(argv[7], &mesh); 
   printf("DONE\n"); fflush(stdout);
   printf("Reading mesh velocity data...      ");
   read_velocity(argv[8], &mesh); 
   printf("DONE\n"); fflush(stdout);

   printf("Computing the flowmap with %d threads... ", omp_get_num_threads());

   /* Set variables where we want to solve the IVP (coords and time) */
   npoints = mesh.nPoints;
   nteval  = atoi(argv[2]);

   double t_eval[nteval];

   t0     = atof(argv[3]);
   tdelta = atof(argv[4]);
   for ( it = 0; it < nteval; it++ )
   {
      t_eval[it] = t0 + tdelta * it;
   }

   /* Set number of steps of each RK4 function call */
   nsteps_rk4 = atoi(argv[9]);

   /* Solve IVPs (using RK4) */
   double *result = malloc( sizeof(double) * npoints * nteval * mesh.nDim );
   if ( policy == 1 )
   {
      gettimeofday(&start, NULL);
      for ( ip = 0; ip < npoints; ip++ )
      {
         for ( it = 0; it < nteval; it++ )
         {
	    result[ip * nteval * mesh.nDim + it * mesh.nDim] = mesh.points[ip].coordinates[0];
	    result[ip * nteval * mesh.nDim + it * mesh.nDim + 1] = mesh.points[ip].coordinates[1];
            if (mesh.nDim == 3) result[ip * nteval * mesh.nDim + it * mesh.nDim + 2] = mesh.points[ip].coordinates[2];
	    itprev = 0;
            while ( itprev+1 < mesh.nTimes && mesh.times[itprev+1] < t_eval[it] )
            {
	       itprev++;
               runge_kutta_4 ( &result[ ip * nteval * mesh.nDim + it * mesh.nDim ], 
                         mesh.times[itprev-1],
                         mesh.times[itprev], 
                         &mesh, &result[ ip * nteval * mesh.nDim + it * mesh.nDim ], 
                         nsteps_rk4);
            }
            runge_kutta_4 ( &result[ ip * nteval * mesh.nDim + it * mesh.nDim ],
                         mesh.times[itprev],
                         t_eval[it],
                         &mesh, &result[ ip * nteval * mesh.nDim + it * mesh.nDim ],
                         nsteps_rk4);
         }
      }
      gettimeofday(&end, NULL);
   }
   else if ( policy == 2 )
   {
      gettimeofday(&start, NULL);
      #pragma omp parallel for default(none) shared(npoints, nteval, mesh, t_eval, nsteps_rk4, result, sched_chunk_size) private(ip, it, itprev) schedule(static, sched_chunk_size)
      for ( ip = 0; ip < npoints; ip++ )
      {
         for ( it = 0; it < nteval; it++ )
         {
            result[ip * nteval * mesh.nDim + it * mesh.nDim] = mesh.points[ip].coordinates[0];
            result[ip * nteval * mesh.nDim + it * mesh.nDim + 1] = mesh.points[ip].coordinates[1];
            if (mesh.nDim == 3) result[ip * nteval * mesh.nDim + it * mesh.nDim + 2] = mesh.points[ip].coordinates[2];
            itprev = 0;
            while ( itprev+1 < mesh.nTimes && mesh.times[itprev+1] < t_eval[it] )
            {
               itprev++;
               runge_kutta_4 ( &result[ ip * nteval * mesh.nDim + it * mesh.nDim ],
                         mesh.times[itprev-1],
                         mesh.times[itprev],
                         &mesh, &result[ ip * nteval * mesh.nDim + it * mesh.nDim ],
                         nsteps_rk4);
            }
            runge_kutta_4 ( &result[ ip * nteval * mesh.nDim + it * mesh.nDim ],
                         mesh.times[itprev],
                         t_eval[it],
                         &mesh, &result[ ip * nteval * mesh.nDim + it * mesh.nDim ],
                         nsteps_rk4);
         }
      }
      gettimeofday(&end, NULL);
   }
   else if ( policy == 3 )
   {
      gettimeofday(&start, NULL);
      #pragma omp parallel for default(none) shared(npoints, nteval, mesh, t_eval, nsteps_rk4, result, sched_chunk_size) private(ip, it, itprev) schedule(dynamic, sched_chunk_size)
      for ( ip = 0; ip < npoints; ip++ )
      {
         for ( it = 0; it < nteval; it++ )
         {
            result[ip * nteval * mesh.nDim + it * mesh.nDim] = mesh.points[ip].coordinates[0];
            result[ip * nteval * mesh.nDim + it * mesh.nDim + 1] = mesh.points[ip].coordinates[1];
            if (mesh.nDim == 3) result[ip * nteval * mesh.nDim + it * mesh.nDim + 2] = mesh.points[ip].coordinates[2];
            itprev = 0;
            while ( itprev+1 < mesh.nTimes && mesh.times[itprev+1] < t_eval[it] )
            {
               itprev++;
               runge_kutta_4 ( &result[ ip * nteval * mesh.nDim + it * mesh.nDim ],
                         mesh.times[itprev-1],
                         mesh.times[itprev],
                         &mesh, &result[ ip * nteval * mesh.nDim + it * mesh.nDim ],
                         nsteps_rk4);
            }
            runge_kutta_4 ( &result[ ip * nteval * mesh.nDim + it * mesh.nDim ],
                         mesh.times[itprev],
                         t_eval[it],
                         &mesh, &result[ ip * nteval * mesh.nDim + it * mesh.nDim ],
                         nsteps_rk4);
         }
      }
      gettimeofday(&end, NULL);
   }
   else if ( policy == 4 )
   {
      gettimeofday(&start, NULL);
      #pragma omp parallel for default(none) shared(npoints, nteval, mesh, t_eval, nsteps_rk4, result, sched_chunk_size) private(ip, it, itprev) schedule(guided, sched_chunk_size)
      for ( ip = 0; ip < npoints; ip++ )
      {
         for ( it = 0; it < nteval; it++ )
         {
            result[ip * nteval * mesh.nDim + it * mesh.nDim] = mesh.points[ip].coordinates[0];
            result[ip * nteval * mesh.nDim + it * mesh.nDim + 1] = mesh.points[ip].coordinates[1];
            if (mesh.nDim == 3) result[ip * nteval * mesh.nDim + it * mesh.nDim + 2] = mesh.points[ip].coordinates[2];
            itprev = 0;
            while ( itprev+1 < mesh.nTimes && mesh.times[itprev+1] < t_eval[it] )
            {
               itprev++;
               runge_kutta_4 ( &result[ ip * nteval * mesh.nDim + it * mesh.nDim ],
                         mesh.times[itprev-1],
                         mesh.times[itprev],
                         &mesh, &result[ ip * nteval * mesh.nDim + it * mesh.nDim ],
                         nsteps_rk4);
            }
            runge_kutta_4 ( &result[ ip * nteval * mesh.nDim + it * mesh.nDim ],
                         mesh.times[itprev],
                         t_eval[it],
                         &mesh, &result[ ip * nteval * mesh.nDim + it * mesh.nDim ],
                         nsteps_rk4);
         }
      }
      gettimeofday(&end, NULL);
   }
   else
   {
      fprintf( stderr, "Error: wrong scheduling policy\n" );
      exit(-1);
   }

   gettimeofday(&end, NULL);
   printf("DONE\n");

   /* Show execution time */
   time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1000000.0;
   printf("Exec time (npoints %d, ntimes %d) = %f\n", npoints, nteval, time);

   /* Print res to file */
   if ( atoi(argv[12]) )
   {
      double fsolfmx = 0, fsolfmy = 0, fsol_sum = 0;
      fp_w = fopen("output_for_ftle.csv", "w");
      for ( ip = 0; ip < npoints; ip++)
      {
         for ( it = 0; it < nteval; it++ )
         {
	    if ( mesh.nDim == 2 )
   	    {
	       fsolfmx += result[ip*nteval*mesh.nDim + it*mesh.nDim];
	       fsolfmy += result[ip*nteval*mesh.nDim + it*mesh.nDim+1];
	       fprintf(fp_w, "%1.14lf\n%1.14lf\n", result[ip*nteval*mesh.nDim + it*mesh.nDim], result[ip*nteval*mesh.nDim + it*mesh.nDim+1]);
	    }
	    else
            {
	       fprintf(fp_w, "%1.14lf\n%1.14lf\n%1.14lf\n", 
			result[ip*nteval*mesh.nDim + it*mesh.nDim], 
			result[ip*nteval*mesh.nDim + it*mesh.nDim+1],
			result[ip*nteval*mesh.nDim + it*mesh.nDim+2]);
            }
         }
      }
      fclose(fp_w);
      printf("Fsol fmx: %f, fsol fmy: %f, fsol sum: %f\n", fsolfmx, fsolfmy, fsolfmx+fsolfmy);
   }
   free (result);

   return 0;
}
