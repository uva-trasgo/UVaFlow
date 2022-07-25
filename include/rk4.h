#ifndef STRUCTS_H
#define STRUCTS_H
typedef struct Point {
   int       index; // Point index in mesh->points
   double   *coordinates;
   double   *velocity;
   int       nFaces;
   int      *faces;
} point_t;

typedef struct Mesh {
   int       nDim;
   int       nPoints;
   int       nFaces;
   int       nVertsPerFace;
   int       nTimes;
   point_t  *points;
   int      *faces;
   double   *times;
} mesh_t;
#endif

/* Function to apply Runge-Kutta (4th order) method */
/* Supported dimensions: 2D, 3D */
/* OUT: result containing the estimation for 2D/3D */
void runge_kutta_4 ( double *Pcoords, double t0, double tend, mesh_t *mesh, double *result, int nsteps); //, int it );
//void runge_kutta_4 ( double *Pcoords, double t0, double *y0, double tend, mesh_t *mesh, double *result, int nsteps );
void runge_kutta_45 ( double *Pcoords, double t0, double tend, mesh_t *mesh, double *result, int nsteps); //, int it );
//void runge_kutta_45 ( double *Pcoords, double t0, double *y0, double tend, mesh_t *mesh, double *result, int nsteps );
