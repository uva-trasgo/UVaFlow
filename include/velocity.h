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

/* Compute velocity [x+t, y+t, z+t] */
/* Supported dimensions: 2D, 3D */
/* OUT: velocity (double *) */
void compute_velocity ( double t, double *coords, mesh_t *mesh, double *velocity );

/* Compute velocity vector as v_i = [x_i+t, y_i+t, z_i+t] */
/* Supported dimensions: 2D, 3D */
/* OUT: velocity (double *) */
void compute_velocity_vector ( mesh_t *mesh );
