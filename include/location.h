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

/* Check if a given point is inside a triangle described by V1, V2, V3 vertices */
/* Supported dimensions: 2D */
/* OUT: int 1 (yes) / int 0 (no) */
int p_inside_triangle ( mesh_t *mesh, double *V1, double *V2, double *V3, double *P );

/* Check if a given point is inside a mesh face */
/* Supported dimensions: 2D */
/* OUT: point index in mesh->points array (-1 if the given point is not the mesh) */
int find_point_simplex_in_mesh ( double *Pcoords, mesh_t *mesh );

/* Check if a given point is or not a mesh point */
/* Supported dimensions: 2D, 3D */
/* OUT: point index in mesh->points array (-1 if the given point is not the mesh) */
int find_point_in_mesh ( double *Pcoords, mesh_t *mesh );

/* Check if a given point is inside a tetrahedron described by V1, V2, V3, V4 vertices */
/* share_side is a support function */
/* Supported dimensions: 3D */
/* OUT: int 1 (yes) / int 0 (no) */
int p_inside_tetrahedron ( mesh_t *mesh, double *V1, double *V2, double *V3, double *V4, double *P );
int share_side ( mesh_t *mesh, double *V1, double *V2, double *V3, double *V4, double *P );
