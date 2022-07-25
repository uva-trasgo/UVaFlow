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

/* Function to set all mesh points' coordinates read from file "filename" */
/* Supported dimensions: 2D, 3D */
/* NOTE: 1st line of given file must contain mesh->nPoints */
/* OUT: mesh->points[0...nPoints].coordinates and mesh->points[0...nPoints].index set */
void read_coordinates ( char *filename, mesh_t *mesh );

/* Function to set all mesh faces read from file "filename" */
/* Supported dimensions: 2D, 3D */
/* NOTE: 1st line of given file must contain mesh->nFaces */
/* OUT: mesh->faces and mesh->nFaces set */
void read_faces ( char *filename, mesh_t *mesh );

/* Function to set all mesh (velocity) times read from file "filename" */
/* Supported dimensions: 2D, 3D */
/* NOTE: 1st line of given file must contain mesh->nTimes */
/* OUT: mesh->times and mesh->nTimes set */
void read_time ( char *filename, mesh_t *mesh );

/* Function to set all mesh points' velocity read from file "filename" */
/* Supported dimensions: 2D, 3D */
/* OUT: mesh->points[0...nPoints].velocity set */
void read_velocity ( char *filename, mesh_t *mesh );

/* Function to set the lists of faces in which each point is a vertex (and nFaces) */
/* Supported dimensions: 2D, 3D */
/* OUT: mesh->points[0...nPoints].faces and mesh->points[0...nPoints].nFaces set */
void assign_faces_to_verts ( mesh_t *mesh );
