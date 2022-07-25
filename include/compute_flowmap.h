#ifndef STRUCTS_H
#define STRUCTS_H
typedef struct Point {
   int       index; // Point index in mesh->points
   double   *coordinates;
   double   *velocity;
   int       nfaces;
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
