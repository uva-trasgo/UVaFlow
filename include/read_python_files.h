/* Function to set all mesh points' coordinates read from file "filename" */
/* Supported dimensions: 2D, 3D */
/* NOTE: 1st line of given file must contain mesh->nPoints */
/* OUT: mesh->points[0...nPoints].coordinates and mesh->points[0...nPoints].index set */
void read_coordinates ( char *filename, int nDim, int nPoints, double *coords );

/* Function to set all mesh faces read from file "filename" */
/* Supported dimensions: 2D, 3D */
/* NOTE: 1st line of given file must contain mesh->nFaces */
/* OUT: mesh->faces and mesh->nFaces set */
void read_faces ( char *filename, int nDim, int nVertsPerFace, int nFaces, int *faces );

/* Function to set all mesh (velocity) times read from file "filename" */
/* Supported dimensions: 2D, 3D */
/* NOTE: 1st line of given file must contain mesh->nTimes */
/* OUT: mesh->times and mesh->nTimes set */
void read_time ( char *filename, int nTimes, double *times );

/* Function to set all mesh points' velocity read from file "filename" */
/* Supported dimensions: 2D, 3D */
/* OUT: mesh->points[0...nPoints].velocity set */
void read_velocities ( char *filename, int nPoints, int nDim, int nTimes, double *coords, double *velocity );

/* Function to set the lists of faces in which each point is a vertex (and nFaces) */
/* Supported dimensions: 2D, 3D */
/* OUT: mesh->points[0...nPoints].faces and mesh->points[0...nPoints].nFaces set */
//void assign_faces_to_verts ( mesh_t *mesh );

void create_nFacesPerPoint_vector ( int nDim, int nPoints, int nFaces, int nVertsPerFace, int *faces, int *nFacesPerPoint );
void create_facesPerPoint_vector ( int nDim, int nPoints, int nFaces, int nVertsPerFace, int *faces, int *nFacesPerPoint, int *facesPerPoint );
