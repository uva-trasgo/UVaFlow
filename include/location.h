/* Check if a given point is inside a triangle described by V1, V2, V3 vertices */
/* Supported dimensions: 2D */
/* OUT: int 1 (yes) / int 0 (no) */
int p_inside_triangle ( int iv1, int iv2, int iv3, double *coords_x, double *coords_y, double *P, int nDim );

/* Check if a given point is inside a mesh face */
/* Supported dimensions: 2D */
/* OUT: point index in mesh->points array (-1 if the given point is not the mesh) */
int find_point_simplex_in_mesh_2D ( double *Pcoords, int nDim, int nFaces, int nVertsPerFace, double *coords_x, double *coords_y, int *faces );
int find_point_simplex_in_mesh_3D ( double *Pcoords, int nDim, int nFaces, int nVertsPerFace, double *coords_x, double *coords_y, double *coords_z, int *faces );

/* Check if a given point is or not a mesh point */
/* Supported dimensions: 2D, 3D */
/* OUT: point index in mesh->points array (-1 if the given point is not the mesh) */
int find_point_in_mesh_2D ( double *Pcoords, int nDim, int nPoints, double *coords_x, double *coords_y );
int find_point_in_mesh_3D ( double *Pcoords, int nDim, int nPoints, double *coords_x, double *coords_y, double *coords_z );

/* Check if a given point is inside a tetrahedron described by V1, V2, V3, V4 vertices */
/* share_side is a support function */
/* Supported dimensions: 3D */
/* OUT: int 1 (yes) / int 0 (no) */
int p_inside_tetrahedron ( int iv1, int iv2, int iv3, int iv4, double *coords_x, double *coords_y, double *coords_z, double *P, int nDim );
int share_side ( double *V1, double *V2, double *V3, double *V4, double *P, int nDim );
