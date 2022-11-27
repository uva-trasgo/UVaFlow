/* Find the index of the closest mesh point */
/* Supported dimensions: 2D, 3D */
/* OUT: mesh closest point index */
int reset_coordinates_2D ( int nDim, int nPoints, double *coords_x, double *coords_y, double *Pcoords );
int reset_coordinates_3D ( int nDim, int nPoints, double *coords_x, double *coords_y, double *coords_z, double *Pcoords );

/* Interpolate velocity in a point P inside a triangle defined by V1-V2-V3 */
/* Based on the barycentric technique */
/* Supported dimensions: 2D */
/* OUT: interpolated_vel (double *) */
void interpolate_triangle (   double *coords_x, double *coords_y, double *velocities,
                              int iv1, int iv2, int iv3, int itime, 
                              int nPoints, double *P, double *interpolated_vel, int nDim );

/* Interpolate velocity in a tetrahedral given its 4 vertices V1-V2-V3-V4 and a point P */
/* Based on the barycentric technique */
/* Supported dimensions: 3D */
/* OUT: interpolated_vel (double *) */
void interpolate_3D_tetrahedral ( double *coords_x, double *coords_y, double *coords_z, double *velocities,
				  int iv1, int iv2, int iv3, int iv4, int itime, 
				  int nPoints, double *P, double *interpolated_vel, int nDim );

/* Linear interpolation main function, performs:  */
/* A) Classical linear interpolation if the P is a mesh point */
/* B) Barycentric technique based interpolation if P is inside a mesh face */
/* C) Classical linear interpolation of the reseted (to the closest mesh point) P if P is outside the mesh */
/* Supported dimensions: 2D, 3D */
/* OUT: interpolation (double *) */
void linear_interpolation_approach2_2D ( double t, double *Pcoords, double *times, double *interpolation, 
					 int nDim, int nPoints, int nTimes, int nVertsPerFace, int nFaces, 
					 int *faces, double *coords_x, double *coords_y, double *velocities);

void linear_interpolation_approach2_3D ( double t, double *Pcoords, double *times, double *interpolation, 
					 int nDim, int nPoints, int nTimes, int nVertsPerFace, int nFaces, 
					 int *faces, double *coords_x, double *coords_y, double *coords_z, double *velocities);
