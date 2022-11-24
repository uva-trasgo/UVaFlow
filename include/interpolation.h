/* Find the index of the closest mesh point */
/* Supported dimensions: 2D, 3D */
/* OUT: mesh closest point index */
int reset_coordinates ( int nDim, int nPoints, double *coords, double *Pcoords );

/* Interpolate velocity in a point P inside a triangle defined by V1-V2-V3 */
/* Based on the barycentric technique */
/* Supported dimensions: 2D */
/* OUT: interpolated_vel (double *) */
void interpolate_triangle (   double *V1, double *velV1, 
                              double *V2, double *velV2, 
                              double *V3, double *velV3,
                              double *P,
                              double *interpolated_vel,
			      int nDim );

/* Interpolate velocity in a tetrahedral given its 4 vertices V1-V2-V3-V4 and a point P */
/* Based on the barycentric technique */
/* Supported dimensions: 3D */
/* OUT: interpolated_vel (double *) */
void interpolate_3D_tetrahedral ( double *V1, double *velV1,
                      double *V2, double *velV2,
                      double *V3, double *velV3, 
                      double *V4, double *velV4,
                      double *P,
                      double *interpolated_vel,
		      int nDim );

/* Linear interpolation main function, performs:  */
/* A) Classical linear interpolation if the P is a mesh point */
/* B) Barycentric technique based interpolation if P is inside a mesh face */
/* C) Classical linear interpolation of the reseted (to the closest mesh point) P if P is outside the mesh */
/* Supported dimensions: 2D, 3D */
/* OUT: interpolation (double *) */
//void linear_interpolation ( double t, double *Pcoords, mesh_t *mesh, double *interpolation );
//void linear_interpolation_approach1 ( double t, double *Pcoords, mesh_t *mesh, double *interpolation); //, int it );
void linear_interpolation_approach2 ( double t, double *Pcoords, double *times, double *interpolation, int nDim, int nPoints, int nTimes, int nVertsPerFace, int nFaces, int *faces, double *coords, double *velocities); //, int it );
