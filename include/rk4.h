/* Function to apply Runge-Kutta (4th order) method */
/* Supported dimensions: 2D, 3D */
/* OUT: result containing the estimation for 2D/3D */
void runge_kutta_4 ( 	double Pcoords_x, double Pcoords_y, double Pcoords_z, double t0, double tend, 
			double *result, int nsteps, int nDim, int nPoints, int nTimes, 
			double *times, int nVertsPerFace, int nFaces, int *faces, 
			double *coords_x, double *coords_y, double *coords_z, double *velocities, void *kdtree,
			int *nFacesPerPoint, int *facesPerPoint);
