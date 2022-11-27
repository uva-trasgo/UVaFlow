/* Binary search of an element in an array (of double type) of length size       */
/* Returns: -1 (element not in array min-max range) | 0 (not found) | 1 (found ) */
/* When returning 0, itime stores the index in array                             */
/* When returning 1, itime stores the index of the closest higher array element  */
int binarySearch(double *array, int size, double element, int *itime);
