#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "search.h"

int binarySearch(double *array, int size, double element, int *itime)
{
    int low  = 0;
    int high = size - 1;
    int middle;
    
    if ( ( array[0] > element ) || ( array[size-1] < element ) )
        return -1;

    int res = 0;
    *itime = array[size-1];

    while (high - low > 1) 
    {
        middle = (high + low) / 2;
        if (array[middle] < element) 
        {
            low = middle;
        }
        else 
        {
            *itime = middle;
            high = middle;
        }
    }
    if (array[low]  == element) 
    {
        *itime = low;
        res = 1;
    }
    else
    {
        if (array[high] == element)
        {
            *itime = high;
            res = 1;
        }
    }
    return res;
}
