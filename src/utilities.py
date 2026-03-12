
import numpy as np

def curvature(x, y):
    """
    Calculates the 2D radisu of curviture for list of points x & Y

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    curvature = np.zeros_like(x)
    for i in range(1,len(x)-1):
        r = circumradius(x[i-1:i+2], y[i-1:i+2])
        
        # If the points are colinear then we assume radius is infinity
        if r == 0:
            curvature[i] = np.finfo(np.float64).tiny
        else:
            curvature[i] = 1.0 / r
    return curvature

def circumradius(x, y):
    """
    Calculate the circumradius given 3 points
    """
    den = 2*((x[1]-x[0])*(y[2]-y[1]) - (y[1]-y[0])*(x[2]-x[1]))
    num = np.sqrt( ((x[1]-x[0])**2 + (y[1]-y[0])**2) * ((x[2]-x[1])**2+(y[2]-y[1])**2) * ((x[0]-x[2])**2 + (y[0]-y[2])**2) )
    if den == 0:
        return np.finfo(np.float64).max
    return np.abs(num/den)

