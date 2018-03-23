from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt

def mesh_setup(npts,ifPlot):

    # define domain boundaries
    xl = yl = 0
    xr = yl = 1

    # compute 'npts' random points in each direction
    points = np.random.rand(2,npts)
    
    # compute the delauney triangulation 
    tess = Delaunay(points)

    # plot if requested
    if ifPlot is False:
        plt.triplot(points[:,0],points[:,1],tess.simplices.copy())
        plt.plot(points[:,0],points[:,1],'o')
        plt.show()
    
    return tess
