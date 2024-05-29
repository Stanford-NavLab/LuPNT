import os
import jax
import jax.numpy as jnp
import numpy as np
import matplotlib.pyplot as plt
import pylupnt as pnt

class GridWorld:
    def __init__(self):
        # constants
        self.N = 100                                    # grid size
        self.params = 3                                 # number of parameters for each cell
        self.grid = self.create_grid()                  # grid
        self.moon_pa_origin = np.array([0, 0, 0])       # moon PA origin (to change)


    def create_grid(self):
        """
        Returns a 2D grid of size N x N where each cell is a tuple (x,y) representing the coordinates of the cell.
        Output:
            grid: N x N numpy array of tuples
        """
        grid = np.zeros((self.N, self.N, self.params))
        return grid
    
    def create_crater(self, radius, max_depth, loc = None) :
        if loc is not None:
            crater = np.array([loc[0], loc[1], radius, -max_depth])
        else:
            crater = np.array([np.random.randint(0, self.N), np.random.randint(0, self.N), radius, -max_depth])
        return crater
    
    def add_elevation(self, loc, param_idx, elevation):
        # add the elevation to a grid cell
        self.grid[loc[0], loc[1], param_idx] = elevation
        # return self.grid
    
    def add_crater(self, crater, slope_factor):
        # add the crater to the grid
        # make it sloped! linear for now
        crater_r0 = crater[2]
        if slope_factor > 1:
            raise ValueError("Slope factor must be less than 1")
        elif slope_factor == 1:
            # no slope
            crater_ri = crater_r0
        else:
            crater_ri = crater_r0*slope_factor
            m = (crater[3])/(crater_ri - crater_r0)
            b = -m*crater_r0

        for i in range(self.N):
            for j in range(self.N):
                if (i - crater[0])**2 + (j - crater[1])**2 < crater_r0**2:
                    # we are inside the crater
                    if (i - crater[0])**2 + (j - crater[1])**2 < crater_ri**2:
                        # we are very much inside the crater
                        self.grid[i, j, 0] = crater[3] + self.grid[i, j, 0]
                    else:
                        # we are going down the crater
                        # radial distance from the center of the crater
                        rad_loc = np.sqrt((i - crater[0])**2 + (j - crater[1])**2)
                        # calculate the elevation
                        elevation = m*rad_loc + b
                        self.grid[i, j, 0] = elevation + self.grid[i, j, 0]

    def plot_grid(self, param_idx = 0):
        fig, ax = plt.subplots(dpi=150)
        heatmap = ax.pcolor(self.grid[:,:,0])
        plt.colorbar(heatmap)
        plt.axis('equal')
        plt.xlim(0, self.N)
        plt.ylim(0, self.N)
        plt.show()