import os
import jax
import jax.numpy as jnp
import numpy as np
import matplotlib.pyplot as plt
import pylupnt as pnt
import scipy as sc

class GridWorld:
    def __init__(self, N = 100, res = 1, Nt = 3600, moon_pa_origin = np.array([0, 0, 0]), mu_noise = 0.0, sigma_noise = 1.0):
        # constants
        self.N = N                                      # grid size
        # TODO: change this up
        self.params = 3                                 # number of parameters for each cell
        self.params_names = ['elevation', 'sat_visibility', 'PDOP']
        self.moon_pa_origin = moon_pa_origin            # moon PA origin (to change)
        self.res = res                                  # resolution of the grid (KILOMETERS)
        self.lent = Nt                                  # length of simulation (s)

        # create the grid (which is time dependent, we will see how that actually plays out)
        self.grid = self.create_grid(mu_noise, sigma_noise)                  # grid

    def create_grid(self, mu_noise, sigma_noise):
        grid = np.zeros((self.N, self.N, self.lent, self.params))
        # TODO: add noise to the elevation of the grid
        elevation_noise = np.random.normal(mu_noise, sigma_noise, size=(self.N,self.N))
        elevation_noise = elevation_noise.reshape(self.N, self.N, 1)
        elevation_noise = np.tile(elevation_noise, (1, 1, self.lent))
        grid[:, :, :, 0] += elevation_noise
        return grid
    
    def get_grid(self):
        # provide the grid, perhaps to populate with satellite metrics
        return self.grid
    
    #TODO --> make sure that the rover movement is REPRESENTATIVE of EAST/NORTH movement!!!

    def grid_PA_coords(self, lat_0, lon_0):
        # lets assume that x is E and y is N
        # converting origin from PA to ENU
        R_Moon = 1737.4  # km
        M_mcmf_2enu = mcmf_to_enu(lat_0, lon_0)
        # print(M_mcmf_2enu)
        origin_enu = M_mcmf_2enu @ self.moon_pa_origin
        # print(self.moon_pa_origin)
        # print(origin_enu)

        # store all the grid coordinates
        grid_coords = np.zeros((self.N, self.N, 3))
        grid_lat_long = np.zeros((self.N, self.N, 2))

        for i in range(self.N):
            for j in range(self.N):
                # move in east and north in the Moon PA frame
                delta_east = i*self.res
                delta_north = j*self.res

                # get change in latitude and longitude
                delta_lat = delta_north/R_Moon
                delta_long = delta_east/R_Moon

                lat = lat_0 + delta_lat; long = lon_0 + delta_long
                cell_elevation = self.get_elevation(i, j)
                
                grid_coords[i, j, :] = latlong_to_MoonPA(lat, long, cell_elevation) # N x N x 3
                grid_lat_long[i, j] = [lat, long] # N x N x 2
        

        return grid_coords, grid_lat_long

    
    def create_crater(self, radius, max_depth, loc = None) :
        if loc is not None:
            crater = np.array([loc[0], loc[1], radius, -max_depth])
        else:
            crater = np.array([np.random.randint(0, self.N), np.random.randint(0, self.N), radius, -max_depth])
        return crater
    
    def add_elevation(self, loc, elevation, param_idx=0):
        # add the elevation to a grid cell
        self.grid[loc[0], loc[1], :, param_idx] = elevation
        # return self.grid
        
    def add_PDOP(self, loc, param_idx, PDOP):
        # add the PDOP to a grid cell
        self.grid[loc[0], loc[1], :, param_idx] = PDOP
        
    def get_elevation(self, x, y):
        elevation = self.grid[:, :, 0, 0] # of size 40 x 40 x time length
        # elevation = elevation[:, :, 0] # of size 40 x 40
        return elevation[x, y]
    
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
                        self.grid[i, j, :, 0] = crater[3] + self.grid[i, j, :, 0]
                    else:
                        # we are going down the crater
                        # radial distance from the center of the crater
                        rad_loc = np.sqrt((i - crater[0])**2 + (j - crater[1])**2)
                        # calculate the elevation
                        elevation = m*rad_loc + b
                        self.grid[i, j, :, 0] = elevation + self.grid[i, j, :, 0]

    def plot_grid_elev(self, t = 0, param_idx = 0):
        fig, ax = plt.subplots(dpi=150)
        heatmap = ax.pcolor(self.grid[:,:, t, param_idx])
        cbar = plt.colorbar(heatmap,fraction=0.046, pad=0.04)
        plt.axis('equal')
        ax.set_xlim(0, self.N)
        ax.set_ylim(0, self.N)
        ax.set_xlabel('y [km]')
        ax.set_ylabel('x [km]')

        if param_idx == 0:
            cbar.set_label('Elevation [km]')
        elif param_idx == 2:
            cbar.set_label('PDOP')
        # plt.show()

        return fig, ax

    def generate_deterministic_trajectory(self, state_0, state_f, p, N = 10):
        # A function that generates trajectories as quadratic approximations of the p-function
        x0, yf, t0 = state_0
        xf, y0, tf = state_f

        eqn = lambda x : ((1**p - (x/xf)**p)**(1/p))*y0 - ((y0/xf)*x)
        alpha = sc.optimize.fsolve(eqn, xf/2)[0]

        A = np.zeros((6,6))
        poly_basis = lambda t_i :[1, t_i, t_i**2]

        A[0,0:3] = poly_basis(t0); A[1,3:] = poly_basis(t0);
        A[2,0:3] = poly_basis(N); A[3,3:] = poly_basis(N);
        A[4,0:3] = poly_basis(N/2); A[5,3:] = poly_basis(N/2);

        q = np.array([x0, y0, xf, yf, alpha, ((y0/xf)*alpha)])
        coeffs = np.linalg.solve(A, q)
        coeffs_x = coeffs[0:3]
        coeffs_y = coeffs[3:]

        # obtain the trajectory in terms of time
        traj = np.zeros((N, 2))
        t = np.linspace(t0, N, N)
        for i in range(len(t)):
            x_t = np.dot(coeffs_x, poly_basis(t[i]))
            y_t = np.dot(coeffs_y, poly_basis(t[i]))

            traj[i,:] = [x_t, y_t]

        traj[:,1] = y0-traj[:,1]

        return coeffs, traj
    
    def discretize_traj(self, traj):
        # linear interpolation between points
        x = traj[:,0]
        y = traj[:,1]
        traj_length = len(x)

        traj_discrete_x = [x[0]]
        traj_discrete_y = [y[0]]
        
        for i in range(traj_length-1):
            x0 = int(x[i])
            y0 = int(y[i])
            x1 = int(x[i+1])
            y1 = int(y[i+1])
            x_vals, y_vals = self.Bresenham_algorithm(x0, y0, x1, y1)
            traj_discrete_x += x_vals[1:]
            traj_discrete_y += y_vals[1:]

        traj_discrete = np.array([traj_discrete_x, traj_discrete_y]).T

        # bring back to the grid any values that are outside the grid
        traj_discrete[traj_discrete < 0] = 0
        traj_discrete[traj_discrete >= self.N-1] = self.N - 1
        
        # A function that discretizes the trajectory
        return traj_discrete
    
    def Bresenham_algorithm(self, x0, y0, x1, y1):
        # A function that discretizes a line segment
        # from Bresenham-line-drawing-algorithm GitHub Repo
        dx = x1 - x0
        dy = y1 - y0
        absdx = np.abs(dx)
        absdy = np.abs(dy)

        x = [x0]; y = [y0]

        if absdx > absdy:
            p = 2*absdy - absdx
            while x0 != x1:
                x0 += np.sign(dx)
                if p < 0:
                    p += 2*absdy
                else:
                    p += 2*absdy - 2*absdx
                    y0 += np.sign(dy)

                # we do not want to go diagonal
                if np.abs(x0 - x[-1]) == 1 and np.abs(y0-y[-1]) == 1:

                    x.append(x0)
                    y.append(y[-1])
                    
                    x.append(x0)
                    y.append(y0)

                else:
                    x.append(x0)
                    y.append(y0)

        else:
            p = 2*absdx - absdy
            while y0 != y1:
                y0 += np.sign(dy)
                if p < 0:
                    p += 2*absdx
                else:
                    p += 2*absdx - 2*absdy
                    x0 += np.sign(dx)
                    
                # we do not want to go diagonal
                if np.abs(x0 - x[-1]) == 1 and np.abs(y0-y[-1]) == 1:
                    x.append(x0)
                    y.append(y[-1])
                     
                    x.append(x0)
                    y.append(y0)

                else:
                    x.append(x0)
                    y.append(y0)

        return x, y

def mcmf_to_enu(lat_user, lon_user):
    # convert from MCMF to ENU coordinates
    # lat_user = rover latitude in radians (1x1)
    # lon_user = rover longitude in radians (1x1)
    phi, lmda = lat_user, lon_user
    M = np.array([
        [-np.sin(lmda), np.cos(lmda), 0],
        [-np.cos(lmda)*np.sin(phi), -np.sin(lmda)*np.sin(phi), np.cos(phi)],
        [np.cos(lmda)*np.cos(phi), np.sin(lmda)*np.cos(phi), np.sin(phi)]
    ])
    return M

def latlong_to_MoonPA(lat, long, cell_elevation):
    # convert latitude and longitude to MoonPA
    # cell_elevation is a float
    R_Moon = 1737.4  # km
    user_mcmf_pos_x = (R_Moon + cell_elevation) * np.cos(lat) * np.cos(long)
    user_mcmf_pos_y = (R_Moon + cell_elevation) * np.cos(lat) * np.sin(long)
    user_mcmf_pos_z = (R_Moon + cell_elevation) * np.sin(lat)

    # 3x1 position vector
    user_mcmf_loc = np.array([user_mcmf_pos_x, user_mcmf_pos_y, user_mcmf_pos_z])
    return user_mcmf_loc