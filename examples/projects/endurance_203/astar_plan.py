import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
np.random.seed(4)

class AStarPlanner(object):

    def __init__(self, statespace_lo, statespace_hi, x_init, x_goal, grid_env, tspan, resolution=1, diag=False, obstacles=None, weights = [0,0]):

        self.statespace_lo = statespace_lo
        self.statespace_hi = statespace_hi
        # this is the "occupancy grid" that we are working with
        self.grid_env = grid_env
        self.grid = grid_env.get_grid()
        self.obstacles = obstacles
        self.resolution = resolution
        self.diag = diag

        self.rover_speed = 1                            # km/hr
        self.t = tspan[0]                               # time 
        self.tspan = tspan  
        self.t_max = tspan[-1]                          # max time 

        self.width = self.grid.shape[0]
        self.height = self.grid.shape[1]

        self.x_offset = x_init
        self.x_init = self.snap_to_grid(x_init)    # initial state
        self.x_goal = self.snap_to_grid(x_goal)    # goal state

        self.closed_set = set()    # the set containing the states that have been visited
        self.open_set = set()      # the set containing the states that are candidate for future exploration

        self.est_cost_through = {}  # dictionary of the estimated cost from start to goal passing through state (often called f score)
        self.cost_to_arrive = {}    # dictionary of the cost-to-arrive at state from start (often called g score)
        self.came_from = {}         # dictionary keeping track of each state's parent to reconstruct the path
        self.travel_time = {}      

        # for path planning
        self.dist_weight = weights[0]
        self.elev_weight = weights[1]
        self.pdop_weight = weights[2]

        self.open_set.add(self.x_init)
        self.cost_to_arrive[self.x_init] = 0
        self.est_cost_through[self.x_init] = self.h_cost(self.x_init,self.x_goal)

        self.path = None        # the final path as a list of states

    def snap_to_grid(self, x):
        """ Returns the closest point on a discrete state grid
        Input:
            x: tuple state
        Output:
            A tuple that represents the closest point to x on the discrete state grid
        """
        return (
            self.resolution * round((x[0] - self.x_offset[0]) / self.resolution) + self.x_offset[0],
            self.resolution * round((x[1] - self.x_offset[1]) / self.resolution) + self.x_offset[1],
        )

    def cost_to_go(self, x1, x2, t=0):
        # scale it according to resolution
        # euc_dist = np.linalg.norm(np.array(x2) - np.array(x1))
        # get the index of the states to relate to the grid
        x1_idx = [int(x1[0]/self.resolution), int(x1[1]/self.resolution)]
        x2_idx = [int(x2[0]/self.resolution), int(x2[1]/self.resolution)]
        #take into account elevation of the grid
        elev1 = self.grid[x1_idx[0], x1_idx[1], t, 0]
        elev2 = self.grid[x2_idx[0], x2_idx[1], t, 0]
        elev_diff = np.abs(elev2 - elev1)

        state1 = np.array([x1[0], x1[1], elev1])
        state2 = np.array([x2[0], x2[1], elev2])


        # take into account the pdop of the grid
        pdop1 = self.grid[x1_idx[0], x1_idx[1], t, 2]
        pdop2 = self.grid[x2_idx[0], x2_idx[1], t, 2]
        pdop_cost = (pdop2 - pdop1)*1e4

        # check for nan and ensure cost is never negative
        if np.isnan(pdop_cost):
            print('here')
            pdop_cost = 5
        else:
            # pdop_cost = min([pdop1, pdop2])*10000

            if pdop_cost >=0:
                # if the cost is positive, neighbor has higher PDOP, we don't want that
                pdop_cost = np.abs(pdop_cost)*10
            else:
                # if the cost is negative, neighbor has lower PDOP, we want that
                pdop_cost = np.abs(pdop_cost)

        # print(f'Distance cost {self.dist_weight*np.linalg.norm(state1-state2)}')
        # print(f'Elevation cost {self.elev_weight*(elev_diff)}')
        # print(f'PDOP cost {self.pdop_weight*pdop_cost}')

        return self.dist_weight*np.linalg.norm(state1-state2) + self.elev_weight*(elev_diff) + self.pdop_weight*pdop_cost
    
    def distance_traveled(self, x1, x2):
        """
        Computes the Euclidean distance between two states.
        Inputs:
            x1: First state tuple
            x2: Second state tuple
        Output:
            Float Euclidean distance
        """
        # scale it according to resolution
        # euc_dist = np.linalg.norm(np.array(x2) - np.array(x1))
        x1_idx = [int(x1[0]/self.resolution), int(x1[1]/self.resolution)]
        x2_idx = [int(x2[0]/self.resolution), int(x2[1]/self.resolution)]
        #take into account elevation of the grid
        elev1 = self.grid[x1_idx[0], x1_idx[1], 0, 0]
        elev2 = self.grid[x2_idx[0], x2_idx[1], 0, 0]
        # elev_diff = self.elev_weight*np.abs(elev2 - elev1)
        state1 = np.array([x1[0], x1[1], elev1])
        state2 = np.array([x2[0], x2[1], elev2])

        return np.linalg.norm(state1-state2)
        # return np.sqrt(euc_dist**2 + elev_diff**2)

    def h_cost(self, x1, x2):
        # scale it according to resolution
        euc_dist = self.resolution*(np.linalg.norm(np.array(x2) - np.array(x1)))
        return euc_dist

    def get_neighbors(self, x):
        """
        Gets the FREE neighbor states of a given state x. Assumes a motion model
        where we can move up, down, left, right, or along the diagonals by an
        amount equal to self.resolution.
        Input:
            x: tuple state
        Ouput:
            List of neighbors that are free, as a list of TUPLES
        """
        neighbors = []
        # get the distance that you can travel with one step
        step = self.resolution
        if self.diag:
            # go up, down, left, right
            neighbors = neighbors + [ (x[0], x[1] + step) , (x[0], x[1] - step), (x[0]-step,x[1]), (x[0] + step,x[1]) ]
            # go diagonal as well (the the cell right along the diagonal)
            neighbors = neighbors + [ (x[0] + step, x[1] + step) , (x[0] + step, x[1] - step),
                                    (x[0] - step, x[1] - step) , (x[0] - step, x[1] + step) ]
        else:
            neighbors = neighbors + [ (x[0], x[1] + step) , (x[0], x[1] - step), (x[0]-step,x[1]), (x[0] + step,x[1]) ]
            
        # go through the neighbors, snap to grid, and check if they are free
        for i in range(len(neighbors)):
            neighbors[i] = self.snap_to_grid(neighbors[i])

        
        # remove the neighbors that are not free (if there are obstacles in the way)
        # TODO: perhaps is_free becomes is_visible!!!!!! as well as checking bounds
        neighbors_copy = neighbors.copy()
        for x_neigh in neighbors_copy:
            if not self.is_free(x_neigh):
                neighbors.remove(x_neigh)
        
        return neighbors

    def find_best_est_cost_through(self):
        """
        Gets the state in open_set that has the lowest est_cost_through
        Output: A tuple, the state found in open_set that has the lowest est_cost_through
        """
        return min(self.open_set, key=lambda x: self.est_cost_through[x])

    def reconstruct_path(self, goal):
        """
        Use the came_from map to reconstruct a path from the initial location to
        the goal location
        Output:
            A list of tuples, which is a list of the states that go from start to goal
        """
        path = [goal]
        current = path[-1]
        while current != self.x_init:
            path.append(self.came_from[current])
            current = path[-1]
        return list(reversed(path))

    def reconstruct_time(self, path):
        time = []
        for x in path:
            time.append(self.travel_time[x])
        return np.array(time)

    def plot_path(self, fig_num=0, show_init_label=True):
        """Plots the path found in self.path and the obstacles"""
        if not self.path:
            return

        fig, ax = self.grid_env.plot_grid_elev()

        solution_path = np.asarray(self.path)/self.resolution
        xi_plot = [self.x_init[1]/self.resolution, self.x_goal[1]/self.resolution]
        yi_plot = [self.x_init[0]/self.resolution, self.x_goal[0]/self.resolution]

        plt.plot(solution_path[:,1],solution_path[:,0], color="red", linewidth=2, label="A* path", zorder=10)
        plt.scatter(xi_plot, yi_plot, color="red", s=30, zorder=10)

        if show_init_label:
            plt.annotate(r"$q_{init}$", np.array([self.x_init[1], self.x_init[0]])/self.resolution + np.array([-2, -2]), fontsize=16)
            
        plt.annotate(r"$q_{goal}$", np.array([self.x_goal[1], self.x_goal[0]])/self.resolution + np.array([1, 1]), fontsize=16)
        plt.legend(loc = 'upper left')
        # plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, ncol=3)

        plt.axis([0, self.grid.shape[0], 0,self.grid.shape[1]])
        # plt.tight_layout()

        if self.obstacles is not None:
            for obs in self.obstacles:
                ax.add_patch(
                patches.Rectangle(
                obs[0],
                obs[1][0]-obs[0][0],
                obs[1][1]-obs[0][1], color = "tab:gray"))
        return fig, ax

    def get_travel_time(self, x1, x2):
        # get the distance between the two states
        dist = self.distance_traveled(x1, x2)
        # get the time it takes to travel that distance
        time = dist/self.rover_speed*3600
        return time
    
    def solve(self):
        """
        Solves the planning problem using the A* search algorithm. It places
        the solution as a list of tuples (each representing a state) that go
        from self.x_init to self.x_goal inside the variable self.path
        Input:
            None
        Output:
            Boolean, True if a solution from x_init to x_goal was found

        HINTS:  We're representing the open and closed sets using python's built-in
                set() class. This allows easily adding and removing items using
                .add(item) and .remove(item) respectively, as well as checking for
                set membership efficiently using the syntax "if item in set".
        """
        ########## Code starts here ##########

        # Initialization is taken care of, this is what it looks like

        ## Initialize the open set with x_init
        # self.open_set.add(self.x_init)
        ## add to cost to arrive ( C(q_I), cost-to-come)
        # self.cost_to_arrive[self.x_init] = 0
        ## add to cost thru ( the estimated cost to get to the goal, cost-to-go )
        # self.est_cost_through[self.x_init] = self.distance(self.x_init, self.x_goal)

        self.travel_time[self.x_init] = 0

        while len(self.open_set) > 0:
            # get the current state by choosing the one with the minimum f score in the open set
            min_val = min( [ self.est_cost_through[k] for k in self.open_set ] )
            x_current = [key for key in self.open_set if self.est_cost_through[key] == min_val][0]

            # check if the current state is the goal state --> if so, update the path and end the search
            if x_current == self.x_goal:
                # make sure to update the path
                self.path = self.reconstruct_path(self.x_goal)
                return True

            # the grid evolves over time, so we need to check if we have visibility or not
            current_time = self.travel_time[x_current]
            # print(current_time/3600)
            # x_idx = [int(x_current[0]/self.resolution), int(x_current[1]/self.resolution)]
            # print(x_idx)

            if current_time > self.t_max:
                print("Time exceeded")
                # find out the closest point to the goal
                dist_from_goal = []
                for x in self.closed_set:
                    dist_from_goal.append(self.h_cost(x, self.x_goal))
                x_current = [key for key in self.closed_set if self.h_cost(key, self.x_goal) == min(dist_from_goal)][0]
                self.path = self.reconstruct_path(x_current)
                return True

            # print(current_time/3600)
            t_idx = (np.abs(self.tspan - (current_time))).argmin()
            pdop_check = self.pdop_check(x_current, t_idx)
            # print(pdop_check)
            if not pdop_check and self.pdop_weight > 0:
                # we do not have any satellites in view, so we are going to stay still for 10 minutes
                print("No satellites in view")
                self.travel_time[x_current] += (60*10)
                continue

            # remove the current state from the queue and add it to the considered state set
            self.open_set.remove(x_current)
            self.closed_set.add(x_current)

            # search tree
            for x_neigh in self.get_neighbors(x_current):

                # ### NEW --> get how long it takes to get to the neighbor based on rover speed
                time_to_get_there = self.get_travel_time(x_current, x_neigh)
                t_idx = (np.abs(self.tspan - (current_time+time_to_get_there))).argmin()
                # print(t_idx)

                # if the neighbor has already been considered as part of the frontier queue, skip it
                # we don't need to check it again, we know what the minimum cost path to get to it is
                if x_neigh in self.closed_set:
                    continue
                
                # set the tentative cost to arrive --> 
                # C_tilde(q') = C(q) [cost to come] + C(q,q') [immediate cost to come]
                tent_cost_to_arrive = self.cost_to_arrive[x_current] + self.cost_to_go(x_current, x_neigh, t_idx)

                # if the neighbor state is not already in the queue, add it
                if x_neigh not in self.open_set:
                    self.open_set.add(x_neigh)
                # if it was in the queue and the tentative cost to arrive is more 
                # than the cost-to-arrive previously estimated, skip to the next neighbor (not a worthy node to get to from current state)
                elif tent_cost_to_arrive > self.cost_to_arrive[x_neigh]:
                    continue
                # keep track of the parent of this node
                self.came_from[x_neigh] = x_current
                self.travel_time[x_neigh] = time_to_get_there + current_time
                # in this case, the neighbor was not in the queue and we are getting to it for the first time
                # set the cost to come to the tentative cost C(q') = C_tilde(q') -- (g score)
                self.cost_to_arrive[x_neigh] = tent_cost_to_arrive
                # set the est cost from start to finish to C(q') + h(q') -- (f score)
                # for now, h(q') = distance(q', goal)
                self.est_cost_through[x_neigh] = tent_cost_to_arrive + self.h_cost(x_neigh, self.x_goal)


        return False
        # raise NotImplementedError("solve not implemented")
        ########## Code ends here ##########

    def pdop_check(self, x, t):
        """check if we have enough satellites in view"""
        x_idx = [int(x[0]/self.resolution), int(x[1]/self.resolution)]
        if np.isnan(self.grid[x_idx[0], x_idx[1], t, 2]):
            return False
        else:
            return True

    def is_free(self, x):
        """
        Checks if a give state x is free, meaning it is inside the bounds of the map and
        is not inside any obstacle.
        Inputs:
            x: state tuple
        Output:
            Boolean True/False
        Hint: self.occupancy is a DetOccupancyGrid2D object, take a look at its methods for what might be
              useful here
        """
        
        # Note to self: yes, this is where you check for collisions AND being outside bounds
        # since there are no obstacles, we just need to check if the point is within the bounds

        x_idx = [int(x[0]/self.resolution), int(x[1]/self.resolution)]

        # get_occup = self.occupancy              # get the corresponding DetOccupancyGrid2D
        if self.obstacles is None:
            obstacle_check = True
        else:
            obstacle_check = self.is_occupied(x)  

        # check within bounds
        bound_check = True
        if x_idx[0] < self.statespace_lo[0] or x_idx[0] >= self.statespace_hi[0] or x_idx[1] < self.statespace_lo[1] or x_idx[1] >= self.statespace_hi[1]:
            bound_check = False

        # pdop_check = True
        # if np.isnan(self.grid[x_idx[0], x_idx[1], 0, 2]):
        #     pdop_check = False

        if bound_check and obstacle_check:
            return True
        else:
            return False
        
    def is_occupied(self, x):
        """Verifies that point is not inside any obstacles by some margin"""
        for obs in self.obstacles:
            if x[0] >= obs[0][0] - self.width * .01 and \
               x[0] <= obs[1][0] + self.width * .01 and \
               x[1] >= obs[0][1] - self.height * .01 and \
               x[1] <= obs[1][1] + self.height * .01:
                return False
        return True

    def get_waypoints(self, num_points):
        """Returns the waypoints that define the path"""
        solution_path = np.asarray(self.path)
        print(solution_path.shape)
        index = np.linspace(0, solution_path.shape[0] - 1, num_points, dtype=int)
        waypoints = solution_path[index, :]
        # waypoints = solution_path[::len(solution_path)//(num_points-2)]
        # print(waypoints)
        # waypoints = np.vstack([waypoints, solution_path[-1]])
        # print(path_plan_red.shape)
        # add the desired heading to the waypoints
        # headings = np.zeros((waypoints.shape[0], 1))
        # stack the headings
        # waypoints = np.hstack([waypoints, headings])

        return waypoints