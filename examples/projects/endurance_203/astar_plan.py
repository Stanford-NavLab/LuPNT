import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class AStarPlanner(object):

    def __init__(self, statespace_lo, statespace_hi, x_init, x_goal, grid_env, resolution=1, diag=False, obstacles=None, elev_weight = 10):

        self.statespace_lo = statespace_lo
        self.statespace_hi = statespace_hi
        # this is the "occupancy grid" that we are working with
        self.grid_env = grid_env
        self.grid = grid_env.get_grid()
        self.obstacles = obstacles
        self.resolution = resolution
        self.diag = diag

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

        # for path planning
        self.elev_weight = elev_weight

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

    # change this to cost-to-go
    def cost_to_go(self, x1, x2):
        """
        Computes the Euclidean distance between two states.
        Inputs:
            x1: First state tuple
            x2: Second state tuple
        Output:
            Float Euclidean distance
        """

        # print(f"x1: {x1}, x2: {x2}")

        # scale it according to resolution
        euc_dist = self.resolution*(np.linalg.norm(np.array(x2) - np.array(x1)))
        #take into account elevation of the grid
        elev1 = self.grid[x1[0], x1[1], 0, 0]
        elev2 = self.grid[x2[0], x2[1], 0, 0]
        # print('current cell')
        # print(elev1)
        # print('neigh cell')
        # print(elev2)

        elev_diff = self.elev_weight*np.abs(elev2 - elev1)
        # print(elev_diff)

        # print(f"euc_dist: {euc_dist}, elev_diff: {elev_diff}")

        # return elev_diff
        # print(elev_diff)
        return elev_diff
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

    def reconstruct_path(self):
        """
        Use the came_from map to reconstruct a path from the initial location to
        the goal location
        Output:
            A list of tuples, which is a list of the states that go from start to goal
        """
        path = [self.x_goal]
        current = path[-1]
        while current != self.x_init:
            path.append(self.came_from[current])
            current = path[-1]
        return list(reversed(path))

    def plot_path(self, fig_num=0, show_init_label=True):
        """Plots the path found in self.path and the obstacles"""
        if not self.path:
            return

        fig, ax = self.grid_env.plot_grid_elev()

        solution_path = np.asarray(self.path)

        plt.plot(solution_path[:,1],solution_path[:,0], color="red", linewidth=2, label="A* path", zorder=10)
        plt.scatter([self.x_init[1], self.x_goal[1]], [self.x_init[0], self.x_goal[0]], color="red", s=30, zorder=10)

        if show_init_label:
            plt.annotate(r"$x_{init}$", np.array([self.x_init[1], self.x_init[0]]) + np.array([-2, -2]), fontsize=16)
            
        plt.annotate(r"$x_{goal}$", np.array([self.x_goal[1], self.x_goal[0]]) + np.array([1, 1]), fontsize=16)
        plt.legend()
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

    # def plot_tree(self, point_size=15):
    #     plot_line_segments([(x, self.came_from[x]) for x in self.open_set if x != self.x_init], linewidth=1, color="blue", alpha=0.2)
    #     plot_line_segments([(x, self.came_from[x]) for x in self.closed_set if x != self.x_init], linewidth=1, color="blue", alpha=0.2)
    #     px = [x[0] for x in self.open_set | self.closed_set if x != self.x_init and x != self.x_goal]
    #     py = [x[1] for x in self.open_set | self.closed_set if x != self.x_init and x != self.x_goal]
    #     plt.scatter(px, py, color="blue", s=point_size, zorder=10, alpha=0.2)

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

        while len(self.open_set) > 0:
            # get the current state by choosing the one with the minimum f score in the open set
            min_val = min( [ self.est_cost_through[k] for k in self.open_set ] )
            x_current = [key for key in self.open_set if self.est_cost_through[key] == min_val][0]

            # check if the current state is the goal state --> if so, update the path and end the search
            if x_current == self.x_goal:
                # make sure to update the path
                self.path = self.reconstruct_path()
                return True

            # remove the current state from the queue and add it to the considered state set
            self.open_set.remove(x_current)
            self.closed_set.add(x_current)

            for x_neigh in self.get_neighbors(x_current):
                # if the neighbor has already been considered as part of the frontier queue, skip it
                # we don't need to check it again, we know what the minimum cost path to get to it is
                if x_neigh in self.closed_set:
                    continue
                
                # set the tentative cost to arrive --> 
                # C_tilde(q') = C(q) [cost to come] + C(q,q') [immediate cost to come]
                tent_cost_to_arrive = self.cost_to_arrive[x_current] + self.cost_to_go(x_current, x_neigh)

                # if the neighbor state is not already in the queue, add it
                if x_neigh not in self.open_set:
                    self.open_set.add(x_neigh)
                # if it was in the queue and the tentative cost to arrive is more 
                # than the cost-to-arrive previously estimated, skip to the next neighbor (not a worthy node to get to from current state)
                elif tent_cost_to_arrive > self.cost_to_arrive[x_neigh]:
                    continue
                # keep track of the parent of this node
                self.came_from[x_neigh] = x_current
                # in this case, the neighbor was not in the queue and we are getting to it for the first time
                # set the cost to come to the tentative cost C(q') = C_tilde(q') -- (g score)
                self.cost_to_arrive[x_neigh] = tent_cost_to_arrive
                # set the est cost from start to finish to C(q') + h(q') -- (f score)
                # for now, h(q') = distance(q', goal)
                self.est_cost_through[x_neigh] = tent_cost_to_arrive + self.h_cost(x_neigh, self.x_goal)
            
        return False
        # raise NotImplementedError("solve not implemented")
        ########## Code ends here ##########

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

        # get_occup = self.occupancy              # get the corresponding DetOccupancyGrid2D
        if self.obstacles is None:
            obstacle_check = True
        else:
            obstacle_check = self.is_occupied(x)  

        # check within bounds
        bound_check = True
        if x[0] < self.statespace_lo[0] or x[0] >= self.statespace_hi[0] or x[1] < self.statespace_lo[1] or x[1] >= self.statespace_hi[1]:
            bound_check = False

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
        waypoints = solution_path[::len(solution_path)//(num_points-2)]
        waypoints = np.vstack([waypoints, solution_path[-1]])
        # print(path_plan_red.shape)
        return waypoints