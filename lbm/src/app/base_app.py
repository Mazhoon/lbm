# Generic imports
import math

# Custom imports
from lbm.src.utils.shapes import *
from lbm.src.core.obstacle import *

###############################################
### Base app
class base_app():
    ### Set inlet fields
    def set_inlets(self, lattice, it):

        pass

    ### Compute observables
    def observables(self, lattice, it):

        pass

    ### Finalize
    def finalize(self, lattice):

        pass

    def add_obstacle(self, lattice, obs, tag):
        shape = generate_shape(obs.n_pts, obs.pos,
                                obs.type,  obs.size,
                                obs.name,  obs.n_spts,
                                lattice.output_dir)
        obs.set_polygon(shape.curve_pts)
        obs.set_tag(tag)
        area, bnd, ibb = lattice.add_obstacle(obs)
        obs.fill(area, bnd, ibb)
        
    ### Add obstacles
    def add_obstacles(self, lattice, obstacles):

        for i in range(len(obstacles)):
            self.add_obstacle(lattice, obstacles[i], (i+1))

    ### Iteration printings
    def printings(self, it):

        if (self.stop == 'it'):
            print('# it = '+str(it)+' / '+str(self.it_max), end='\r')
        if (self.stop == 'obs'):
            str_d  = "{:10.6f}".format(self.drag_buff.obs)
            str_l  = "{:10.6f}".format(self.lift_buff.obs)

            print('# it = '+str(it)+
                  ', avg drag ='+str_d+', avg lift ='+str_l, end='\r')

    ### Check stopping criterion
    def check_stop(self, it):

        compute = True

        if (self.stop == 'it'):
            if (it >= self.it_max):
                compute = False
                print('\n')
                print('# Computation ended: it>it_max')

        if (self.stop == 'obs'):
            if (self.drag_buff.obs_cv and self.lift_buff.obs_cv):
                compute = False
                print('\n')
                print('# Computation ended: converged')

        return compute
