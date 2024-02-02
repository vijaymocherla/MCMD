#!/usr/bin/env python
#
#

import numpy as np

from initialise import get_positions, get_velocities, read_xyz

class moldyn:
    """ An NVE Molecular Dynamics Simulation
    """
    def __init__(self, N, rho, T, rcut=2.50):
        """
        """
        self.N = N
        self.rho = rho 
        self.box = np.array([(N/rho)**(1/3)]*3)
        self.T = T
        self.rcut = rcut
        self.ecut = 4.0*((1/rcut)**6 - (1/rcut)**12)
        self.x = np.empty((N,3), dtype=np.float64)
        self.v = np.empty((N,3), dtype=np.float64)
        self.f = np.empty((N,3), dtype=np.float64)
        
    def initialize(self, fill_type="random", sample_type="uniform"):
        try:
            self.x = get_positions(self.N, self.rho, fill_type)
            self.v = get_velocities(self.N, self.T, sample_type)
        except:
            raise RuntimeError("Failed to initialise NVE MD simulation!")
        return 1
    
    @staticmethod
    def force(x, box, rcut=2.50, ecut=0.016316891136):
        """
        """
        N = x.shape[0]
        for i in range(N):
            for j in range(i+1, N):
                xr = x[i,:] - x[j,:]
                rij = np.sqrt(sum(xr**2))
                if (rij < rcut):
                    f

        return 1

    def calc_potential(self):
        """
        """
        potential = 0.0
        return potential
    
    def verlet(ti, tf, dt):
        
        return 1

    
