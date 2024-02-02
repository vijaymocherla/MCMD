#!/usr/bin/env python
#  
# 

import numpy as np
from itertools import product

def get_positions(N, rho, fill_type="random"):
    vol = N/rho
    l = vol**(1/3)
    box = np.array([l]*3)
    if fill_type == "random":
        x = random_fill(N, box)
    else:
        x = lattice_fill(N, box, fill_type)
    return x

def get_velocities(N, T, sample_type='uniform'):
    v = np.random.random((N,3))
    if (sample_type=='uniform'):
        v = v - 0.5
    elif(sample_type=='boltzmann'):
        #  sampling v in (-4*sigma, 4*sigma)
        v = 8*v 
        v = 1/(np.sqrt(2*np.pi)) * np.exp((v)**2 / 2)            
    else:
        raise NameError("sample type '{}' not recognised!".format(sample_type))
    v_com = np.sum(v, axis=0)/N
    v[:,0] =(v[:,0]-v_com[0])
    v[:,1] =(v[:,1]-v_com[1])
    v[:,2] =(v[:,2]-v_com[2])
    # Scaling velocities to given temperature
    alpha = np.sqrt((3*N*T)/(np.sum(v**2)))
    v = alpha*v
    return v

def random_fill(N, box):
    x = np.zeros((N,3))
    i = 0
    def recursive_fill(i):
        xi = box*np.random.random(3)
        for j in range(i):
            xr = (x[j,:] - xi)
            xr = xr - np.rint(xr/box)*box
            dr = np.sqrt(np.sum(xr**2))
            if (dr < 1.0):
                xi = recursive_fill(i)
        return xi
    while (i<N):
        xi = recursive_fill(i)
        x[i, :] = xi
        i += 1
    return x
    
def lattice_fill(N, box, lattice):
    rvecs = { 'cubic': np.array([[0.25, 0.25, 0.25]]),
              'bcc'  : np.array([[0.25, 0.25, 0.25],
                                 [0.75, 0.75, 0.75]]),
              'fcc'  : np.array([[0.25, 0.25, 0.25],
                                 [0.25, 0.75, 0.75],
                                 [0.75, 0.25, 0.75],
                                 [0.75, 0.75, 0.25]])}
    x = np.empty((N,3))
    if lattice =='cubic':
        n = int((N)**(1/3))
    elif lattice == 'fcc':
        n = int((N/4)**(1/3))+1
    elif lattice == 'bcc':
        n = int((N/2)**(1/3))+1
    else: 
        raise NameError("lattice type '{}' not recognised!".format(lattice)) 
    cell_size = box/n
    indices = product(range(n), repeat=3)
    for m, index in enumerate(indices):
        for r in rvecs[lattice]:
            x[m,:] = cell_size*(np.array(index) + r)
    if m == N:
        raise Exception("no. of particles error!")
    return x
    
def check_overlap(x, box):
    N = x.shape[0]
    overlap = False
    for i in range(N):
        for j in range(i+1,N):
            xr = x[i,:] - x[j,:]
            xr = xr - np.rint(xr/box)*box
            dr = np.sqrt(np.sum(xr**2))
            if (dr < 1.0):
                print("Particles {i}, {j} overlap".format(i=i,j=j))
                overlap = True
    if overlap:
        return 0
    else:
        return 1


# experimental 
# random fill algorithm using minimization/newton-raphson 

from scipy.optimize import minimize

def opt_random_fill(N, box):
    x0 = box*np.random.random((N,3))
    x0 = x0.reshape(N*3)
    def potential(x, rcut=2.50, ecut=0.016316891136):
        en = 0
        n = x.shape[0]
        n = int(n/3)
        x = x.reshape(n,3)    
        for i in range(n):
            for j in range(i+1,n):
                xr  = x[i,:] - x[j,:]
                xr = xr - np.array(xr/box)*box
                r2 = np.sum(xr**2)
                if (r2 < rcut**2):
                    r2i = 1/r2
                    en = 4.0*(r2i**(3)-r2i**(6)) - ecut
        return en
    res = minimize(potential, x0, method='BFGS')
    overlap = check_overlap(x, box)
    if (overlap == 0):
        raise Exception("Minimization Failed!")
    x = res.x.reshape(N,3)
    return x


def read_xyz(filename, alpha=3.405):
    file = open(filename)
    n = int(file.readline())
    comment = file.readline()
    x = np.empty((n,3))
    for i in range(n):
        line = file.readline().split()
        element = line[0]
        x[i,0] = float(line[1])
        x[i,1] = float(line[2])
        x[i,2] = float(line[3])
    x = x/alpha
    return x
