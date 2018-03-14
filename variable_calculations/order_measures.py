import numpy as np
from copy import copy
from scipy.ndimage.interpolation import shift

'''
    series of functions which will extract all major order parameters from the 
    Ising simulation
'''

def magnetization(lattice):
    N_s = np.prod(lattice.shape)
    return abs((1/N_s) * np.sum(lattice));

def magnetization_2(grid):
    ''' calculate <m^2>'''
    N_s = np.prod(grid.shape);
    return ((1/N_s))*np.sum(grid*grid); #if we do grid*grid, this is some constant number...

def energy(grid, J=1):
    '''
    this function must be able to generallize to arbitrary dimensions
    :param grid:
    :param J:
    :return:
    '''
    N = grid.shape;
    dimension = len(N); #dimension of system being simulated
    NN = copy(grid);
    E = 0;
    neighbors = 0;

    for i in range(dimension):
        for j in [-1,1]:
            neighbors += np.roll(NN, shift=j, axis = i);
            E+=J*np.sum(grid*np.roll(NN, shift=j, axis = i));
    DeltaE = J * (grid* neighbors)/(np.prod(N));
    return  -np.sum(DeltaE)/2;#-(E/np.prod(N))/2; #return is avg energy per site

def energy_2(grid,J=1):
    ''' measurement of energy^2 per site, which should be constant'''
    N_s = np.prod(grid.shape);
    N = grid.shape;
    dimension = len(N); #dimension of system being simulated
    NN = copy(grid);
    E = 0;
    for i in range(dimension):
        for j in [-1,1]:
            E+=J*(grid*np.roll(NN, shift=j, axis = i));
    return np.sum(E*E)/N_s/2;

def s0sr(grid, r):
    '''
    :param grid:
    :param r:
    :return:
    '''
    '''there is no special point so just pick one origin for o'''
    N_s = grid.shape;
    s0 = grid[0,0];
    ssr = 0;
    for i in range(N_s[0]):
        for j in range(N_s[1]):
            ssr += s0*grid[i,j];
    return ssr/np.prod(N_s)

def susceptibility(grid, beta):
    '''
    formula chi = 1/T(<m^2> - <m>^2)
    :param grid:
    :param beta:
    :return:
    '''
    m = magnetization(grid);
    m_2 = magnetization_2(grid);
    chi = beta*(m_2 - m**2);
    return chi;

def heat_capacity(grid,K):
    E = energy(grid, K);
    E2 = energy_2(grid, K);
    cv = K**2*(E2 - E**2);
    return cv;
