import numpy as np

'''
    series of functions which will extract all major order parameters from the 
    Ising simulation
'''

def magnetization(grid):
    N_s = np.prod(grid.shape)
    return abs((1/N_s)*np.sum(grid));

def magnetization_2(grid):
    N_s = np.prod(grid.shape);
    return(1/N_s)*np.sum(grid*grid);

def energy(grid,J):
    N_s = np.prod(grid.shape);
    #displace the grid in four directions
    top = np.roll(grid, shift = 1, axis = [0,1])
    bottom = np.roll(grid, shift = 1, axis = [0,1])
    right = np.roll(grid, shift = 1, axis = [0,1])
    left = np.roll(grid, shift = 1, axis = [0,1])
    #factor of 0.5 is to account for double counting
    E = 0.5*np.sum(J*(grid*top + grid*bottom + grid*right + grid*left));
    return -E/N_s;

def energy_2(grid,J):
    ''' measurement of energy^2 per site, which should be constant'''
    N_s = np.prod(grid.shape);
    #displace the grid in four directions
    top = np.roll(grid, shift = 1, axis = [0,1])
    bottom = np.roll(grid, shift = 1, axis = [0,1])
    right = np.roll(grid, shift = 1, axis = [0,1])
    left = np.roll(grid, shift = 1, axis = [0,1])
    E  = 0.5*J**2*np.sum(np.square((grid*top) + np.square(grid*bottom) +\
                 np.square(grid*right) + np.square(grid*left)));
    return -E/N_s

def s0sr(grid):
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
    chi = beta*(m_2 - np.square(m));
    return chi;

def heat_capacity(grid,J, beta):
    E = energy_2(grid, J);
    E2 = energy_2(grid, J);
    cv = beta**2*(E2 - E**2);
    return cv;
