import numpy as np
from copy import copy
from core_functions.getNN import *

def select_neighbour(Lattice, visited, cluster, K, N, s_i):
    '''
    :param Lattice:
    :param visited: visited or not?
    :param cluster: list of elements belonging to cluster
    :param K: coupling constant strength (really beta*J)
    :param N: size of lattice
    :param s_i: current site or root
    :return:
    '''
    new_sites = []
    p = 1-np.exp(-2*K)
    num_NN = 1
    NN = getNN(s_i, N, num_NN)
    for NN_coords in NN:
        # accepting condition:
        NN_spin = Lattice[tuple(NN_coords)]
        if(NN_spin == Lattice[tuple(s_i)]):
            r = np.random.rand;
            if all([r < p, NN_spin not in cluster, NN_spin not in new_sites]):
                new_sites.append(NN_spin)   # mark it as new member of cluster
                cluster.append(NN_spin)     # add to cluster
                visited[NN_spin] = 1;
    return new_sites