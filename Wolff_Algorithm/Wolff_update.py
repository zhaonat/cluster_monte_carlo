import numpy as np
from copy import copy
from core_functions.getNN import *
import matplotlib.pyplot as plt
from variable_calculations.order_measures import *
''' executes the SW algorithm by creating and operating on a single cluster at a time'''


def add_NN(Lattice, bonded, clusters, K, N, s_i):
    '''
    check and add nearest neighbors to a cluster
    :param Lattice:
    :param bonded: visited or not?
    :param clusters: list of elements belonging to cluster
    :param K: coupling constant strength (really beta*J)
    :param N: size of lattice
    :param s_i: current site or root
    :return:
    '''
    new_cluster_sites = []
    p = 1-np.exp(-2*K)
    num_NN = 1
    NN = getNN(s_i, N, num_NN);
    bonded[tuple(s_i)] = -1;
    for NN_coords in NN:
        # accepting condition:
        NN_spin = Lattice[tuple(NN_coords)]
        if(NN_spin == Lattice[tuple(s_i)]):
            r = np.random.rand();
            if( r < p and bonded[tuple(NN_coords)] == 1):
                new_cluster_sites.append(tuple(NN_coords))   # mark it as new member of cluster
                clusters.append(tuple(NN_coords))     # add to cluster
                bonded[tuple(NN_coords)] = -1;
    return new_cluster_sites, bonded

## we need an efficient way of generating the clusters
def generate_cluster(Lattice, visited, K, root, search_cutoff = 500):
    '''
    this is the graph search like function
    :param root: coordinate to generate cluster from
    :param Lattice: ising grid
    :param K: J*beta
    :return: coordinates of a single new cluster in as a list of tuples
    '''
    N = tuple(Lattice.shape);
    bonded = np.ones(N)
    #start = (0, 0)
    if(visited[root] == 0):
        queue = [root];
        cluster = [root];
        visited[root] = 1;c = 0;
        while (len(queue) > 0): ## cluster formation can become expensive when the clusters get large on large grid
            s_i = queue.pop(0);
            new_sites, bonded = add_NN(Lattice, bonded, cluster, K, N, s_i);
            for site in new_sites:
                queue.append(site)
            c+=1;
            if(c > search_cutoff):
                break;
        return cluster, bonded;
    else:
        return [], bonded;


## open question, how do we search the grid to search for viable clusters?
def run_cluster_sim(epochs, N, J, beta, disp_cutoff = 100):
    '''
    run cluster epoch simulations
    :param epochs:
    :param N:
    :param K:
    :return:
    '''
    relax_time =100;
    Lattice = 2 * np.random.randint(0, 2, N) - 1;
    bond_history = np.zeros(N);
    data = list();
    K = J*beta;
    for t in range(epochs):
        #select site to start cluster
        root = [];
        for i in range(len(N)):
            root.append(np.random.randint(0,N[i],1)[0])
        root = tuple(root);
        visited = np.zeros(N);
        cluster, bonded = generate_cluster(Lattice, visited, K, root)
        bond_history += bonded;
        Lattice = Lattice*bonded;
        if(t > relax_time):
            m = abs(magnetization(Lattice))
            E = energy(Lattice,J);
            sosr = abs(s0sr(Lattice));
            chi = susceptibility(Lattice, beta);
            cv = heat_capacity(Lattice,J,beta)
            data.append([m, E, sosr, chi, cv])

        if(t%disp_cutoff ==0):
            print('epoch: '+str(t))
            plt.subplot('121')
            plt.imshow(Lattice[:,:])
            plt.colorbar();
            plt.subplot('122')
            plt.imshow(bond_history[:,:]);
            bond_history = np.zeros(N);
            plt.colorbar()
            plt.draw();
    data = np.array(data)
    data = np.mean(data, axis = 0);
    return Lattice, data


# # basic tests of the cluster generator
# N = (5,5); N_s = np.prod(N);
# Lattice = 2*np.random.randint(0,2,N)-1;
# bonded = np.ones(N);
# # print(visited.shape);
# # print(Lattice.shape);
#
# #sites_to_visit;
# #this is a list of indices corresponding to sites that we can still visit
# # every time we visit a site, we pop its index from the list
# sites_to_visit = list(range(N_s));
# clusters = list();
# K = 2;
#
# # generate a site to visit
# start = (0,0)
# queue = [start];
# cluster = [start];
# counter = 0;
# while(len(queue)>0):
#     # if(counter> 0):
#     #     break;
#     s_i = queue.pop(0);
#     new_sites, bonded = add_NN(Lattice, bonded, cluster, K, N, s_i);
#     print(new_sites)
#     for site in new_sites:
#         queue.append(site)
#     counter+=1;
# print(cluster)
# cluster = np.array(cluster);
# plt.scatter(cluster[:,0], cluster[:,1]);
#
#
# print(bonded)
# print(Lattice)
# plt.figure();
# plt.imshow(np.rot90(Lattice))
# plt.show()
#

