import numpy as np
import matplotlib.pyplot as plt
from copy import copy;
from cluster_functions.graph_search import *
from variable_calculations import order_measures as OM
#try to write it in general so it is applicable for any dimension on cubic lattice
#also need periodic boundary condition
plt.close('all')

#specify initial 2D grid
Nx = 100;
Ny = 100;

#initialize grid of -1, and 1's, this is our starting lattice
import multiprocessing as mp

Lattice = 2*np.random.randint(0,2,(Nx,Ny))-1

#specify temperature, which we will do via beta
beta = 5;
J = 1; #up spin preferred

#simulation parameters
epoch = 400;
nearest_neighbors = 1;
#begin iterating over epochs
for t in range(epoch):
    print('epoch: '+str(t))
    #scan through every element of the lattice

    #propose a random lattice site to generate a cluster
    bonded = np.zeros((Nx,Ny));
    clusters = dict();  # keep track of bonds
    ## iterate through the entire lattice to assign bonds
    ## and clusters
    for i in range(Nx):
        for j in range(Ny):
            start_spin = Lattice[i,j];
            ## at this point, we do a BFS search to create the cluster
            bonded, clusters, visited = SW_BFS(Lattice, bonded, clusters, [i,j], beta, J);
    #print(bonded);
    if(t%40==0):
        plt.imshow(Lattice);
        plt.show()
    #print(clusters)
    #iterate through clusters
    #for each cluster, flip all the spins with probability 1/2

    for cluster_index in clusters.keys():
        [x0, y0] = np.unravel_index(cluster_index, (Nx,Ny));
        r = np.random.rand();
        if(r < 0.5):
            for coords in clusters[cluster_index]:
                [x,y] = coords;
                #print(Lattice[x,y], end=', '); #check clusters
                Lattice[x,y] = -1*Lattice[x,y];

    #calculate magnetization
    m = OM.magnetization(Lattice);
    print(m);

    plt.imshow(Lattice)
plt.show()