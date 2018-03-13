import numpy as np
import matplotlib.pyplot as plt
from copy import copy;
#try to write it in general so it is applicable for any dimension on cubic lattice
#also need periodic boundary condition
def getNN(site_indices, site_ranges, num_NN):
    Nearest_Neighbors = list();
    for i in range(len(site_indices)):
        for j in range(-num_NN,num_NN+1): #of nearest neighbors to include
            if(j == 0): continue;
            NN = copy(site_indices); #don't want to overwite;
            NN[i] +=j;
            if(NN[i] >= site_ranges[i]):
                NN[i] = NN[i] - site_ranges[i];
            if(NN[i] < 0):
                NN[i] = site_ranges[i]+NN[i];
            Nearest_Neighbors.append(NN)
    return Nearest_Neighbors;
#specify initial 2D grid
Nx = 200;
Ny = 200;

#initialize grid of -1, and 1's, this is our starting lattice
import multiprocessing as mp

Lattice = 2*np.random.randint(0,2,(Nx,Ny))-1

#specify temperature, which we will do via beta
beta = 0.5;
J = 1; #up spin preferred

#simulation parameters
epoch = 100;
nearest_neighbors = 1;
#begin iterating over epochs
for t in range(epoch):
    print('epoch: '+str(t))
    #scan through every element of the lattice
    for i in range(Nx):
        for j in range(Ny):
            site = Lattice[i,j];
            #propose a change
            proposal = -site;
            #calculate energy change from neighbors in this:
            NN = getNN([i,j], [Nx, Ny], nearest_neighbors);
            deltaE = 0;
            for a in NN:
                neighbor = Lattice[a[0], a[1]];
                deltaE += J*(neighbor*(site-proposal));

            if(deltaE<0): #accept immediately
                Lattice[i,j] = -1*Lattice[i,j];
            else:
                #calculate Boltzmann Weight
                Boltzmann = np.exp(-deltaE*beta);
                #generate random number
                p = np.random.rand();
                if(p < Boltzmann):
                    Lattice[i, j] = -1 * Lattice[i, j];
    if(t%20 == 0):
        plt.imshow(Lattice);
        plt.show()
