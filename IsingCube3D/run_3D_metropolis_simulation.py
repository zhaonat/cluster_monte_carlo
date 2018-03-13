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
Nx = 20;
Ny = 100;
Nz = 100;
#initialize grid of -1, and 1's, this is our starting lattice
import multiprocessing as mp

Lattice = 2*np.random.randint(0,2,(Nx,Ny,Nz))-1

#specify temperature, which we will do via beta
beta = 0.4;
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
            for k in range(Nz):
                site = Lattice[i,j, k];
                #propose a change
                proposal = -site;
                #calculate energy change from neighbors in this:
                NN = getNN([i,j, k], [Nx, Ny, Nz], nearest_neighbors);
                deltaE = 0;
                for a in NN:
                    neighbor = Lattice[a[0], a[1], a[2]];
                    deltaE += J*(neighbor*(site-proposal));

                if(deltaE<0): #accept immediately
                    Lattice[i,j,k] = -1*Lattice[i,j,k];
                else:
                    #calculate Boltzmann Weight
                    Boltzmann = np.exp(-deltaE*beta);
                    #generate random number
                    p = np.random.rand();
                    if(p < Boltzmann):
                        Lattice[i, j, k] = -1 * Lattice[i, j, k];
    if(t%10 == 0):
        plt.imshow(Lattice[5,:,:]);
        plt.show()
