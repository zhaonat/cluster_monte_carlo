import numpy as np
import matplotlib.pyplot as plt
from core_functions.getNN import getNN
from variable_calculations.order_measures import *

def metropolis_simulation(Lattice, K, epochs, thermalization = 100, num_views = 1, view_lattice = False):
    N = Lattice.shape;
    N_sites = np.prod(N);
    nearest_neighbors = 1;
    # begin iterating over epochs
    data = list();
    for t in range(epochs):
        computations_tracker = 0;
        # scan through every element of the lattice
        for i in range(N_sites):
            #convert i to a coordinate
            r = np.unravel_index(i,N);
            computations_tracker += 1;
            site = Lattice[r];
            # propose a change
            proposal = -site;
            # calculate energy change from neighbors in this:
            NN = getNN(r, N, nearest_neighbors);
            deltaE = 0;
            for nn_site in NN:
                neighbor = Lattice[nn_site];
                deltaE += K * (neighbor * (site - proposal));

            if (deltaE < 0):  # accept immediately
                Lattice[r] = -1 * Lattice[r];
            else:
                # calculate Boltzmann Weight
                Boltzmann = np.exp(-deltaE);
                # generate random number
                p = np.random.rand();
                if (p < Boltzmann):
                    Lattice[r] = -1 * Lattice[r];
        if (t % int(epochs/num_views) == 0):
            print('epoch: '+str(t))
        if (t % int(epochs/num_views) == 0 and view_lattice):
            plt.imshow(Lattice);
            plt.show()
        if(t > thermalization):
            data.append(magnetization(Lattice))
    return Lattice, data

def metropolis_sim_epoch(Lattice, K, nearest_neighbors = 1):
    '''
    :param Lattice:
    :param K: beta*J, temperature/coupling strength
    :param nearest_neighbors: number of nearest neighbors
    :return:
    '''
    N = Lattice.shape;
    N_sites = np.prod(N);
    for i in range(N_sites):
        #convert i to a coordinate
        r = np.unravel_index(i,N);
        site = Lattice[r];
        # propose a change
        proposal = -site;
        # calculate energy change from neighbors in this:
        NN = getNN(r, N, nearest_neighbors);
        deltaE = 0;
        for nn_site in NN:
            neighbor = Lattice[nn_site];
            deltaE += K * (neighbor * (site - proposal));
        if (deltaE < 0):  # accept immediately
            Lattice[r] = -1 * Lattice[r];
        else:
            # calculate Boltzmann Weight
            Boltzmann = np.exp(-deltaE);
            # generate random number
            p = np.random.rand();
            if (p < Boltzmann):
                Lattice[r] = -1 * Lattice[r];
    return Lattice

''' work in progress'''
def metropolis_sim_vectorized(Lattice, K, epochs, num_views = 1):
    N = Lattice.shape;
    N_sites = np.prod(N);
    epoch = 100;
    nearest_neighbors = 1;
    # begin iterating over epochs
    for t in range(epoch):
        print('epoch: ' + str(t))
        computations_tracker = 0;
        # scan through every element of the lattice

        #for every site, propose a change as a matrix
        proposal_matrix = -Lattice;

        #calculate energy change using the original Lattice
        # need a shift operator...

        for i in range(N_sites):
            # convert i to a coordinate
            r = np.unravel_index(i, N);
            site = Lattice[r];
            # propose a change
            proposal = -site;


            # calculate energy change from neighbors in this:
            NN = getNN(r, N, nearest_neighbors);
            deltaE = 0;
            for a in NN:
                neighbor = Lattice[a[0], a[1]];
                deltaE += K * (neighbor * (site - proposal));

            if (deltaE < 0):  # accept immediately...we initiate a change with every flip so that NN may be affected
                              ## in a vectorized version, we would actually have to do this on a checkerboard...
                Lattice[r] = -1 * Lattice[r];
            else:
                # calculate Boltzmann Weight
                Boltzmann = np.exp(-deltaE);
                # generate random number
                p = np.random.rand();
                if (p < Boltzmann):
                    Lattice[r] = -1 * Lattice[r];
        if (t % int(epochs / num_views) == 0):
            plt.imshow(Lattice);
            plt.show()



