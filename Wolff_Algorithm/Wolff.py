import numpy as np
from core_functions.getNN import getNN;
import matplotlib.pyplot as plt
from variable_calculations.order_measures import *

## Wolf function
def Wolff_simulation(Lattice,K, epochs, thermalization_epochs = 100, num_views = 10):
    '''
    :param Lattice:
    :param K:
    :param epochs:
    :param thermalization_epochs:
    :param num_views:
    :return:
    '''
    plt.ion();
    ## wolff test using Frontier idea
    N = Lattice.shape;
    # generate random particle
    Lattice_History = list();
    p = 1 - np.exp(-2 * K)
    data = list();
    for t in range(epochs):
        change_tracker = np.ones(N);
        visited = np.zeros(N);
        root = []; # generate random coordinate by sampling from uniform random...
        for i in range(len(N)):
            root.append(np.random.randint(0, N[i], 1)[0])
        root = tuple(root);
        visited[root]=1;
        C = [root];  # denotes cluster coordinates
        F_old = [root];  # old frontier
        change_tracker[root] = -1;
        while (len(F_old) != 0):
            F_new = [];
            for site in F_old:
                site_spin = Lattice[tuple(site)]
                # get neighbors
                NN_list = getNN(site, N, num_NN=1);
                for NN_site in NN_list: ## if we do the full search, this is bad, because
                    nn = tuple(NN_site)
                    if (Lattice[nn] == site_spin and visited[nn] == 0):
                        if (np.random.rand() < p):
                            F_new.append(nn); visited[nn] = 1;
                            C.append(nn);
                            change_tracker[nn] = -1;
            F_old = F_new;

        # update the cluster
        Lattice = Lattice*change_tracker;
        if(t > thermalization_epochs):
            data.append(magnetization(Lattice));
        # for site in C:
        #     Lattice[site] = -1 * Lattice[site]
        # if (t % int(epochs/num_views) == 0):
        #     print('epoch: ' + str(t));
        #     plt.imshow(Lattice);
        #     plt.pause(0.05)

    # plt.imshow(Lattice);
    # plt.show()
    return Lattice, data;

def run_Wolff_epoch(Lattice, N, p):
    '''
    run one wolff epoch
    :param Lattice:
    :param N:
    :param p: 1-exp(-2K);
    :return:
    '''

    change_tracker = np.ones(N);
    visited = np.zeros(N);
    root = [];  # generate random coordinate by sampling from uniform random...
    for i in range(len(N)):
        root.append(np.random.randint(0, N[i], 1)[0])
    root = tuple(root);
    visited[root] = 1;
    C = [root];  # denotes cluster coordinates
    F_old = [root];  # old frontier
    change_tracker[root] = -1;
    while (len(F_old) != 0):
        F_new = [];
        for site in F_old:
            site_spin = Lattice[tuple(site)]
            # get neighbors
            NN_list = getNN(site, N, num_NN=1);
            for NN_site in NN_list:  ## if we do the full search, this is bad, because
                nn = tuple(NN_site)
                if (Lattice[nn] == site_spin and visited[nn] == 0):
                    if (np.random.rand() < p):
                        F_new.append(nn);
                        visited[nn] = 1;
                        C.append(nn);
                        change_tracker[nn] = -1;
        F_old = F_new;
    Lattice = Lattice * change_tracker;

    return Lattice;
# ## wolff test using Frontier idea
# N = (200,200);
# #lattice
# Lattice = 2*np.random.randint(0,2,N)-1;
# #Lattice = 2*np.ones(N);
# #generate random particle
# change_tracker = np.ones(N);
#
# K = 0.8;
# p = 1-np.exp(-2*K)
# epochs = 800; #epochs scales with number of sites that must be probed...
# for t in range(epochs):
#     root = [];
#     for i in range(len(N)):
#         root.append(np.random.randint(0,N[i],1)[0])
#     root = tuple(root);
#     #print('root: '+str(root))
#     C = [root]; # denotes cluster coordinates
#     visited = np.zeros(N);
#     visited[root] = 1;
#     F_old = [root]; #old frontier
#     change_tracker[root] = -1;
#     computations_tracker = 0;
#     while(len(F_old) != 0): #there are situations where the cluster search slows down when the clusters get large
#         F_new = [];
#         for site in F_old:
#             site_spin = Lattice[tuple(site)]
#             #get neighbors
#             NN_list = getNN(site, N, num_NN= 1);
#             for NN_site in NN_list: ## this has to probe through 2*d nearest neighbors...
#                 computations_tracker+=1;
#                 nn = tuple(NN_site)
#                 if(Lattice[nn] == site_spin and visited[nn] == 0): # original codenn not in C):
#                     ## nn not in C is an expensive computation if C is just a simple array
#                     if(np.random.rand() < p):
#                         F_new.append(nn);
#                         C.append(nn); visited[nn] = 1;
#                         change_tracker[nn] = -1;
#         F_old = F_new;
#
#     #at the end we have a cluster
#     #print(C)
#     #print(len(C))
#     print('computations: '+str(computations_tracker))
#     #update the cluster
#     Lattice = Lattice*change_tracker;
#     # for site in C:
#     #     Lattice[site] = -1*Lattice[site]
#     if(t%80 == 0):
#         plt.imshow(Lattice);
#         plt.show();
#
# plt.imshow(Lattice);
# plt.show()
