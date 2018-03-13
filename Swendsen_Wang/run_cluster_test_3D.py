from Wolff_Algorithm.Wolff_update import *

#try to write it in general so it is applicable for any dimension on cubic lattice
#also need periodic boundary condition
plt.close('all')

#specify initial 2D grid
Nx = 50;
Ny = 50;
Nz = 50;
N = (Nx,Ny,Nz)
#initialize grid of -1, and 1's, this is our starting lattice

Lattice = 2*np.random.randint(0,2,(Nx,Ny,Nz))-1

#specify temperature, which we will do via beta
beta = 0.6;
J = 1; #up spin preferred
K = beta*J;
#simulation parameters
epochs = 1000;
nearest_neighbors = 1;

nearest_neighbors = 1;

Lattice = run_cluster_sim(epochs, N, K, disp_cutoff=200)

# #begin iterating over epochs
# for t in range(epoch):
#     print('epoch: '+str(t))
#     #scan through every element of the lattice
#
#     #propose a random lattice site to generate a cluster
#     bonded = np.zeros(N);
#     clusters = dict();  # keep track of bonds
#     ## iterate through the entire lattice to assign bonds
#     ## and clusters
#     for i in range(Nx):
#         for j in range(Ny):
#             for k in range(Nz):
#                 start_spin = Lattice[i,j,k];
#                 ## at this point, we do a BFS search to create the cluster
#                 bonded, clusters, visited = \
#                     SW_BFS(Lattice, bonded, clusters, [i,j,k], beta, J);
#     #print(bonded);
#     if(t%30==0):
#         plt.imshow(Lattice[:,:,2]);
#         plt.show()
#     #print(clusters)
#     #iterate through clusters
#     #for each cluster, flip all the spins with probability 1/2
#
#     for cluster_index in clusters.keys():
#         r = np.unravel_index(cluster_index, N);
#         p = np.random.rand();
#         if(p < 1):
#             for coords in clusters[cluster_index]:
#
#                 #print(Lattice[x,y], end=', '); #check clusters
#                 Lattice[tuple(coords)] = -1*Lattice[tuple(coords)];
#
#     #calculate magnetization
#     m = OM.magnetization(Lattice);
#     print(m);
#
#
# plt.show()