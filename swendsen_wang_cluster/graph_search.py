from core_functions.getNN import getNN
import numpy as np


# def SW_BFS(lattice, bonded, clusters, start, beta, J, nearest_neighbors = 1):
#     '''
#     :param lattice: lattice
#     :param bonded: 1 or 0, indicates whether a site has been assigned to a cluster
#            or not
#     :param clusters: dictionary containing all existing clusters, keys are an integer
#             denoting natural index of root of cluster
#     :param start: root node of graph (x,y)
#     :param beta: temperature
#     :param J: strength of lattice coupling
#     :param nearest_neighbors: number or NN to probe
#     :return:
#     '''
#     p = 1 - np.exp(-2 * beta * J); #bond forming probability
#     [Nx,Ny] = lattice.shape;
#     visited = np.zeros((Nx,Ny)); #indexes whether we have visited nodes during
#                                  #this particular BFS search
#     queue = list();
#     if(bonded[start[0], start[1]] != 0):
#         return bonded, clusters, visited;
#
#     queue.append(start);
#     cluster_spin = lattice[start[0], start[1]]
#     color = np.max(bonded) + 1;
#     index = sub2ind((Nx, Ny), start[1], start[0]); #cluster is built from
#     clusters[index] = list();
#     #whatever the input coordinates are
#     while(len(queue) > 0):
#         #print(queue)
#         [x,y] = queue.pop(0);
#         ##print(x,y)
#         if(visited[x, y] == 0): #if not visited
#             visited[x,y] = 1;
#             clusters[index].append([x, y]);
#             #to see clusters, always use different numbers
#             bonded[x, y] = color;
#             NN = getNN([x,y],[Nx,Ny], nearest_neighbors);
#             for nn_coords in NN:
#                 [xn,yn] = nn_coords;
#                 if(lattice[xn, yn] == cluster_spin and bonded[xn, yn] == 0\
#                    and visited[xn,yn] == 0): #require spins to be aligned
#                     random = np.random.rand();
#                     #print(str(p)+', '+str(random))
#                     if (random < p):  # accept bond proposal
#                         queue.append([xn,yn]); #add coordinate to search
#                         clusters[index].append([xn,yn])
#                         bonded[xn,yn] = color;
#     return bonded, clusters, visited;


## basic tests of cluster generation
# import matplotlib.pyplot as plt
#
# Nx = 10; Ny = 10;
# lattice = 2*np.random.randint(0,2,(Nx,Ny))-1
# lattice[0:5,0:5] = 1;
# #lattice = np.ones((Nx,Ny))
# #let's make the lattice clustered
# clusters = dict();
# plt.imshow(lattice);
#
# bonded = np.zeros((Nx, Ny));
#
# beta = 0.89; J = 1;
# p = 1-np.exp(-2*J*beta);
#
# # for i in range(Nx):
# #     for j in range(Ny):
# #         visited, clusters = SW_BFS(lattice, visited, clusters, [i,j], beta, J);
#
# i = 0; j = 0;
# bonded, clusters, visited = SW_BFS(lattice, bonded, clusters, [i, j], beta, J);
#
# print(len(clusters))
# print(bonded)
# xr = list(range(0,Nx));
# yr = list(range(0,Ny));
# [Xr, Yr] = np.meshgrid(xr,yr)
# plt.figure()
# plt.scatter(Xr.flatten(), Yr.flatten(), c = bonded.ravel(), cmap ='jet');
#
# plt.show()

## this should be dimension independent
def SW_BFS(lattice, bonded, clusters, start, beta, J, nearest_neighbors = 1):
    '''
    function currently cannot generalize to dimensions higher than 2...
    main idea is that we populate a lattice with clusters according to SW using a BFS from a root coord
    :param lattice: lattice
    :param bonded: 1 or 0, indicates whether a site has been assigned to a cluster
           or not
    :param clusters: dictionary containing all existing clusters, keys are an integer
            denoting natural index of root of cluster
    :param start: root node of graph (x,y)
    :param beta: temperature
    :param J: strength of lattice coupling
    :param nearest_neighbors: number or NN to probe
    :return:
    '''
    p = 1 - np.exp(-2 * beta * J); #bond forming probability
    N = lattice.shape;
    visited = np.zeros(N); #indexes whether we have visited nodes during
                                 #this particular BFS search
    queue = list();
    if(bonded[tuple(start)] != 0): #cannot construct a cluster from this site
        return bonded, clusters, visited;

    queue.append(start);
    cluster_spin = lattice[tuple(start)]
    color = np.max(bonded) + 1;

    ## need to make sub2ind work in arbitrary dimensions
    index = np.ravel_multi_index(tuple(start), dims=tuple(N), order='C')
    clusters[index] = list();
    #whatever the input coordinates are
    while(len(queue) > 0):
        #print(queue)
        r = tuple(queue.pop(0));
        ##print(x,y)
        if(visited[r] == 0): #if not visited
            visited[r] = 1;
            clusters[index].append(r);
            #to see clusters, always use different numbers
            bonded[r] = color;
            NN = getNN(r,N, nearest_neighbors);
            for nn_coords in NN:
                rn = tuple(nn_coords);
                if(lattice[rn] == cluster_spin and bonded[rn] == 0\
                   and visited[rn] == 0): #require spins to be aligned
                    random = np.random.rand();
                    if (random < p):  # accept bond proposal
                        queue.append(rn); #add coordinate to search
                        clusters[index].append(rn) #add point to the cluster
                        bonded[rn] = color; #indicate site is no longer available
    return bonded, clusters, visited;

## basic tests of cluster generation
#
# Nx = 10; Ny = 10; Nz = 10;
# N = (10,10,10)
# lattice = 2*np.random.randint(0,2,N)-1
# #lattice[0:5,0:5] = 1;
# #lattice = np.ones((Nx,Ny))
# #let's make the lattice clustered
# clusters = dict();
# plt.imshow(lattice[:,:,0]);
#
# bonded = np.zeros(N);
#
# beta = 2; J = 1;
# p = 1-np.exp(-2*J*beta);
#
# # for i in range(Nx):
# #     for j in range(Ny):
# #         visited, clusters = SW_BFS(lattice, visited, clusters, [i,j], beta, J);
#
# i = 0; j = 0;
# bonded, clusters, visited = SW_BFS(lattice, bonded, clusters, [i, j, 0], beta, J);
#
# print(len(clusters))
# print(bonded)
# xr = list(range(0,Nx));
# yr = list(range(0,Ny));
# [Xr, Yr] = np.meshgrid(xr,yr)
# plt.figure()
# plt.scatter(Xr.flatten(), Yr.flatten(), c = bonded[:,:,0].ravel(), cmap ='jet');
#
# plt.show()

# ## check 3D nearest neighbor search
# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# nn = np.array(getNN((1,1,1),(10,10,10), 1));
# ax.scatter(nn[:,0], nn[:,1], nn[:,2])
#
