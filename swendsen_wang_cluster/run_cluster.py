from swendsen_wang_cluster.graph_search import *


#try to write it in general so it is applicable for any dimension on cubic lattice
#also need periodic boundary condition


def run_cluster_epoch(N, Lattice, beta, J, NN):
    Nx = N[0]; Ny = N[1];
    #simulation parameters
    nearest_neighbors = NN;
    #scan through every element of the lattice

    #propose a random lattice site to generate a cluster
    bonded = np.zeros((Nx,Ny));
    clusters = dict();  # keep track of bonds
    ## iterate through the entire lattice to assign bonds
    ## and clusters
    for i in range(Nx):
        for j in range(Ny):
            ## at this point, we do a BFS search to create the cluster
            bonded, clusters, visited = SW_BFS(Lattice, bonded, clusters, [i,j], beta, J,nearest_neighbors=NN);

    for cluster_index in clusters.keys():
        [x0, y0] = np.unravel_index(cluster_index, (Nx,Ny));
        r = np.random.rand();
        if(r < 0.5):
            for coords in clusters[cluster_index]:
                [x,y] = coords;
                #print(Lattice[x,y], end=', '); #check clusters
                Lattice[x,y] = -1*Lattice[x,y];

    return Lattice;


## function which runs cluster_epoch in D dimensions
def run_cluster_epoch(N, Lattice, beta, J, NN):
    #simulation parameters
    nearest_neighbors = NN;
    N_s = np.prod(N);
    #scan through every element of the lattice

    #propose a random lattice site to generate a cluster
    bonded = np.zeros(N);
    clusters = dict();  # keep track of bonds
    ## iterate through the entire lattice to assign bonds
    ## and clusters
    for i in range(N_s):
        r = np.unravel_index(i,N)
        #convert natural index into coordinates

        ## at this point, we do a BFS search to create the cluster
        bonded, clusters, visited = SW_BFS(Lattice, bonded, clusters, r, beta, J,nearest_neighbors=NN);

    for cluster_index in clusters.keys():
        coord_0 = np.unravel_index(cluster_index, N);
        r = np.random.rand();
        if(r < 0.5):
            for coords in clusters[cluster_index]:
                #print(Lattice[x,y], end=', '); #check clusters
                Lattice[coords] = -1*Lattice[coords];

    return Lattice;