from copy import copy

def getNN(site_indices, site_ranges, num_NN):
    '''
        site_indices: [i,j], site to get NN of
        site_ranges: [Nx,Ny], boundaries of the grid
        num_NN: number of nearest neighbors, usually 1

        function which gets NN on any d dimensional cubic grid
        with a periodic boundary condition
    '''

    Nearest_Neighbors = list();
    for i in range(len(site_indices)):
        for j in range(-num_NN,num_NN+1): #of nearest neighbors to include
            if(j == 0): continue;
            NN = list(copy(site_indices)); #don't want to overwite;
            NN[i] +=j;
            if(NN[i] >= site_ranges[i]):
                NN[i] = NN[i] - site_ranges[i];
            if(NN[i] < 0):
                NN[i] = site_ranges[i]+NN[i];
            Nearest_Neighbors.append(tuple(NN))
    return Nearest_Neighbors;