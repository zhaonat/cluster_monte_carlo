import numpy as np

#in order to keep track of bonds, we create the dual lattice representing bonds
#and their nearest neighbors

Nx = 10;
Ny = 10;
Lattice = np.random.randint(0,2, (Nx,Ny));
print(Lattice)

bond_dual = np.zeros((Nx,Ny));
print(bond_dual)

#if I consider a lattice site indexed by (0,0), bonds are (1,0), (0,1), (Nx,0),(0,Ny)
# each site indexes a CENTER of the dual lattice

#dicsovered is an array marking whether a site is part of a cluster
discovered = np.zeros((Nx,Ny))