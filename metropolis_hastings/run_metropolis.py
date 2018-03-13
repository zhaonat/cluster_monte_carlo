import numpy as np
import matplotlib.pyplot as plt
from copy import copy;
from metropolis_hastings.metropolis import  *
#try to write it in general so it is applicable for any dimension on cubic lattice
#specify initial 2D grid
Nx = 50;
Ny = 50;
N = (Nx,Ny)
#initialize grid of -1, and 1's, this is our starting lattice
import multiprocessing as mp

Lattice = 2*np.random.randint(0,2,N)-1

#specify temperature, which we will do via beta
beta = 2;
J = 1; #up spin preferred

K = J*beta;
epochs = 200;
metropolis_simulation(Lattice, K, epochs, num_views = 1);

## 3D
# Nx = 40;
# Ny = 40;
# Nz = 40
# N = (Nx,Ny,Nz)
# #initialize grid of -1, and 1's, this is our starting lattice
# import multiprocessing as mp
#
# Lattice = 2*np.random.randint(0,2,N)-1
#
# #specify temperature, which we will do via beta
# beta = 2;
# J = 1; #up spin preferred
#
# K = J*beta;
# epochs = 100;
# metropolis_simulation(Lattice, K, epochs, num_views = 20);