import numpy as np
from Wolff_Algorithm.Wolff import *
import matplotlib.pyplot as plt
from variable_calculations.order_measures import *
from metropolis_hastings.metropolis import *
import settings
N = (100,100);
Lattice  = 2*np.random.randint(0,2,N)-1;

figure_dir = settings.ROOT_DIR+'figure_dir\\'
K = 0.5;
epochs = 1000;
print(K)
p = 1 - np.exp(-2 * K);
plt.figure();
counter = 1
for t in range(epochs):
    Lattice = run_Wolff_epoch(Lattice, N, p);
    if(t%300 == 0):
        print(t);
        plt.subplot(2,2,(counter))
        plt.imshow(Lattice)
        title = 'Grid Visualization at epoch='+str(t)+' and K = '+str(K);
        plt.title(title)
        counter+=1;
plt.show();

plt.savefig(figure_dir+'\\'+'Ising_Visualization and K='+str(K)+'.png')
