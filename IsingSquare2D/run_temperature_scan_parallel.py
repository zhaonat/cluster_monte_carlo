import numpy as np
from Wolff_Algorithm.Wolff import *
import matplotlib.pyplot as plt
from variable_calculations.order_measures import *
from metropolis_hastings.metropolis import *
N = (50,50);
Lattice  = 2*np.random.randint(0,2,N)-1;


mag_v_temp =list();
beta_scan = np.linspace(0.3, 0.34, 10);
beta_near_crit = np.linspace(0.35, 0.54, 50);
beta_scan = np.append(beta_scan, beta_near_crit);
beta_scan = np.append(beta_scan, np.linspace(0.55, 1, 10));

for K in beta_scan:
    epochs = 1000;
    print(K)
    p = 1 - np.exp(-2 * K);
    magn = list();
    if(K> 0.4): #as temperature get higher, the lattice will equilibrate faster
        epochs = 200;
    if(K> 0.6):
        epochs = 100;
    if(K > 0.5):
        print('metropolis')
        epochs = 200;
        #run a metropolis hastings simulation
        for t in range(epochs):
            Lattice = metropolis_sim_epoch(Lattice, K, nearest_neighbors = 1);
            magn.append(magnetization(Lattice));
            if(t%100 == 0):
                print('epoch: '+str(t))
    else:
        for t in range(epochs):
            Lattice = run_Wolff_epoch(Lattice, N, p);
            if(t%400 == 0):
                print(t);
                # plt.imshow(Lattice)
                # plt.show();
            if(t > 100):
                magn.append(magnetization(Lattice));

    M = np.mean(magn);
    mag_v_temp.append(M);

plt.figure();
plt.plot(1/beta_scan, mag_v_temp, '.-')
plt.xlabel('temperature')
plt.show()


# vars = ['m per site', 'E per site', 'spin_corr', 'chi', 'heat_capacity'];
# for i in range(len(vars)):
#     plt.figure();
#     plt.plot(1/beta_scan, avg_data[:,i]);
#     plt.title(vars[i])
#     plt.xlabel('temperature')
#     plt.ylabel(vars[i])
#     plt.savefig(vars[i]+'_vs_T.png');
#
# plt.show()
