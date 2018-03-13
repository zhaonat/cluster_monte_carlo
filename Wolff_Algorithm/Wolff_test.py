import numpy as np
from Wolff_Algorithm.Wolff import *
import matplotlib.pyplot as plt
from variable_calculations.order_measures import *
from metropolis_hastings.metropolis import *
N = (50,50);
Lattice  = 2*np.random.randint(0,2,N)-1;


mag_v_temp =list();
for K in np.linspace(0.1, 0.6, 30):
    epochs = 1000;
    print(K)
    p = 1 - np.exp(-2 * K);
    magn = list();
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
plt.plot(np.linspace(0.1, 0.6, 30), mag_v_temp)
plt.show()

## 3D
# # 3D K_crit is 1/4.511
# N = (10,10,10);
# Lattice  = 2*np.random.randint(0,2,N)-1;
#
#
# mag_v_temp =list();
# beta_scan = np.linspace(0.05, 0.4, 30)
# for K in beta_scan:
#     epochs = 1000;
#     print(K)
#     p = 1 - np.exp(-2 * K);
#     magn = list();
#     if(K > 0.2217):
#         epochs = 100;
#
#     for t in range(epochs):
#         Lattice = run_Wolff_epoch(Lattice, N, p);
#         if(t%200 == 0):
#             print(t);
#             # plt.imshow(Lattice)
#             # plt.show();
#         if(t > 50):
#             magn.append(magnetization(Lattice));
#
#     M = np.mean(magn);
#     mag_v_temp.append(M);
#
# plt.figure();
# plt.plot(beta_scan, mag_v_temp)
# plt.show()
#
