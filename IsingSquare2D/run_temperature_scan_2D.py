import numpy as np
from Wolff_Algorithm.Wolff import *
import matplotlib.pyplot as plt
from variable_calculations.order_measures import *
from metropolis_hastings.metropolis import *
from fitting_module.critical_exponents import fit_power_law
import pickle

N = (50,50);
lattice_size_name = str(N[0])+'x'+str(N[1]);

mag_v_temp =list();
beta_scan = np.linspace(0.3, 0.34, 20);
beta_near_crit = np.linspace(0.35, 0.54, 60);
beta_scan = np.append(beta_scan, beta_near_crit);
beta_scan = np.append(beta_scan, np.linspace(0.55, 1, 10));

simulation_data = dict(); #keys will be betas

## ===========LATTICE INITIALIZATION ============##
# we claim that it is better to do the initialization once and then run the simulation
# as we run to the next temperature, the previous lattice will actually be not too far from thesteady state lattice of the next
Lattice = 2 * np.random.randint(0, 2, N) - 1;

## ===============================================

for K in beta_scan:
    lattice_history = list();

    epochs = 1500;
    print(K)
    p = 1 - np.exp(-2 * K);
    magn = list();
    if(K> 0.4): #as temperature get higher, the lattice will equilibrate faster
        epochs = 500;
    if(K> 0.6):
        epochs = 200;

    ## simulation runner
    if(K > 0.5):
        print('metropolis')
        #run a metropolis hastings simulation
        for t in range(epochs):
            Lattice = metropolis_sim_epoch(Lattice, K, nearest_neighbors = 1);
            magn.append(magnetization(Lattice));
            if(t%100 == 0):
                print('epoch: '+str(t))
            lattice_history.append(Lattice);

    else:
        for t in range(epochs):
            Lattice = run_Wolff_epoch(Lattice, N, p);
            if(t%400 == 0):
                print(t);
                # plt.imshow(Lattice)
                # plt.show();
            if(t > 100):
                magn.append(magnetization(Lattice));
            lattice_history.append(Lattice);

    simulation_data[K] = lattice_history;
    M = np.mean(magn);
    mag_v_temp.append(M);

## ============= SAVE LATTICE HISTORY DATA =====================##
pickle.dump([simulation_data, epochs, beta_scan, N], open(lattice_size_name+'_Ising_2D_Lattice_Temp_Scan.p', 'wb'));
## ==============================================================

plt.figure();
plt.plot(1/beta_scan, mag_v_temp)
plt.show()

T_c = 1/2.269;
## fit power laws
critical_mag = fit_power_law(1/beta_scan, mag_v_temp, T_c);
print(critical_mag);



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
