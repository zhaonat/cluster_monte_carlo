import numpy as np
from Wolff_Algorithm.Wolff import *
import matplotlib.pyplot as plt
from variable_calculations.order_measures import *
from metropolis_hastings.metropolis import *
from fitting_module.critical_exponents import fit_power_law
import pickle

'''
illustrates that the metropolis algorithm is incomparably slow relative to the cluster algorithm
most particularly at high temperatures
'''

N = (40,40);
lattice_size_name = str(N[0])+'x'+str(N[1]);

mag_v_temp =list();
E_v_temp = list();
beta_scan = np.linspace(0.2, 0.34, 30);
beta_near_crit = np.linspace(0.35, 0.54, 60);
beta_scan = np.append(beta_scan, beta_near_crit);
beta_scan = np.append(beta_scan, np.linspace(0.55, 1, 20));

simulation_data = dict(); #keys will be betas

## ===========LATTICE INITIALIZATION ============##
# we claim that it is better to do the initialization once and then run the simulation
# as we run to the next temperature, the previous lattice will actually be not too far from thesteady state lattice of the next
Lattice = 2 * np.random.randint(0, 2, N) - 1;

## ===============================================
K_counter = 0;
for K in beta_scan:
    Lattice = 2 * np.random.randint(0, 2, N) - 1;

    lattice_history = list();
    epochs = 1500;
    print(K)
    p = 1 - np.exp(-2 * K);
    magn = list(); ene = list();
    if(K> 0.4): #as temperature get higher, the lattice will equilibrate faster
        epochs = 1000;
    if(K> 0.6):
        epochs = 1000;

    #run a metropolis hastings simulation
    for t in range(epochs):
        Lattice = metropolis_sim_epoch(Lattice, K, nearest_neighbors = 1);
        magn.append(magnetization(Lattice));
        if(t%100 == 0):
            print('epoch: '+str(t))
        lattice_history = list();

    # #run a wolff simulation
    # for t in range(epochs):
    #     Lattice = run_Wolff_epoch(Lattice, N, p);
    #     if(t%100 == 0):
    #         print(t);
    #     if(t > 100):
    #         magn.append(magnetization(Lattice));
    #         ene.append(energy(Lattice, 1)); #J = 1
    #     lattice_history.append(Lattice);

    simulation_data[K_counter] = lattice_history;
    M = np.mean(magn);
    E = np.mean(ene);
    mag_v_temp.append(M);
    E_v_temp.append(E);
    K_counter+=1;

## ============= SAVE LATTICE HISTORY DATA =====================##
pickle.dump([simulation_data, epochs, beta_scan, N], open(lattice_size_name+'_Ising_2D_random_Lattice_Temp_Scan_metropolis.p', 'wb'));
## ==============================================================

plt.figure();
plt.plot(1/beta_scan, mag_v_temp)

plt.figure();
plt.plot(1/beta_scan, E_v_temp)
plt.show()

T_c = 1/2.269;
# ## fit power laws
# critical_mag = fit_power_law(1/beta_scan, mag_v_temp, T_c);
# print(critical_mag);

# vars = ['m per site', 'E per site', 'spin_corr', 'chi', 'heat_capacity'];
# for i in range(len(vars)):
#     plt.figure();
#     plt.plot(1/beta_scan, avg_data[:,i]);
#     plt.title(vars[i])
#     plt.xlabel('temperature')
#     plt.ylabel(vars[i])
#     plt.savefig(vars[i]+'_vs_T.png');
# plt.show()
