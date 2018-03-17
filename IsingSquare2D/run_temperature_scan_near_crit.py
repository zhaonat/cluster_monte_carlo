import numpy as np
from Wolff_Algorithm.Wolff import *
import matplotlib.pyplot as plt
from variable_calculations.order_measures import *
from metropolis_hastings.metropolis import *
from fitting_module.critical_exponents import fit_power_law
import pickle

N = (20,20);
lattice_size_name = str(N[0])+'x'+str(N[1]);

mag_v_temp =list();
E_v_temp = list();
beta_scan = np.linspace(0.35, 0.54, 50);

simulation_data = dict(); #keys will be betas

## ===========LATTICE INITIALIZATION ============##
# we claim that it is better to do the initialization once and then run the simulation
# as we run to the next temperature, the previous lattice will actually be not too far from thesteady state lattice of the next

##
np.random.seed(1); #lattice alwasy initializes to the same thing when we run again
Lattice = 2 * np.random.randint(0, 2, N) - 1;

## ===============================================
K_counter = 0; avg_data = list();
for K in beta_scan:
    #Lattice = 2 * np.random.randint(0, 2, N) - 1; --
    Lattice = np.ones(N);

    lattice_history = list();
    epochs = 20000;
    print(K)
    p = 1 - np.exp(-2 * K);
    magn = list(); ene = list();
    m_avg = 0;
    m_2 = 0;
    E = 0;
    E_2 = 0;
    c_avg = 0;
    #run a wolff simulation
    for t in range(epochs):
        Lattice = run_Wolff_epoch(Lattice, N, p);
        if(t%100 == 0):
            print(t);
        if(t > 1000):
            magn.append(magnetization(Lattice));
            ene.append(energy(Lattice, 1)); #J = 1
            m_avg += magnetization(Lattice);
            m_2 += magnetization(Lattice)**2;
            E += energy(Lattice, 1); #constant should not be beta....
            E_2 += energy(Lattice,1)**2;
            c_avg+=1;
        lattice_history.append(Lattice);
    chi = K*(m_2/c_avg - (m_avg/c_avg)**2);
    cv = K**2*(E_2/c_avg - (E/c_avg)**2)
    avg_data.append([m_avg/c_avg, E/c_avg, chi, cv])

    simulation_data[K_counter] = lattice_history;
    M = np.mean(magn);
    E = np.mean(ene);
    mag_v_temp.append(M);
    E_v_temp.append(E);
    K_counter+=1;

avg_data = np.array(avg_data);

## ============= SAVE LATTICE HISTORY DATA =====================##
pickle.dump([simulation_data, epochs, beta_scan, N, avg_data], open(lattice_size_name+'_Ising_2D_near_crit_Temp_Scan.p', 'wb'));
## ==============================================================



vars = ['Magnetization', 'Energy', 'Susceptibility', 'Heat Capacity'];
plt.figure(figsize = (20,15))
# fig, axes = plt.subplots(nrows=2, ncols=2)
# fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
d = str(2)
for i in range(len(vars)):
    plt.subplot(2,2,i+1);
    plt.plot(1/beta_scan[:], avg_data[:,i], linewidth=3);
    plt.title(vars[i]+ ' '+d+' Grid')
    plt.xlabel('temperature (k_bT)')
    plt.ylabel(vars[i]+' value')

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
