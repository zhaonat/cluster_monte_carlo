import numpy as np
from Wolff_Algorithm.Wolff import *
import matplotlib.pyplot as plt
from variable_calculations.order_measures import *
from metropolis_hastings.metropolis import *
from fitting_module.critical_exponents import fit_power_law
import pickle
import pandas as pd

N = (5,5,5,5,5);
d = str(len(N));
lattice_size_name = '';
for i in range(len(N)):
    lattice_size_name +=str(N[i])+'x';
    if(i == len(N)-1):
        lattice_size_name += str(N[i]);

print('simulation at dimension: '+str(d));
mag_v_temp =list();
E_v_temp = list();

## beta_scans for different dimensions
beta_scan = np.linspace(0.35, 0.54, 50); #2D
#beta_scan = np.linspace(0.2, 0.4, 40); #3D
beta_scan = np.linspace(0.12, 0.25, 50); #4D
# beta_scan = np.linspace(0.1, 0.2, 50);

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
    epochs = 5000;
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
        if(t%1000 == 0):
            print(str(K)+', '+str(t));
        magn.append(magnetization(Lattice));
        ene.append(energy(Lattice, 1));  # J = 1
        m_avg += magnetization(Lattice);
        m_2 += magnetization(Lattice) ** 2;
        E += energy(Lattice, 1);  # constant should not be beta....
        E_2 += energy(Lattice, 1) ** 2;
        c_avg += 1;
        if(t > epochs - 100): #only save the last 100 epochs...
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
pickle.dump([simulation_data, epochs, beta_scan, N, avg_data], open(lattice_size_name+'_Ising_'+d+'D_near_crit_Temp_Scan.p', 'wb'));
## ==============================================================

## save data as dataframe, eventually move to the temp scan
temp_scan_data = pd.DataFrame(avg_data, columns = ['magnetization per site', 'E per site', 'chi', 'heat capacity'])
temp_scan_data['beta'] = beta_scan;
temp_scan_data.to_csv(d+'D_temp_data.csv')

vars = ['Magnetization', 'Energy', 'Susceptibility', 'Heat Capacity'];
plt.figure(figsize = (20,15))
# fig, axes = plt.subplots(nrows=2, ncols=2)
# fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
for i in range(len(vars)):
    plt.subplot(2,2,i+1);
    plt.plot(1/beta_scan[:], avg_data[:,i], linewidth=3);
    plt.title(vars[i]+ ' '+d+' Grid')
    plt.xlabel('temperature (k_bT)')
    plt.ylabel(vars[i]+' value')


