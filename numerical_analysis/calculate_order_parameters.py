import numpy as np
import pickle
import settings
from variable_calculations.order_measures import *
import matplotlib.pyplot as plt
plt.close("all")

'''
    compact script to load and analyze MCMC data and produce order measures
'''
dir = settings.ROOT_DIR+'numerical_analysis\\';

#file = dir+'20x20x20_Ising_3D_Lattice_Temp_Scan.p';
#file = dir+'Ising_4D_Lattice_Temp_Scan.p';
file = dir+'30x30_Ising_2D_random_Lattice_Temp_Scan.p';

[simulation_data, epochs, beta_scan, N] = pickle.load(open(file, 'rb'));
print(str(len(simulation_data)) + ', ' + str(epochs))

beta_data = list();
thermalization = 0;
for key in simulation_data.keys():
    beta = key;
    print('beta: '+str(beta_scan[beta]));

    lattice_history = simulation_data[key];
    avg_data = list();
    c = 0; K = beta;
    N_s = np.prod(lattice_history[0].shape)
    n_epochs = len(lattice_history)-thermalization;
    n1, n2 = 1.0 / (n_epochs), 1.0 / (n_epochs **2)
    for Lattice in lattice_history:
        #calculate order parameters
        #magnetization
        M = magnetization(Lattice);
        m_2 = magnetization_2(Lattice);
        E = energy(Lattice, 1); #constant should not be beta....
        E_2 = energy_2(Lattice,1);
        chi = (m_2 - M*M)*beta_scan[beta];
        HC = (E_2 - E*E)*beta_scan[beta]**2;
        #mag^2
        if(c > thermalization):
            avg_data.append([M, E, chi, HC]);
        c = c+1;
    avg_data = np.array(avg_data);
    beta_data.append(np.mean(avg_data, axis = 0));

# plt.plot(1/beta_scan, beta_data);
# plt.show();
beta_data = np.array(beta_data);
vars = ['m per site', 'E per site', 'chi', 'heat capacity'];
for i in range(len(vars)):
    plt.figure();
    plt.plot(1/beta_scan, beta_data[:,i]);
    plt.title(vars[i])
    plt.xlabel('temperature')
    plt.ylabel(vars[i])

    # plt.figure();
    # plt.plot(beta_scan, beta_data[:,i]);
    # plt.title(vars[i])
    # plt.xlabel('temperature')
    # plt.ylabel(vars[i])

plt.show()

