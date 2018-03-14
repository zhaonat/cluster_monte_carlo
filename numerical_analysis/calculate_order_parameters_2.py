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
file = dir+'Ising_4D_Lattice_Temp_Scan.p';
file = dir+'50x50_Ising_2D_Lattice_Temp_Scan.p';

[simulation_data, epochs, beta_scan, N] = pickle.load(open(file, 'rb'));
print(str(len(simulation_data)) + ', ' + str(epochs))

thermalization = 0;
avg_data = list();

for key in simulation_data.keys():
    beta = key;
    print('beta: '+str(beta));

    lattice_history = simulation_data[key];
    c = 0; K = beta;
    N_s = np.prod(lattice_history[0].shape)
    m_avg = 0; m_2 = 0; E = 0; E_2 = 0;
    for Lattice in lattice_history:
        #calculate order parameters
        #magnetization
        m_avg += magnetization(Lattice);
        m_2 += magnetization(Lattice)**2;
        E += energy(Lattice, 1); #constant should not be beta....
        E_2 += energy(Lattice,1)**2;
        # chi = (m_2 - M*M)*beta_scan[beta];
        # HC = (E_2 - E*E)*beta_scan[beta]**2;
        #mag^2

        c = c+1;
    chi = beta*(m_2/c - (m_avg/c)**2);
    cv = beta**2*(E_2/c - (E/c)**2)
    avg_data.append([m_avg/c, E/c, chi, cv])
avg_data = np.array(avg_data);

# plt.plot(1/beta_scan, beta_data);
# plt.show();
vars = ['m per site', 'E per site', 'chi', 'heat capacity'];
for i in range(len(vars)):
    plt.figure();
    plt.plot(1/beta_scan, avg_data[:,i]);
    plt.title(vars[i])
    plt.xlabel('temperature')
    plt.ylabel(vars[i])

    # plt.figure();
    # plt.plot(beta_scan, beta_data[:,i]);
    # plt.title(vars[i])
    # plt.xlabel('temperature')
    # plt.ylabel(vars[i])

plt.show()

