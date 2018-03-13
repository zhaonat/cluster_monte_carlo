import numpy as np
import pickle
import settings
from variable_calculations.order_measures import *
import matplotlib.pyplot as plt

'''
    compact script to load and analyze MCMC data and produce order measures
'''
dir = settings.ROOT_DIR+'\\numerical_analysis\\';

file = dir+'Ising_2D_Lattice_Temp_Scan.p';

[simulation_data, epochs, beta_scan, N] = pickle.load(open(file, 'rb'));
print(str(len(simulation_data)) + ', ' + str(epochs))

beta_data = list();
thermalization = 100;
for key in simulation_data.keys():
    lattice_history = simulation_data[key];
    avg_data = list();
    c = 0;
    for Lattice in lattice_history:
        #calculate order parameters
        #magnetization
        M = magnetization(Lattice);
        #mag^2
        if(c > thermalization):
            avg_data.append(M);
        c = c+1;
    beta_data.append(np.mean(avg_data));

plt.plot(1/beta_scan, beta_data);
plt.show();


