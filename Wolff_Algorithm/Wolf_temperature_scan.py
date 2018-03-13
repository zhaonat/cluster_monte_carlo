import numpy as np
from Wolff_Algorithm.Wolff import  *
import matplotlib.pyplot as plt

''' function to verify the wolff algorithm gives good results on magnetization'''

K_scan = np.linspace(0.1, 0.9, 30)
N = (100,100);
Lattice  = 2*np.random.randint(0,2,N)-1;

#epochs = 1000;
# temp_data = list();
# for K in K_scan:
#     print('K: '+str(K))
#     Lattice, data = Wolff_simulation(Lattice,K, epochs, thermalization_epochs = 300);
#     m = np.mean(data);
#     temp_data.append(m);
#
# plt.figure()
# plt.plot(K_scan, temp_data);
# plt.show();

## smart wolff temperature scan
''' 
we'll run the simulation from multiple starting points and attempt to average it to resolve the behavior
'''
trials = 20;
epochs_0 = 300;
overall_mag = list();
for K in K_scan:
    if(K > 0.4):
        epochs = 100;
    else:
        epochs = epochs_0;
    temp_data = list();
    for i in range(trials):
        Lattice = 2 * np.random.randint(0, 2, N) - 1;
        print('K: '+str(K))
        Lattice, data = Wolff_simulation(Lattice,K, epochs, thermalization_epochs = 90);
        m = np.mean(data);
        temp_data.append(m);
    overall_mag.append(np.mean(temp_data));

plt.plot(K_scan, overall_mag);
plt.show();
