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
figure_dir = settings.ROOT_DIR+'figure_dir\\'

#file = dir+'20x20x20_Ising_3D_Lattice_Temp_Scan.p';
file = dir+'Ising_4D_Lattice_Temp_Scan.p';
file = dir+'70x70_Ising_2D_Lattice_Temp_Scan.p';
#file = dir+'5x5x5x5x5_Ising_5D_Lattice_Temp_Scan.p';

[simulation_data, epochs, beta_scan, N] = pickle.load(open(file, 'rb'));
print(str(len(simulation_data)) + ', ' + str(epochs))

thermalization = 60;
avg_data = list();

for key in simulation_data.keys():
    beta = key;
    print('beta: '+str(beta));

    if(beta_scan[beta] > 1/2.4):
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
            #mag^2
            if(c%100 == 0):
                plt.figure();
                plt.imshow(Lattice);
                title = 'Ising Lattice at Epoch='+str(c)+'_beta='+str(np.round(beta_scan[beta],2))+'.png'
                print(title)
                plt.title(title);
                plt.savefig(figure_dir+'\\'+title)
                plt.xlabel('x direction')
                plt.ylabel('y direction')
                plt.show()

            c = c+1;



