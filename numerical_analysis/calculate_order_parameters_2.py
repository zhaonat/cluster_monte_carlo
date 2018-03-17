import numpy as np
import pickle
import settings
from variable_calculations.order_measures import *
import matplotlib.pyplot as plt
# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)

plt.close("all")
'''
    compact script to load and analyze MCMC data and produce order measures
    actually, this is bad if we save TOO MANY LATTICE HISTORIES
'''

dir = settings.ROOT_DIR+'raw_data\\';
#file = dir+'20x20x20_Ising_3D_Lattice_Temp_Scan.p';
file = dir+'Ising_4D_Lattice_Temp_Scan.p';
#file = dir+'100x100_Ising_2D_Lattice_Temp_Scan.p';
#file = dir+'5x5x5x5x5_Ising_5D_Lattice_Temp_Scan.p';
file = dir+'30x30x30_Ising_3D_random_Lattice_Temp_Scan.p';
file = dir+'Ising_5D_Lattice_Temp_Scan.p';

#file = dir+'7x7x7x7x7_Ising_5D_Lattice_Temp_Scan.p';

d = '3D'
[simulation_data, epochs, beta_scan, N] = pickle.load(open(file, 'rb'));
print(str(len(simulation_data)) + ', ' + str(epochs))

thermalization = 10;
avg_data = list();

for key in simulation_data.keys():
    beta = key;
    print('beta: '+str(beta));

    lattice_history = simulation_data[key];
    c = 0; K = beta;
    N_s = np.prod(lattice_history[0].shape)
    m_avg = 0; m_2 = 0; E = 0; E_2 = 0;
    c_avg = 0;
    for Lattice in lattice_history:
        #calculate order parameters
        #magnetization
        if(c > thermalization):
            m_avg += magnetization(Lattice);
            m_2 += magnetization(Lattice)**2;
            E += energy(Lattice, 1); #constant should not be beta....
            E_2 += energy(Lattice,1)**2;
            c_avg +=1;
            #mag^2
        c = c+1;
    chi = beta*(m_2/c_avg - (m_avg/c_avg)**2);
    cv = beta**2*(E_2/c_avg - (E/c_avg)**2)
    avg_data.append([m_avg/c_avg, E/c_avg, chi, cv])
avg_data = np.array(avg_data);

# plt.plot(1/beta_scan, beta_data);
# plt.show();

## ==============  Plot parameters ========================== ##

SMALL_SIZE = 20
MEDIUM_SIZE = 20
BIGGER_SIZE = 20
plt.rc('axes', linewidth=2)
plt.rc('axes', )
plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)  # fontsize of the tick labels
plt.rc('ytick', labelsize=18)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

## =========================================================##

vars = ['Magnetization', 'Energy', 'Susceptibility', 'Heat Capacity'];
plt.figure(figsize = (20,15))
# fig, axes = plt.subplots(nrows=2, ncols=2)
# fig.tight_layout() # Or equivalently,  "plt.tight_layout()"

for i in range(len(vars)):
    plt.subplot(2,2,i+1);
    plt.plot(1/beta_scan[0:], avg_data[:,i], linewidth=3);
    plt.title(vars[i]+ ' '+d+' Grid')
    plt.xlabel('temperature (k_bT)')
    plt.ylabel(vars[i]+' value')

    # plt.figure();
    # plt.plot(beta_scan, beta_data[:,i]);
    # plt.title(vars[i])
    # plt.xlabel('temperature')
    # plt.ylabel(vars[i])

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
plt.show()

