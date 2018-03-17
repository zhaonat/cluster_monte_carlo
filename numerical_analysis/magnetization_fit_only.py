import numpy as np
import pickle
import settings
from variable_calculations.order_measures import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)

plt.close("all")

'''
    compact script to load and analyze MCMC data and produce order measures
'''
dir = settings.ROOT_DIR+'raw_data\\';
dir2 = settings.ROOT_DIR+'IsingCube\\';
#file = dir+'20x20x20_Ising_3D_Lattice_Temp_Scan.p';
file = dir+'Ising_4D_Lattice_Temp_Scan.p';
file = dir+'10x10X10_Ising_3D_Lattice_Temp_Scan.p';
file = dir+'20x20_Ising_2D_near_crit_Temp_Scan.p';
file = dir+'20x20_Ising_2D_near_crit_Temp_Scan.p';
#file = dir+'5x5x5x5x5_Ising_5D_Lattice_Temp_Scan.p';

[simulation_data, epochs, beta_scan, N, avg_data] = pickle.load(open(file, 'rb'));
print(str(len(simulation_data)) + ', ' + str(epochs))

## save data as dataframe, eventually move to the temp scan
temp_scan_data = pd.DataFrame(avg_data, columns = ['magnetization per site', 'E per site', 'chi', 'heat capacity'])
temp_scan_data['beta'] = beta_scan;
temp_scan_data.to_csv('2D_temp_data_another_run.csv')

plt.figure();
for i in range(4):
    plt.plot(1/beta_scan, avg_data[:,i]);
    plt.show()

thermalization = 60;
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
        #mag^2
        c = c+1;
    chi = beta*(m_2/c - (m_avg/c)**2);
    cv = beta**2*(E_2/c - (E/c)**2)
    avg_data.append([m_avg/c, E/c, chi, cv])
avg_data = np.array(avg_data);

# plt.plot(1/beta_scan, beta_data);
# plt.show();

## ==============  Plot parameters ========================== ##

SMALL_SIZE = 25
MEDIUM_SIZE = 25
BIGGER_SIZE = 25
plt.rc('axes', linewidth=2)
plt.rc('axes', )
plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=28)  # fontsize of the tick labels
plt.rc('ytick', labelsize=28)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

## =========================================================##
T_c = 4.451;

beta_crit = 1/T_c;
#get the beta_crit of the beta_scan
(bc_index, bc) = min(enumerate(beta_scan), key=lambda x: abs(x[1]-beta_crit))



vars = ['m per site', 'E per site', 'chi', 'heat capacity'];
for i in range(len(vars)):
    order_parameter = avg_data[:,i];
    plt.figure(figsize = (15,10));
    plt.plot(1/beta_scan, avg_data[:,i], linewidth=3);
    plt.title(vars[i])
    plt.xlabel('temperature (k_bT)')
    plt.ylabel(vars[i])


    # plt.figure();
    # plt.plot(beta_scan, beta_data[:,i]);
    # plt.title(vars[i])
    # plt.xlabel('temperature')
    # plt.ylabel(vars[i])

plt.show()

def fit_func(x, a, b, d):
    return a *(T_c- x) ** b + d
initial_guess = [1,0.1, 0];
magnetization = avg_data[bc_index:-15,0]
chi = avg_data[bc_index:-15,2]
beta_fit = beta_scan[bc_index:-15]


plt.plot(1/beta_fit, magnetization)
plt.plot(1/beta_fit, (1/beta_crit-1/beta_fit)**0.2)
y = fit_func(1/beta_fit, 1, 0.1, 0);
plt.plot(1/beta_fit,y)
params =curve_fit(fit_func, 1/beta_fit,magnetization, p0 = initial_guess,  maxfev=10000);
