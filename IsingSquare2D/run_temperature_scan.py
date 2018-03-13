from cluster_functions.run_cluster import *
from variable_calculations import order_measures as OM
#try to write it in general so it is applicable for any dimension on cubic lattice
#also need periodic boundary condition
plt.close('all')

#initialize grid of -1, and 1's, this is our starting lattice
#specify initial 2D grid
Nx = 50;
Ny = 50;
N = [Nx,Ny];
Lattice = 2*np.random.randint(0,2,(Nx,Ny))-1
NN = 1
#specify temperature, which we will do via beta

# initialize all parameters
epochs = 400;
mags = list();
J = 1;
avg_data = list();
relax_time = 100;

#notes: beta critical should be at beta ~0.55

beta_scan = np.linspace(0.3, 1.8,50)
for beta in beta_scan:
    print('beta: '+str(beta))
    data = list();
    for t in range(epochs):
        Lattice = run_cluster_epoch(N, Lattice, beta, J, NN)
        if(t%100 == 0):
            print('epoch: '+str(t));
        #     plt.imshow(Lattice)
        #     plt.show()
        if(t > relax_time):
            m = abs(OM.magnetization(Lattice))
            E = OM.energy(Lattice,J);
            sosr = abs(OM.s0sr(Lattice));
            chi = OM.susceptibility(Lattice, beta);
            cv = OM.heat_capacity(Lattice,J,beta)
            data.append([m, E, sosr, chi, cv])
    ## once the epoch if over, we need to average the quantities
    data = np.array(data);
    avg_data.append([np.mean(data, axis = 0)])
avg_data = np.squeeze(np.array(avg_data));

vars = ['m per site', 'E per site', 'spin_corr', 'chi', 'heat_capacity'];
for i in range(len(vars)):
    plt.figure();
    plt.plot(1/beta_scan, avg_data[:,i]);
    plt.title(vars[i])
    plt.xlabel('temperature')
    plt.ylabel(vars[i])
    plt.savefig(vars[i]+'_vs_T.png');

plt.show()
