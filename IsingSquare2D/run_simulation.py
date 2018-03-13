from cluster_functions.run_cluster import *
from variable_calculations import order_measures as OM
import numpy as np

def run_simulation_2D(Lattice, beta_scan, J, epochs, relax_time, NN):
    N = Lattice.shape;
    avg_data = list();
    for beta in beta_scan:
        print('beta: ' + str(beta))
        data = list();
        for t in range(epochs):
            Lattice = run_cluster_epoch(N, Lattice, beta, J, NN)
            if (t % 100 == 0):
                print('epoch: ' + str(t));
            # plt.imshow(Lattice)
            #     plt.show()
            if (t > relax_time):
                m = abs(OM.magnetization(Lattice))
                E = OM.energy(Lattice, J);
                sosr = abs(OM.s0sr(Lattice));
                chi = OM.susceptibility(Lattice, beta);
                cv = OM.heat_capacity(Lattice, J, beta)
                data.append([m, E, sosr, chi, cv])
        ## once the epoch if over, we need to average the quantities
        data = np.array(data);
        avg_data.append([np.mean(data, axis=0)])
    avg_data = np.squeeze(np.array(avg_data));
    return avg_data;

