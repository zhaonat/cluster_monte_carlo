import numpy as np
from cluster_functions.getNN import *
import time
import matplotlib.pyplot as plt
'''
check how slow a DFS is in the case that the cluster is on the order of the entire grid
'''

N = (400,400)
visited = np.zeros(N);
Lattice = np.ones(N);
queue = list();
queue.append((0,0));
#visited[(0,0)] = 1;
t0 = time.time()
c = 0;
while(len(queue)>0):
    r = queue.pop(0);
    if(visited[r] == 0):
        visited[r] = 1;
        nearest_neighbors = getNN(r,N,1);
        for NN in nearest_neighbors:
            queue.append(NN);
    # if(c %10 == 0):
    #     plt.imshow(visited);
    #     plt.show()
    c+=1;

plt.imshow(visited);
plt.colorbar()
plt.show()
t1 = time.time();
print('time: '+str(t1-t0));