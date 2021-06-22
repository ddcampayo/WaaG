#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np

LL= 1

dt=np.loadtxt( 'traj.dat' )
x=dt[:,1]; y=dt[:,2];
#vol=dt[:,3]
#w=dt[:,4];
#    vx=dt[:,5]; vym=dt[:,6];
p=dt[:,3]

plt.axis('scaled')

plt.scatter( x , y , 10, c=p )

plt.xlabel(r'$x$')
plt.ylabel(r'$y$')

plt.xlim([-LL/2.0 , LL/2.0 ])
plt.ylim([-LL/2.0 , LL/2.0 ])
plt.clim(-2,2)
plt.colorbar()

plt.savefig( 'traj.png' )
