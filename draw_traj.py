#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import sys

import sizes

if(len(sys.argv) == 1) :
    infile = 'traj.dat'
else:
    infile = sys.argv[1]

LL= 1

dt=np.loadtxt( infile )
x=dt[:,1]; y=dt[:,2];
#vol=dt[:,3]
#w=dt[:,4];
#    vx=dt[:,5]; vym=dt[:,6];
p=dt[:,3]

p[0] = p[1] # initial pressure = 0 spoils the color range

plt.axis('scaled')

plt.scatter( x , y , 10, c=p )

plt.xlabel(r'$x$')
plt.ylabel(r'$y$')

plt.xlim([-LL/8.0 , LL/8.0 ])
plt.ylim([-LL/8.0 , LL/8.0 ])
#plt.clim(-1,1)
plt.colorbar()

outfile = infile + '.png'
plt.savefig( outfile , dpi=300, bbox_inches = "tight")
