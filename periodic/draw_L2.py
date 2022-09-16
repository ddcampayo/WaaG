#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np

LL= 1

dt=np.loadtxt( 'main.log', skiprows=1 )

T= 2*np.pi*0.2

t=dt[:,1] / T;

L2=dt[:,4]


#plt.tick_params(axis='y', which='minor')

plt.plot( t , L2)


plt.xlabel(r'$t/T$')
plt.ylabel(r'$L_2$')
plt.yscale('log')
#plt.ylim([-0.4 , 1 ])

plt.savefig( 'L2.png' )
