#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
#import matplotlib.cm as cm
#from matplotlib.colors import Normalize

import sys

#print "This is the name of the script: ", sys.argv[0]
#print "Number of arguments: ", len(sys.argv)
#print "The arguments are: " , str(sys.argv)

if(len(sys.argv) == 1) :
    n = str(0)
else:
    n = sys.argv[1]

plt.figure(figsize=(8,8))

path='./'

LL= 1

# Gallouet & Merrigot
#omega = 2 * np.pi / (31 * 0.01)

T_spring = 10 * 0.005
omega = 2 * np.pi / T_spring

print('omega = ' , omega)

Delta_t = 1  #optional prefactor

dt=np.loadtxt(path + n +'/particles.dat')

x=dt[:,0]; y=dt[:,1];
#    vol=dt[:,3]
w=dt[:,4];
#    vx=dt[:,5]; vym=dt[:,6];
#p=dt[:,9] / Delta_t**2
s=dt[:,10]
#  I=dt[:,11];

# Gallouet & Merrigot
#p = 0.5*omega**2 * w

r = np.sqrt( x**2 + y**2 )

plt.xlabel(r'$r$')
plt.ylabel(r'$s$')

plt.plot( r , s , 'x' )
   
#plt.xlim([-LL/2.0 , LL/2.0 ])
#plt.ylim([-LL/2.0 , LL/2.0 ])
    #    pl.colorbar(ticks=[0.45,0.55])

#print( 'step no ' + n )

plt.savefig( 'spring_cts_' + n + '.png')
plt.show()
