#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
#import matplotlib.cm as cm
#from matplotlib.colors import Normalize

import sys

if(len(sys.argv) == 1) :
    n = str(0)
else:
    n = sys.argv[1]

plt.figure(figsize=(8,8))

path='./'

LL= 1

dt=np.loadtxt(path + n +'/particles.dat')

x=dt[:,0]; y=dt[:,1];
#    vol=dt[:,3]
w=dt[:,4];
#    vx=dt[:,5]; vym=dt[:,6];
gradp_x =dt[:,12]
gradp_y =dt[:,13]
#  s=dt[:,10]
#  I=dt[:,11];

# Gallouet & Merrigot
#p = 0.5*omega**2 * w

r     = np.sqrt(       x**2 +        y**2 )

# radial:
gradp = ( gradp_x * x +  gradp_y * y ) / r

plt.plot( r , gradp , 'o' )

def gg(r) : # analytic solution for Gresho's vortex pressure gradient
    r0 = 0.2

    x = r/r0

    if( x < 1)  : return x/r0
    if( x < 2)  : return (2-x)**2/r


    return 0;




v_gg = np.vectorize( gg )
rr = np.linspace( 0 , max(r) , 200 )
plt.plot( rr , v_gg(rr)  )

plt.xlabel(r'$r$')
plt.ylabel(r'$\nabla p/\rho$')


#plt.plot( r , w , 'x' )
   
#plt.xlim([-LL/2.0 , LL/2.0 ])
#plt.ylim([-LL/2.0 , LL/2.0 ])
    #    pl.colorbar(ticks=[0.45,0.55])

#print( 'step no ' + n )

plt.savefig( 'grad_pressure_' + n + '.png')
plt.show()

