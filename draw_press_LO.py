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

LL= 10

Re = 1000

t = float(n)

dt=np.loadtxt(path + n +'/particles.dat')

x=dt[:,0]; y=dt[:,1];
#    vol=dt[:,3]
#w=dt[:,4];
#vx=dt[:,5]; vy=dt[:,6];
p=dt[:,9]
#  s=dt[:,10]
#  I=dt[:,11];

#p = 0.5*omega**2 * w

r = np.sqrt( x**2 + y**2 )
#v = np.sqrt( vx**2 + vy**2 )

#make furthest pressure value 0

rm = np.argmax(r)

# p -= p[ rm ] #  np.min( p )

#plt.plot( r , v , 'o' )

def ff(t) :
    return 100.0 # Re / ( 4 * np.pi * t + Re )

def gg(r,t) :
    return ff(t) * r**2 # Re / ( 4 * np.pi * t + Re )

def pp(r,t) : # analytic solution for Lamb Oseen vortex velocity

    if(r < 1e-10) :
        return 0.0
    f = ff(t)
    g = gg(r,t)
    
    return  \
            f / 8 * (  \
            2 * np.log(2) + (2*np.exp(-g ) - np.exp(-2*g ) -1)/g ) - \
            f  * np.log(2)  /4




v_pp = np.vectorize( pp )
rr = np.linspace( 0 , max(r) , 200 )
plt.plot( rr , v_pp(rr,t)  )

plt.xlabel(r'$r$')
plt.ylabel(r'$p/\rho$')


#plt.plot( r , w , 'x' )
   
#plt.xlim([-LL/2.0 , LL/2.0 ])
#plt.ylim([-LL/2.0 , LL/2.0 ])
    #    pl.colorbar(ticks=[0.45,0.55])

#print( 'step no ' + n )

plt.savefig( 'pressure_LO_' + n + '.png')
plt.show()

