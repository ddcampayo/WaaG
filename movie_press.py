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
    init_t = 0
else:
    init_t = int( sys.argv[1] )

#import pylab as pl


plt.figure(figsize=(8,8))

skip=1
#path='timings_full/'
path='./'

LL= 1

omega = 2 * np.pi / (10 * 0.005)
print('omega = ' , omega)

def pp(r) : # analytic solution for Gresho's vortex pressure
    r0 = 0.2

    x = r/r0

    if( x < 1)  : return 5/2 - 4*np.log(2) + (x**2 -1)/2
    if( x < 2)  : return 4*np.log(x/2) - 4*(x-2) + (x**2 - 4)/2

    return 0;

v_pp = np.vectorize( pp )

plt.figure(figsize=(8,8))


for n in range( init_t ,2000000+skip,skip):

    strn = str(n)
    
    plt.clf()

    dt=np.loadtxt(path + strn +'/particles.dat')

    x=dt[:,0]; y=dt[:,1];
    #    vol=dt[:,3]
    w=dt[:,4];
    #    vx=dt[:,5]; vym=dt[:,6];
    p=dt[:,9]
    #  s=dt[:,10]
    #  I=dt[:,11];

    p = 0.5*omega**2 * w

    r = np.sqrt( x**2 + y**2 )

    #make furthest pressure value 0

    rm = np.argmax(r)

    p -= p[ rm ] #  np.min( p )

    plt.plot( r , p , 'o' )

    rr = np.linspace( 0 , max(r) , 200 )

    plt.plot( rr , v_pp(rr)  )

    print( 'pressure_{:03d}'.format( int(n/skip)  ) )

    plt.savefig( 'pressure_{:03d}'.format( int(n/skip)  ) )

#    plt.show()

    
