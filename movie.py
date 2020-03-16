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

for n in range( init_t ,2000000+skip,skip):
    plt.clf()
    dt=np.loadtxt(path+str(n)+'/particles.dat')
 
    x=dt[:,0]; y=dt[:,1];
#    vol=dt[:,3]
#    w=dt[:,4];
#    vx=dt[:,5]; vym=dt[:,6];
    p=dt[:,9]
#    p=dt[:,10]

    w=dt[:,11];

    p += p - np.min( p )

    #    r = np.sqrt( x**2 + y**2 )
    #    plt.plot( r , p , 'o' )

#    plt.scatter( x , y , c=p )
    plt.scatter( x , y , c=w )

    di = np.loadtxt(path+str(n)+'/diagram.dat')

    xd=di[:,0]; yd=di[:,1];

    for i in range( 0 , xd.size , 2) :
        plt.plot( [ xd[i] , xd[i+1] ] , [ yd[i] , yd[i+1] ] , c='k')

    
    plt.xlim([-LL/2.0 , LL/2.0 ])
    plt.ylim([-LL/2.0 , LL/2.0 ])
    #    pl.colorbar(ticks=[0.45,0.55])

    print( 'snap{:03d}'.format( int(n/skip)  ) )

    plt.savefig( 'snap{:03d}'.format( int(n/skip)  ) )

    
