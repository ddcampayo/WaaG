#!/usr/bin/python

import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('GTK3Agg')
matplotlib.use('Cairo')
import numpy as np
#import matplotlib.cm as cm
#from matplotlib.colors import Normalize

import sys
import glob


#print "This is the name of the script: ", sys.argv[0]
#print "Number of arguments: ", len(sys.argv)
#print "The arguments are: " , str(sys.argv)

if(len(sys.argv) == 1) :
    all_times = True
    init_t = 0
else:
    all_times = False
    init_t = sys.argv[1]

#import pylab as pl


plt.figure(figsize=(8,8))

skip=1
#path='timings_full/'
path='./'

LL= 1


T_spring = 10 * 0.005
omega = 2 * np.pi / T_spring


def by_number(elem):
    return float(elem)

if all_times:
    times = sorted(glob.glob('[0-9]*'), key=by_number)
else:
    times = [init_t]

for time in times:

    part_file = time+'/particles.dat'
    diag_file = time+'/diagram.dat'

#for n in range( init_t ,2000000+skip,skip):
    plt.clf()

    #dt=np.loadtxt(path+str(n)+'/particles.dat')
    dt=np.loadtxt( part_file )
 
    x=dt[:,0]; y=dt[:,1];
    vol=dt[:,3]
    w=dt[:,4];
#    vx=dt[:,5]; vym=dt[:,6];

    p=dt[:,9]
    s=dt[:,10]

    # Gallouet & Merrigot
#    p = 0.5*omega**2 * w

    
    r = np.sqrt( x**2 + y**2 )

#make furthest pressure value 0

    rm = np.argmax(r)

    p -= p[ rm ] #  np.min( p )

    
#    I=dt[:,11];  #  I
#    d2=dt[:,12];  #  d2
#    om=dt[:,13];  #  ang velocity

#    I=dt[:,14];  #  eccentricity

#    p += p - np.min( p )

    #    r = np.sqrt( x**2 + y**2 )
    #    plt.plot( r , p , 'o' )

    plt.axis('scaled')
    plt.scatter( x , y , 10, c=p )
#    plt.scatter( x , y , 10, c=s )
#    plt.scatter( x , y , 20, c= vol - 0.000380805861735379 ) # , vmin=0.0022, vmax=0.0028 )
#    plt.scatter( x , y , 10, c=w )
#    plt.scatter( x , y , 10, c=I )
#    plt.scatter( x , y , 80,  c= I , vmin= 1.02e-6, vmax= 1.06e-6 )
#    plt.scatter( x , y , 80,  c= np.log( d2 + 1e-18 ) )
#    plt.scatter( x , y , 10, c=om )

#    di = np.loadtxt(path+str(n)+'/diagram.dat')
    di=np.loadtxt( diag_file )

    xd=di[:,0]; yd=di[:,1];

    for i in range( 0 , xd.size , 2) :
        plt.plot( [ xd[i] , xd[i+1] ] , [ yd[i] , yd[i+1] ] , c='k')

    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')

    plt.xlim([-LL/2.0 , LL/2.0 ])
    plt.ylim([-LL/2.0 , LL/2.0 ])
    plt.colorbar()

    #    pl.colorbar(ticks=[0.45,0.55])

    
    #print( 'snap{:03d}'.format( int(n/skip)  ) )
    #plt.savefig( 'snap{:03d}'.format( int(n/skip)  ) )

    formatted_time = '{:.3f}'.format( float( time ) )
    print( formatted_time )
    plt.savefig( 'snap'+formatted_time+'.png' )
    
