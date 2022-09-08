#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
#import matplotlib.cm as cm
#from matplotlib.colors import Normalize

import sys
import glob

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

def pp(x , y) : # analytic solution for TG's vortex pressure

    return  ( np.cos( 4 * np.pi * x ) +   np.cos( 4 * np.pi * y ) ) / 4

v_pp = np.vectorize( pp, otypes=[float] )

def by_number(elem):
    return float(elem)

times = sorted(glob.glob('[0-9]*'), key=by_number)

first_time=times[0]

den=1

for time in times:

    part_file = time+'/particles.dat'

    dt=np.loadtxt( part_file )

    x=dt[:,0]; y=dt[:,1];
    #    vol=dt[:,3]
    w=dt[:,4];
    #    vx=dt[:,5]; vym=dt[:,6];
    p=dt[:,9]
    #  s=dt[:,10]
    #  I=dt[:,11];

    # Gallouet & Merrigot
    #p = 0.5*omega**2 * w

    r = np.sqrt( x**2 + y**2 )

    #make furthest pressure value 0

#    rm = np.argmax(r)

 #   p -= p[ rm ] #  np.min( p )

    p -= np.mean( p )

    L2 = np.linalg.norm(  p  - v_pp( x , y ) )

    if time == first_time :
        den = np.linalg.norm(  v_pp( x , y) )

#    print(v_pp(r))
        
    print(time, L2 / den )
