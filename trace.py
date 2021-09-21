#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
#import matplotlib.cm as cm
#from matplotlib.colors import Normalize

#import pylab as pl
import glob

import sys

if(len(sys.argv) == 1) :
    idx = 126
else:
    idx = int( sys.argv[1] )

#print( idx )

skip=1
#path='timings_full/'
path='./'

LL= 1


T_spring = 20 * 0.005
omega = 2 * np.pi / T_spring

list_time = glob.glob('[0-9]*')

#  sort
def by_number(elem):
    return float(elem)

times = sorted(glob.glob('[0-9]*'), key=by_number)

for time in times:

    dt=np.loadtxt(time+'/particles.dat')

    x=dt[:,0]
    y=dt[:,1]

    w=dt[:,4];
    #    vx=dt[:,5]; vym=dt[:,6];
    #  s=dt[:,10]
    #  I=dt[:,11];

    p=dt[:,9]
    # Gallouet & Merrigot only
    #p = 0.5*omega**2 * w

    r = np.sqrt( x**2 + y**2 )

    #make furthest pressure value 0

    rm = np.argmax(r)

    p -= p[ rm ] #  np.min( p )

    print( time , x[idx] , y[idx] , p[idx])

    
