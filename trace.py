#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
#import matplotlib.cm as cm
#from matplotlib.colors import Normalize

#import pylab as pl


skip=1
#path='timings_full/'
path='./'

LL= 1

idx = 126

for n in range(0,2000000+skip,skip):

    dt=np.loadtxt(path+str(n)+'/particles.dat')

    print( 0.005*n , dt[idx,0] ,  dt[idx,1] )
