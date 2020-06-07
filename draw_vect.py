#!/usr/bin/python

import sys

#print ("This is the name of the script: " ,sys.argv[0] )
#print ("Number of arguments: ", len(sys.argv) )
#print "The arguments are: " , str(sys.argv)

if ( len(sys.argv) != 2) :
    print('Usage: ' , sys.argv[0] , '  file_name ' )
    sys.exit()

file = sys.argv[1] + '/particles.dat'

print("From file " , file)

#import pylab as pl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import Normalize

dt = np.loadtxt( file )

x=dt[:,0]
y=dt[:,1]
vol=dt[:,3]
w=dt[:,4]
ux=dt[:,5]
uy=dt[:,6] 

colors = w
norm = Normalize()
norm.autoscale(colors)

colormap = cm.inferno

plt.figure(figsize=(6, 6))

plt.quiver(x,y, ux, uy , color=colormap(norm(colors)) ,angles='xy', 
           scale_units='xy', scale=10, pivot='mid' )
#plt.xlim(-.4,.4)
#plt.ylim(-.4,.4)

#plt.colorbar()

plt.show()

