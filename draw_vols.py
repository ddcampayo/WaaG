#!/usr/bin/python

import sys


#print ("This is the name of the script: " ,sys.argv[0] )
#print ("Number of arguments: ", len(sys.argv) )
#print "The arguments are: " , str(sys.argv)

if ( len(sys.argv) != 2) :
    print('Usage: ' , sys.argv[0] , '  directory for particles.dat and diagram.dat' )
    sys.exit()

path = sys.argv[1]

file = path + '/particles.dat'

print("From file " , file)

#import pylab as pl
import matplotlib.pyplot as plt
import numpy as np

dt = np.loadtxt( file )

x=dt[:,0]
y=dt[:,1]
iid=dt[:,2]
vol=dt[:,3]
w=dt[:,4] 
I=dt[:,11] 
dd2=dt[:,12] 
om=dt[:,13] 
I4=dt[:,14] 

#plt.scatter(x,y, 100*np.sqrt(vol) , c= w)
#plt.scatter(x,y, 100*np.sqrt(vol) , c= I)
plt.scatter(x,y, c= I)
#plt.scatter(x,y, c= np.sqrt( dd2 ) )
#plt.scatter(x,y, 10,  c= om)
#plt.scatter(x,y, 10,  c= I4)
#plt.scatter(x,y, 20,  c= vol)
#plt.scatter(x,y, 20,  c= iid)

ll = 0.6

plt.xlim(-ll,ll)
plt.ylim(-ll,ll)

plt.colorbar()



di = np.loadtxt(path+'/diagram.dat')

xd=di[:,0]; yd=di[:,1];

for i in range( 0 , xd.size , 2) :
    plt.plot( [ xd[i] , xd[i+1] ] , [ yd[i] , yd[i+1] ] , c='k')



plt.show()

