import sys

#print ("This is the name of the script: " ,sys.argv[0] )
#print ("Number of arguments: ", len(sys.argv) )
#print "The arguments are: " , str(sys.argv)

if ( len(sys.argv) != 2) :
    print('Usage: ' , sys.argv[0] , '  file_name ' )
    sys.exit()

file = sys.argv[1]

print("From file " , file)

#import pylab as pl
import matplotlib.pyplot as plt
import numpy as np

dt = np.loadtxt( file )

x=dt[:,0]
y=dt[:,1]
vol=dt[:,3]
w=dt[:,4] 

plt.scatter(x,y, 100*np.sqrt(vol) , c= w)

ll = 0.6

plt.xlim(-ll,ll)
plt.ylim(-ll,ll)

plt.colorbar()


plt.show()

