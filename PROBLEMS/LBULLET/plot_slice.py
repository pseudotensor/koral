#Plot slice.dat data

import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate as intp
from math import *

slicefile = open('slice.dat','r')

slicedata = []

for line in slicefile:
	slicedata.append(map(float,line.split()))
	if len(slicedata[-1])!=5:
		slicedata.pop(-1)

points = [row[0:2] for row in slicedata]
points_x = [row[0] for row in slicedata]
points_y = [row[1] for row in slicedata]
values = [log(row[2]) if row[2] > 0 else 0 for row in slicedata]  


grid_x, grid_y = np.mgrid[1:500:100j,1:500:100j]


grid = intp.griddata(points,values,(grid_x,grid_y),method='linear')
print values

plot = plt.pcolormesh(grid_x,grid_y,grid,vmin = -10, vmax = 0)
plt.colorbar(plot)
#plt.plot(values)
plt.show()
