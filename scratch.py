import numpy as np

# initial parameters
w = 2
theta = - np.pi/2
kernel_ra = (2,2)

# bounding box
x,y = np.meshgrid(np.arange(-kernel_ra[0], kernel_ra[0]+1), np.arange(-kernel_ra[1], kernel_ra[1]+1))

# define parameters based on desired detectable channel width
W = 2*w + 1
sigmax = sigmay = W/(2*(2*np.log(2))**0.5)
f0 = 1/W

# rotation of coordinates
x_ = x*np.cos(theta) + y*np.sin(theta)
y_ = y*np.cos(theta) + x*np.sin(theta)

# real Gabor filter
gb = (1/(2*sigmax*sigmay))*np.exp(-0.5*(x_**2+y_**2)/(sigmax*sigmay))*np.cos(2*np.pi*f0*x_)


import matplotlib.pyplot as plt
fig, ax = plt.subplots()
im = ax.imshow(gb)
fig.colorbar(im)
fig.savefig('test.png')

import os
os.system('nomacs test.png')
