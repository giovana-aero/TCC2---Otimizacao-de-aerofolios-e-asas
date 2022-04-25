from naca import naca4
import numpy as np
import matplotlib.pyplot as plt
from airfoil_interpolation import airfoil_interpolation
from airfoil_rotation import airfoil_rotation

naca1 = '0012'
naca2 = '9530'

x1,y1 = naca4(naca1,100)
x2,y2 = naca4(naca2,100)

coo1 = np.hstack((np.array([x1]).transpose(),np.array([y1]).transpose()))
# coo2 = np.hstack((np.array([x2]).transpose(),np.array([y2]).transpose()))
# coo3 = airfoil_interpolation(coo1,coo2,0.5)
coo2 = airfoil_rotation(coo1,-5)

plt.figure()
plt.plot(x1,y1)
# plt.plot(x2,y2)
# plt.plot(coo3[:,0],coo3[:,1])
plt.plot(coo2[:,0],coo2[:,1])
plt.axis('equal');plt.grid('true')


