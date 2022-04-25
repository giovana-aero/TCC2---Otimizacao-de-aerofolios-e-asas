# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 19:04:51 2022

@author: Guga Weffort
"""

import numpy as np
import matplotlib.pyplot as plt
from cosspace import cosspace
from cosspace_half import cosspace_half


x1 = 0
x2 = 1
# N = 100;

# x = cosspace_half(x1,x2)
# x = cosspace_half(x1,x2,N)

# x = cosspace(x1,x2)
x = cosspace(x1,x2,20)

# x = np.array(x)
y = x**2 + 2*x + 2
print(x)

plt.clf()
plt.plot(x,y)
plt.grid(True)
# plt.hold() # Desnecessário, parece?
plt.scatter(x,y)