# -*- coding: utf-8 -*-
"""
Created on Sun Dec  5 00:23:55 2021

@author: Guga Weffort
"""

# import csv
import numpy as np
import matplotlib.pyplot as plt

# Ler coordenadas
name = 'coordenadas.dat'
coo = np.genfromtxt(name)

# Fazer o gráfico
fig = plt.figure()
ax1 = plt.axes()
ax1.plot(coo[:,0],coo[:,1],'k')
ax1.axis('equal')
ax1.grid(True)