#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 10:20:17 2018

@author: Leo
"""

import matplotlib.pyplot as plt
import numpy as np

dx, dy = 0.01, 0.01
y, x = np.mgrid[slice(0, 1 + dy, dy),
                slice(0, 2 + dx, dx)]
#print(x[0,0])
f = open('solution.dat', 'r')
data = np.genfromtxt(f, delimiter=',')
data=np.delete(data,-1,1) # Erases the last column
data=np.delete(data,-1,1) # Erases the last column
data=np.delete(data,-1,0) # Erases the last row
f.close()
fig, ax = plt.subplots()
plt.pcolor(x, y, np.transpose(data), cmap='RdBu_r')#reversed colormap
plt.title('Steady state heat distribution')
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.colorbar()
plt.savefig("Heated tile.pdf",format='pdf',bbox_inches='tight')
plt.show()