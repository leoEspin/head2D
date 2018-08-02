#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Script for plotting solution
Created on Wed Aug  1 10:20:17 2018

@author: Leo
"""

import matplotlib.pyplot as plt
import numpy as np

dx, dy = 0.01, 0.01
y, x = np.mgrid[slice(0, 2 + dx, dx),
                slice(0, 1 + dy, dy)]
f = open('solution.dat', 'r')
data = np.genfromtxt(f, delimiter=',')
data=data[0:-1,0:-2]
f.close()
fig, ax = plt.subplots()
plt.pcolor(x, y, data, cmap='RdBu_r')#reversed colormap
plt.title('Steady state heat distribution')
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.colorbar()
plt.savefig("Heated tile.pdf",format='pdf',bbox_inches='tight')
plt.show()
