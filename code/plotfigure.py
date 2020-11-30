#!/usr/bin/python3

import glob
import os
import numpy as np
from matplotlib import pyplot as plt

fig, ax = plt.subplots()

x, y = [ ], [ ]
for line in open ('DynamicPolyPI6VIII.dat', 'r'):
    values = [float(s) for s in line.split()]

    x.append (values[0])
    y.append (values[2])

ax.plot(x,y, label="VIII")

x2, y2 = [ ], [ ]
for line in open ('dynpol.dat', 'r'):
    values = [float(s) for s in line.split()]
    x2.append (values[0]/2.7)
    y2.append (values[1]*(-2.7)  )

ax.plot(x2, y2, label='Tipsi'  )
ax.legend()
#ax.plot(x,y, label="VIII", x2/2.7, y2*(-2.7), label='Tipsi' )


plt.savefig ('data.eps', format="eps", bbox_inches='tight')
#plt.legend()
plt.show()

