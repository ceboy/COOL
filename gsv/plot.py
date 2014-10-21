#!/usr/bin/python
# -*- coding: utf-8 -*-
# pour lire les resultats
# from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d import axes3d, Axes3D
#from matplotlib.colors import colorConverter
import sys
#if len(sys.argv) != 2:
  #sys.stderr.write("Usage : python %s filetoplot" % sys.argv[0])
  #raise SystemExit(1)
  
if len(sys.argv) == 2:
  mystring = sys.argv[1]
else:
  mystring = 'Nx800/h'

filetoplot = mystring+'.res'

f = open(filetoplot, 'r')
a=(f.readline()).split()
Nx = len(a) # nombre de cellules, en abcisse
t = 0 # compteur temps (= nombre de lignes)
verts = []
mymin = 0
mymax = 0
zz = []
while(len(a)>1): # while: lecture des lignes (debut)
  u = [float(x) for x in a]
  #verts.append(zip(np.arange(mystep*.5,1.,mystep),u))
  verts.append(zip(range(0,Nx),u))
  zz.append(u)
  a=(f.readline()).split()
  t+=1 # on a lu une nouvelle ligne
  mymin = min([mymin,min(u)])
  mymax = max([mymax,max(u)])

f.close() # while: lecture des lignes (fin)  

Z = np.array(zz)
X = np.empty((t,Nx))
Y = np.empty((t,Nx))
tt = 0
mystep = 10./float(Nx)
while(tt<t):
  X[tt] = np.arange(-5.+mystep*.5,5.,mystep)
  Y[tt] = np.ones((1,Nx))*tt
  tt += 1

#plt.plot(*zip(*verts[1]))

from matplotlib import cm
fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
ax = Axes3D(fig)
#surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1)
#surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.get_cmap(1),
#linewidth=0, antialiased=False)
#fig.colorbar(surf, shrink=0.5, aspect=5)
#ax.plot_wireframe(X, Y, Z)
ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
# ax = fig.gca(projection='3d')
#ax.bar(range(0,Nx), rang(0,t), zs=z, zdir='y', color=cs, alpha=0.8)
#cc = lambda arg: colorConverter.to_rgba(arg, alpha=0.6)
# poly = PolyCollection(verts, facecolors = [cc('g')]) # cc('r'), cc('g'), cc('g'), cc('g')]
#poly = PolyCollection(verts)
#poly.set_alpha(0.7)
#ax.add_collection3d(poly, zs=range(0,t), zdir='y')
ax.set_xlabel('X')
ax.set_xlim3d(-5.,5.)
ax.set_ylabel('Y')
ax.set_ylim3d(0, t-1)
ax.set_zlabel('Z')
#ax.set_zlim3d(-1, 1)
ax.set_zlim3d(mymin, mymax)

plt.savefig(mystring+'3d.png')

plt.show() # parfois inutile

#   reset


