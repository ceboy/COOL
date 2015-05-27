#! /usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import pylab
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

#if len(sys.argv)!=2:
  #sys.stderr.write("Usage : python %s filetoplot" % sys.argv[0])
  #raise SystemExit(1)
#filetoplot = sys.argv[1]

Nx = 30
T = 20
K0 = 1
sigma = 0.1
Kx = 1
filetoplot = 'res'+str(Nx)+'T'+str(T)+'K0'+str(K0)+'sigma'+str(sigma)+'Kx'+str(Kx)

mytimes = np.loadtxt(filetoplot+'.log')
#print(mytimes)
Nt = mytimes.shape[0]

mysol = np.loadtxt(filetoplot+'.txt')

fig, ax = plt.subplots()

#x = np.arange(0, 2*np.pi, 2*np.pi/Nx)
#line, = ax.plot(x, np.sin(x))

x = np.arange(np.pi/Nx, 2*np.pi, 2*np.pi/Nx)

line, = ax.plot(x,mysol[0],'*')

#plt.show()

plt.savefig('initialstate.pdf')

fig3d = plt.figure()
ax3d = Axes3D(fig3d)
# ax3d = fig3d.gca(projection='3d')
for i in range(Nt):
    xx=list([])
    yy=list([])
    zz=list([])
    for j in range(Nx):
        xx.append(mytimes[i][0])
        yy.append(x[j])
        zz.append(mysol[i,j])
    ax3d.plot(xx,yy,zz,'b-')

ax3d.set_xlabel('t')
ax3d.set_xlim3d(0.,T)
ax3d.set_ylabel('x')
ax3d.set_ylim3d(0.,1.)
ax3d.set_zlabel('u')
ax3d.set_zlim3d(-1., 1.)

plt.savefig('allstates3d.pdf')

def animate(i):
#    line.set_ydata(np.sin(x+i/10.0))  # update the data
    line.set_ydata(list(mysol[i]))
    return line,

#Init only required for blitting to give a clean slate.
def init():
    line.set_ydata(np.ma.array(x, mask=True))
    return line,

ani = animation.FuncAnimation(fig, animate, np.arange(1, Nt), init_func=init, interval=35, blit=True)


#ax = fig.add_subplot(111, projection='3d')

plt.show()

# reset




