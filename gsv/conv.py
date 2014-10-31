#!/usr/bin/python
# -*- coding: utf-8 -*-
# use :reset :run in ipython

import numpy as np
import matplotlib.pyplot as plt

mystring = 'Q'
mylist1 = [1, 2, 4, 8]
myfiles = map(lambda x:str(x)+'00/'+mystring+'.res',mylist1)
#myfiles = map(lambda x:'Nx'+str(x)+'00/Q.res',mylist1)
table = []
nbfiles = 0
nblines = 0
for myfile in myfiles:
  f = open(myfile,'r')
  tabletemp = []
  for myline in f:
    tabletemp.append( map(float,myline.split()) )
    if(nbfiles==0):
      nblines += 1 # to count the number of snapshots in time
  #print tabletemp # nblines: counter in list f.readlines()
  table.append(tabletemp) # table: list of lists
  nbfiles += 1 # nbfiles: counter in lists mylist1 and myfiles
  #print nbfiles

print '"table" contains ', nbfiles, ' files of ', nblines, ' lines = snapshots in time'

# table[i][j] is the vector solution for cfl=i and times=j

#x = np.linspace(0, 20, 1000)
#plt.xlim(5, 15)
#plt.ylim(-1.2, 1.2)
#plt.xlabel('this is x!')
#plt.ylabel('this is y!')
#plt.title('My First Plot')
#plt.title(r'$\sin(2 \pi x)$') # the `r` before the string indicates a "raw string"
#plt.plot(x, y, '-r', label='sine')  # solid red line ('r' comes from RGB color scheme)
#Other options for the color characters are:
 #'r' = red
 #'g' = green
 #'b' = blue
 #'c' = cyan
 #'m' = magenta
 #'y' = yellow
 #'k' = black
 #'w' = white
#Options for line styles are
 #'-' = solid
 #'--' = dashed
 #':' = dotted
 #'-.' = dot-dashed
 #'.' = points
 #'o' = filled circles
 #'^' = filled triangles
 #and many, many more.
#For more information, view the documentation of the plot function. In IPython, this can be accomplished using the ? functionality:
#plt.plot?
#http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
#plt.legend(loc='upper right')

ntime = nblines-1
# ntime = nblines/2
tabshow = []
for ii in range(0,nbfiles):
  myline = table[ii][ntime]
  mylen = len(myline)
  mystep = 10./float(mylen)
  # map(lambda x: mystep*.5+x*mystep, range(0, mylen))
  tabshow.append(zip(np.arange(mystep*.5,1.,mystep),myline))
  myarray = np.array(myline)
  # plt.plot(np.arange(mystep*.5,1.,mystep),myline,label = 'Fluid depth h using '+mylist1(ii)+'00 cells')
  plt.plot(np.arange(-5.+mystep*.5,5.,mystep),myline)
  plt.axis([-5., 5., 0.,np.amax(myarray)])

#plt.show() # in this figure, we hope to see the approximation order with respect to "dx", on superimposing various time snapshots
plt.savefig(mystring+'.png')

plt.clf()

errtab = []
reftab = np.array(table[nbfiles-1]) # reference solution with 800 points = columns at each times = lines
Reftab = np.transpose(reftab) 
ttab = np.empty([nbfiles-1,nblines])
n_factor = 1
for nbfile in range(nbfiles-2,-1,-1):
  n_factor *= 2
  curtab = np.array(table[nbfile])
  Curtab = np.transpose(curtab) 
  n_rows = Curtab.shape[0] # points
  n_columns = Curtab.shape[1] # times
  Curtab2 = np.empty([n_rows*n_factor,n_columns])
  i = 0
  for n_row in range(n_rows):
    for factor in range(n_factor):
      Curtab2[i] = Curtab[n_row]
      i +=1
  Errtab = Reftab-Curtab2 # compute absolute pointwise error
  #print Errtab # errtab[i][j] is error at point i and time j 
  ttab[nbfile] = np.sum(np.abs(Errtab),axis=0)*.01 # L1 space-norm of absolute error: sum over points i
  print ttab[nbfile] # print the absolute errors for file #(nbfile) at various times

Ttab = np.transpose(ttab)
#fig = plt.figure()


#for nbline in range(1,nblines-1):
  #plt.plot(np.log(Ttab[nbline][nbfiles-2:0:-1]))

#plt.show() # in this figure, we hope to see the approximation order with respect to "dx", on superimposing various time snapshots


#for nbline in range(1,3):
  #plt.plot(np.log(Ttab[nbline][nbfiles-2:0:-1]))

#plt.show() # in this figure, we see order 1 with respect to "dx" for initial times (1,3) only


#for nbline in range(11,20):
  #plt.plot(np.log(Ttab[nbline][nbfiles-2:0:-1]))

#plt.show() # in this figure, we see order 1 with respect to "dx" for intermediate times (11,20) 


#for nbline in range(25,99):
  #plt.plot(np.log(Ttab[nbline][nbfiles-2:0:-1]))

#plt.show() # in this figure, we see order 1/2 with respect to "dx" for intermediate times (25,99)


#for nbline in [0,51]:
  #plt.plot(np.log(Ttab[nbline][nbfiles-2:0:-1]))

#plt.show() # in this figure, we see order 1/2 with respect to "dx" for times 0 and 51

plt.plot(np.log(np.sum(Ttab,axis=0)/len(Ttab))[nbfiles-2:0:-1])
plt.show() # in this figure, we hope to see the approximation order with respect to "dx", on superimposing various time snapshots

# we observe 1-order convergence 

#raise SystemExit
