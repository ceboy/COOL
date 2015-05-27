set autoscale
set logscale
# phi L2 error
plot 'potential.dat' us 1:2,\
'potentialmixedprimal.dat' us 1:2
replot x**2 with lines
pause 3
# u1=dx(phi) L2 error
plot 'potential.dat' us 1:3,\
 'potentialmixedprimal.dat' us 1:3
replot x with lines
pause 3
# u2=dy(phi) L2 error
plot 'potential.dat' us 1:4,\
 'potentialmixedprimal.dat' us 1:4
replot x with lines
pause 3
# tr(phi) L2 error
plot  'potential.dat' us 1:5,\
'potentialmixedprimal.dat' us 1:5
replot x**(3./2.),x**2. with lines
pause 3
# grad(u).N L2 error
plot 'potential.dat' us 1:6,\
'potentialmixedprimal.dat' us 1:6
replot x**(1./2.),x with lines
pause 3
# varphi L2 error
plot 'potential.dat' us 1:7,\
'potentialmixedprimal.dat' us 1:7
replot x**(3./2.),x**2. with lines
pause 3
# DNO L2 error
plot 'potential.dat' us 1:8,\
'potentialmixedprimal.dat' us 1:8
replot x**(1./2.),x with lines
pause 3


#unset logscale

plot 'Result/geo1mesh1/dt0.01/energymp.txt',\
'Result/geo1mesh2/dt0.01/energymp.txt',\
'Result/geo1mesh3/dt0.01/energymp.txt',\
'Result/geo1mesh4/dt0.01/energymp.txt',\
'Result/geo1mesh1/dt0.01/energy.txt',\
'Result/geo1mesh2/dt0.01/energy.txt',\
'Result/geo1mesh3/dt0.01/energy.txt',\
'Result/geo1mesh4/dt0.01/energy.txt'
pause 3

set autoscale
plot 'Result/geo1mesh1/dt0.01/errphimp.txt',\
'Result/geo1mesh2/dt0.01/errphimp.txt',\
'Result/geo1mesh3/dt0.01/errphimp.txt',\
'Result/geo1mesh4/dt0.01/errphimp.txt'
pause 3

set autoscale
plot 'Result/geo1mesh1/dt0.01/errvarmp.txt',\
'Result/geo1mesh2/dt0.01/errvarmp.txt',\
'Result/geo1mesh3/dt0.01/errvarmp.txt',\
'Result/geo1mesh4/dt0.01/errvarmp.txt'
pause 3

set autoscale
plot 'Result/geo1mesh1/dt0.01/errdnomp.txt',\
'Result/geo1mesh2/dt0.01/errdnomp.txt',\
'Result/geo1mesh3/dt0.01/errdnomp.txt',\
'Result/geo1mesh4/dt0.01/errdnomp.txt'
pause 3

set autoscale
plot 'Result/geo1mesh1/dt0.01/erretamp.txt',\
'Result/geo1mesh2/dt0.01/erretamp.txt',\
'Result/geo1mesh3/dt0.01/erretamp.txt',\
'Result/geo1mesh4/dt0.01/erretamp.txt'
pause 3
