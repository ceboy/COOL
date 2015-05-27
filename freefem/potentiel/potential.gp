# phi L2 error
set autoscale
set logscale
plot 'potential.dat' us 1:2
replot x**2 with lines
pause 3
# u1=dx(phi) L2 error
set autoscale
set logscale
plot 'potential.dat' us 1:3
replot x with lines
pause 3
# u2=dy(phi) L2 error
set autoscale
set logscale
plot 'potential.dat' us 1:4
replot x with lines
pause 3
# tr(phi) L2 error
set autoscale
set logscale
plot 'potential.dat' us 1:5
replot x**(3./2.),x**2. with lines
pause 3
# grad(u).N L2 error
set autoscale
set logscale
plot 'potential.dat' us 1:6
replot x**(1./2.),x with lines
pause 3
# varphi L2 error
set autoscale
set logscale
plot 'potential.dat' us 1:7
replot x**(3./2.),x**2. with lines
pause 3
# DNO L2 error
set autoscale
set logscale
plot 'potential.dat' us 1:8
replot x**(1./2.),x with lines
pause 3



set autoscale
plot 'Result/geo1mesh1/dt0.01/energy.txt',\
'Result/geo1mesh2/dt0.01/energy.txt',\
'Result/geo1mesh3/dt0.01/energy.txt',\
'Result/geo1mesh4/dt0.01/energy.txt'
pause 3

set autoscale
plot 'Result/geo1mesh1/dt0.01/errphi.txt',\
'Result/geo1mesh2/dt0.01/errphi.txt',\
'Result/geo1mesh3/dt0.01/errphi.txt',\
'Result/geo1mesh4/dt0.01/errphi.txt'
pause 3

set autoscale
plot 'Result/geo1mesh1/dt0.01/errvar.txt',\
'Result/geo1mesh2/dt0.01/errvar.txt',\
'Result/geo1mesh3/dt0.01/errvar.txt',\
'Result/geo1mesh4/dt0.01/errvar.txt'
pause 3

set autoscale
plot 'Result/geo1mesh1/dt0.01/errdno.txt',\
'Result/geo1mesh2/dt0.01/errdno.txt',\
'Result/geo1mesh3/dt0.01/errdno.txt',\
'Result/geo1mesh4/dt0.01/errdno.txt'
pause 3

set autoscale
plot 'Result/geo1mesh1/dt0.01/erreta.txt',\
'Result/geo1mesh2/dt0.01/erreta.txt',\
'Result/geo1mesh3/dt0.01/erreta.txt',\
'Result/geo1mesh4/dt0.01/erreta.txt'
pause 3
