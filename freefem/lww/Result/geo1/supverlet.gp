# phi L2 error
set autoscale
set logscale
plot 'supverlet.dat' us 1:2
replot x**2 with lines
pause 2
set terminal png
set output 'hcvl2phi.png'
replot
set terminal x11
# u1=dx(phi) L2 error
set autoscale
set logscale
plot 'supverlet.dat' us 1:3
replot x with lines
pause 2
set terminal png
set output 'hcvl2u1.png'
replot
set terminal x11
# u2=dy(phi) L2 error
set autoscale
set logscale
plot 'supverlet.dat' us 1:4
replot x with lines
pause 2
set terminal png
set output 'hcvl2u2.png'
replot
set terminal x11
# tr(phi) L2 error
set autoscale
set logscale
plot 'supverlet.dat' us 1:5
replot x**(3./2.),x**2. with lines
pause 2
set terminal png
set output 'hcvl2tra.png'
replot
set terminal x11
# grad(u).N L2 error
set autoscale
set logscale
plot 'supverlet.dat' us 1:6
replot x**(1./2.),x with lines
pause 2
set terminal png
set output 'hcvl2dnu.png'
replot
set terminal x11
# varphi L2 error
set autoscale
set logscale
plot 'supverlet.dat' us 1:7
replot x**(3./2.),x**2. with lines
pause 2
set terminal png
set output 'hcvl2var.png'
replot
set terminal x11
# DNO L2 error
set autoscale
set logscale
plot 'supverlet.dat' us 1:8
replot x**(1./2.),x with lines
pause 2
set terminal png
set output 'hcvl2dno.png'
replot
set terminal x11



unset logscale

set autoscale
plot './mesh1/dt0.01/energy.txt',\
'./mesh2/dt0.01/energy.txt',\
'./mesh3/dt0.01/energy.txt',\
'./mesh4/dt0.01/energy.txt'
pause 2

set terminal png
set output 'energydt0.01.png'
replot
set logscale
set autoscale
set output 'energydt0.01log.png'
replot
unset logscale
set terminal x11

set autoscale
plot './mesh1/dt0.01/errphi.txt',\
'./mesh2/dt0.01/errphi.txt',\
'./mesh3/dt0.01/errphi.txt',\
'./mesh4/dt0.01/errphi.txt'
pause 2

set terminal png
set output 'errphidt0.01.png'
replot
set logscale
set autoscale
set output 'errphidt0.01log.png'
replot
unset logscale
set terminal x11

set autoscale
plot './mesh1/dt0.01/errvar.txt',\
'./mesh2/dt0.01/errvar.txt',\
'./mesh3/dt0.01/errvar.txt',\
'./mesh4/dt0.01/errvar.txt'
pause 2

set terminal png
set output 'errvardt0.01.png'
replot
set logscale
set autoscale
set output 'errvardt0.01log.png'
replot
unset logscale
set terminal x11

set autoscale
plot './mesh1/dt0.01/errdno.txt',\
'./mesh2/dt0.01/errdno.txt',\
'./mesh3/dt0.01/errdno.txt',\
'./mesh4/dt0.01/errdno.txt'
pause 2

set terminal png
set output 'errdnodt0.01.png'
replot
set logscale
set output 'errdnodt0.01log.png'
replot
unset logscale
set terminal x11


set autoscale
plot './mesh1/dt0.01/erreta.txt',\
'./mesh2/dt0.01/erreta.txt',\
'./mesh3/dt0.01/erreta.txt',\
'./mesh4/dt0.01/erreta.txt'
pause 2

set terminal png
set output 'erretadt0.01.png'
replot
set logscale
set output 'erretadt0.01log.png'
replot
unset logscale
set terminal x11

