#u L2 error
set autoscale
set logscale
plot 'hcv.dat' us 1:2
replot x**2 with lines
pause 2
# dx(u) L2 error
set autoscale
set logscale
plot 'hcv.dat' us 1:3
replot x with lines
pause 2
# dy(u) L2 error
set autoscale
set logscale
plot 'hcv.dat' us 1:4
replot x with lines
pause 2
# tr(u) L2 error
set autoscale
set logscale
plot 'hcv.dat' us 1:5
replot x**(3./2.),x**2. with lines
pause 2
# s.N L2 error
set autoscale
set logscale
plot 'hcv.dat' us 1:6
replot x**(1./2.),x with lines
pause 2
