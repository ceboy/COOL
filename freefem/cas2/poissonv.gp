#u1 L2 error
set autoscale
set logscale
plot 'hcv.dat' us 1:2
replot x**2 with lines
pause 2
#u2 L2 error
set autoscale
set logscale
plot 'hcv.dat' us 1:3
replot x**2 with lines
pause 2
# dx(u1) L2 error
set autoscale
set logscale
plot 'hcv.dat' us 1:4
replot x with lines
pause 2
# dy(u1) L2 error
set autoscale
set logscale
plot 'hcv.dat' us 1:5
replot x with lines
pause 2
# dx(u2) L2 error
set autoscale
set logscale
plot 'hcv.dat' us 1:6
replot x with lines
pause 2
# dy(u2) L2 error
set autoscale
set logscale
plot 'hcv.dat' us 1:7
replot x with lines
pause 2
# tr(u1) L2 error
set autoscale
set logscale
plot 'hcv.dat' us 1:8
replot x**(3./2.),x**2. with lines
pause 2
# tr(u2) L2 error
set autoscale
set logscale
plot 'hcv.dat' us 1:9
replot x**(3./2.),x**2. with lines
pause 2
# s1.N L2 error
set autoscale
set logscale
plot 'hcv.dat' us 1:10
replot x**(1./2.),x with lines
pause 2
# s2.N L2 error
set autoscale
set logscale
plot 'hcv.dat' us 1:11
replot x**(1./2.),x with lines
pause 2
