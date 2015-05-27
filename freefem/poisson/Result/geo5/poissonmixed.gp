#u L2 error
set autoscale
set logscale
plot 'hcvmixeddual.dat' us 1:2, \
'hcvmixedprimal.dat' us 1:2 title 'potential'
replot x**2 with lines
pause 3

# dx(u)+dy(u) L2 error
set autoscale
set logscale
plot 'hcvmixeddual.dat' us 1:($3+$4), \
'hcvmixedprimal.dat' us 1:($3+$4) title 'flux'
replot x with lines
pause 3

# tr(u) L2 error
set autoscale
set logscale
plot 'hcvmixeddual.dat' us 1:5, \
'hcvmixedprimal.dat' us 1:5 title 'trace'
replot x with lines
pause 3

# s.N L2 error
set autoscale
set logscale
plot 'hcvmixeddual.dat' us 1:6, \
'hcvmixedprimal.dat' us 1:6 title 'normaltrace'
replot x**2 with lines
pause 3

