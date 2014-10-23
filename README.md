COOL
====

/gsv/ Fortran code for generalized Saint-Venant equations

  :$ cp /home/SB03743S/COOL/gsv/* .; make cleanall; make

  :$  ./sv.exe 100; ./plot.py phi; display phi3d.png

  :$ for NX in 100 200 400 800; do mkdir ${NX}; cp ./sv.exe ${NX}; cd ${NX}; ./sv.exe ${NX} > sv.log; cd ..; done &

  :$ rm -rf; for NX in 100 200 400 800; do mkdir ${NX}; cd ${NX}; cp /home/SB03743S/COOL/gsv/* .; make cleanall; make; ./sv.exe ${NX} > sv.log; cd ..; done &

  :$ cd 800; ./plot.py phi; display phi3d.png

