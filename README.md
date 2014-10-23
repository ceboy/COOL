COOL
====

/gsv/ Fortran code for generalized Saint-Venant equations

  :$ cp /home/SB03743S/COOL/gsv/* .; make cleanall; make

  :$  ./sv.exe 100; ./plot.py phi; display phi3d.png

  :$ for NX in 100 200 400 800; do mkdir ${NX}; cp ./sv.exe ${NX}; cd ${NX}; ./sv.exe ${NX}; cd ..; done > sv.log &

  :$ for NX in 100 200 400 800; do rm -rf ${NX}; mkdir ${NX}; cp ./sv.exe ${NX}; cd ${NX}; ./sv.exe ${NX}; cd ..; done > sv.log &

  :$ cd 800; ./plot.py phi; display phi3d.png

