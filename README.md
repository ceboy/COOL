COOL
====

/gsv/ :: Fortran code for generalized Saint-Venant equations with a few testcases

 ** sv.f90 :: main with the time loop
 ** svini.f90 :: to initialize the global data accessible through m_data (see testcases)

>> to upload the files in a directory where to run a testcase

  :$ cp /home/SB03743S/COOL/gsv/* .; make cleanall; make

>> to run a testcase and display the component phi of the solution

  :$  ./sv.exe 100; ./plot.py phi; display phi3d.png

>> to run a convergence testcase in subdirectories

  :$ for NX in 100 200 400 800; do mkdir ${NX}; cp ./sv.exe ${NX}; cd ${NX}; ./sv.exe ${NX} > sv.log; cd ..; done ; ./conv.py &

testcases ::

  (Augmented) standard Saint-Venant equations

  ** 1 :: Initially-at-rest dam-break with Stoker solution on a "wet" bed (1-rarefaction + 2-shock)
  ** 2 :: Initially-at-rest dam-break with Ritter solution on a "dry" bed (1-rarefaction + 2-shock without intermediate)
  ** 3 :: Initially-at-rest doubled dam-break (opposed, one on the left on on the right, with vacuum in between)
  ** 4 :: Initially-at-rest column of water 
  ** 5 :: Initially-at-rest tripled dam-break (intercation of shocks: curves)
  ** 6 :: Initially-at-rest sinusoidal wave of the free-surface (straightening, no dipersion)
  ** 7 :: Initially-at-rest dam-break \sim Stoker solution on a "wet" bed + viscoelastic components (a G=1; b G=10)
  ** 8 :: Initially-at-rest dam-break \sim Ritter solution on a "dry" bed + viscoelastic components (a G=1; b G=10)


