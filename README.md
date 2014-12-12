COOL
====

/gsv/ contains : **Fortran code for generalized Saint-Venant equations with a few testcases**

 * sv            **Main containing the time loop**
 * svini         **Subroutine initializing data for test case**
 * m_zeromachine **Module where the (double precision)myzeromachine "global" variable is declared (without dependency) ; it is assigned in main (sv)**
 * m_physics     **Module where some "global" variables are declared (constants: Pi, gravity; friction, wave number) and assigned (without dependency)**
 * m_cell        **Module for the declaration of the new structure of type t_cell (without dependency)**
 * m_data        **Module where "global" data is declared (requires m_cell)**

>> to upload the files in a directory where to run a testcase

  :$ cp /home/SB03743S/COOL/gsv/* .; make cleanall; make

>> to run a testcase and display the component phi of the solution

  :$  ./sv.exe 100; ./plot.py phi; display phi3d.png

>> to run a convergence testcase in subdirectories

  :$ for NX in 100 200 400 800; do mkdir ${NX}; cp ./sv.exe ${NX}; cd ${NX}; ./sv.exe ${NX} > sv.log; cd ..; done ; ./conv.py &

testcases ::

  (Augmented) standard Saint-Venant equations

  1. :: Initially-at-rest dam-break with Stoker solution on a "wet" bed (1-rarefaction + 2-shock)
  2. :: Initially-at-rest dam-break with Ritter solution on a "dry" bed (1-rarefaction + 2-shock without intermediate)
  3. :: Initially-at-rest doubled dam-break (opposed, one on the left on on the right, with vacuum in between)
  4. :: Initially-at-rest column of water 
  5. :: Initially-at-rest tripled dam-break (intercation of shocks: curves)
  6. :: Initially-at-rest sinusoidal wave of the free-surface (straightening, no dipersion)
  7. :: Initially-at-rest dam-break \sim Stoker solution on a "wet" bed + viscoelastic components (a G=1; b G=10)
  8. :: Initially-at-rest dam-break \sim Ritter solution on a "dry" bed + viscoelastic components (a G=1; b G=10)
  9. :: Initially-at-rest dam-break \sim Stoker solution on a "dry" bed + viscoelastic components including relaxation to equilibrium 
          (a G=1; b G=10) lambda = 1/1 (c G=1; d G=10) lambda = 1/10 (e G=1; f G=10) lambda = 1/100 (g G=1; h G=10) lambda = 1/1000
  10. :: Initially-at-rest dam-break \sim Ritter solution on a "dry" bed + FENEP components including relaxation to equilibrium 
  11. :: Test to check HWNP (limit G->0 as lambda->0) beware 800: until 1. s !!


TO DO:

* add dynamics:
  - FENE P (PEC, Giesekus, PTT)
  - turbulent
* add source terms:
  - topography, introduce theta in svini and use g0=g*cos(theta) in sv.f90 etc...)
  - friction
  - surface tension
  - wind
* 2D version
  - structured
  - unstructured (mesh definition !)
