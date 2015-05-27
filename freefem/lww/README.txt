Linear Water Wave problem "d_t eta = DNO(phi)" for eta surface of a fluid domain where "-Delta phi = 0", "d_t phi = -gravity*eta"

  BVP for the scalar Laplace equation in a fixed domain with mixed Dirichlet-Neumann-Periodic BC

  solved in various polygonal geometries (also various meshes for each geomtry) with

  initial data (from various analytical solutions) and

  essential BC (Dirichlet/Neumann and half-Periodic) imposed by substitution of the corresponding DOF in various variational formulations, conforming :

  primal (using P1: CG for SPD matrix can be improved by preconditionning),

=>primal-mixed (P1-P0^2: undefinite saddle-point equivalent to primal after static condensation can be solved by UMFPACK or GMRES without penalization)

  dual-mixed  (P0-RT0: undefinite saddle-point can be solved by UMFPACK or GMRES without penalization, it is equivalent to non-conforming P1nc)

 PRE:
  * input.idp: mesh, BC (type and value), FE, solver and linalg parameters

 RUN:
  * verlet*.edp: according to the variational formulation of the problem retained

 POST:
  * verlet*.gp: gnuplot script to show *.dat

 TODO:
  * see Poisson for the Poisson sub-problem
  * other time-algorithms
  * BC moving in time ?
