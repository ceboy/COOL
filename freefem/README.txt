/freefem/ contains

/poisson/:

  mixed Dirichlet-Neumann BVP for the scalar poisson equation "-Delta u = f" solved in various polygonal geometries with data from various analytical solutions

  essential BC imposed by substitution of the corresponding DOF

  primal (P1: CG can be improved by preconditionning),
  primal-mixed (P1-P0^2: equivalent to primal after static condensation) and
  dual-mixed formulations (P0-RT0: GmRes for SADDLE-POINT could be replaced by STATIC CONDENSATION -- equal to P1nc (Crouzeix-Raviart) here ! -- or UZAWA)

/cas2/:

  vector Laplace equation solved in
    //geo1 (a square) uniformly meshed
    geo3 (a triangle) uniformly meshed
  given
    //Dirichlet BC for geo1
    Dirichlet BC for geo3
  and a source term
    //f for geo1
    f for geo3

/cas3/:

  vector Laplace equation solved in
    geo3 (a triangle) uniformly meshed
  given
    mixed Neumann-Dirichlet BC for geo3
  and a source term
    f for geo3

/potential/:

  linear water waves problem solved with Velocity Verlet and various variational formulations in space (see poisson above)

/lww/:

  linear water waves problem solved with Velocity Verlet in primal formulation







for(int i=0;i<ARGV.n;++i)
  {
    cout << ARGV[i] << endl;
  }


