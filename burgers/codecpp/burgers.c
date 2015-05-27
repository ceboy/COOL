#include "mtrand.cpp"  // Mersenne-Twister
#include <cmath>       // floor M_PI
#include <fstream>
#include <sstream>
#include <cstdio>
//#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstring>
//#include <string.h>
#include <string>
#include <sstream>
#include <cstdlib>    // atoi atof
//#include <stdlib.h>
// #include <vector>
// #include <list>
// #include <assert.h>     /* assert */
// #include <mpi.h>
#include <ctime>  // clock()

using namespace std;
//std::ios::sync_with_stdio(false);// to accelerate iostream
// #ifndef M_PI
//    #define M_PI 3.14159265358979323846
// #endif
// #define MPI_MASTER 0

double FLUX(double localstate)
{
  return (localstate*localstate)/2.;
}

int cmp(const void *x, const void *y)
{
  double xx = *(double*)x, yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
}

int main(int argc, char* argv[])
{
  setprecision(1.1);           // set parameter of cout
/*
  int rank, size;
  MPI_Init(&argc,&argv); //Initialize MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Current Process ID
  MPI_Comm_size(MPI_COMM_WORLD, &size); // Total Number of Processes
  MPI_Status status;
*/
  // -- Space discretization
  double xleft  = 0.;
  double xright = 1.;
  int        Nx = 200;         // Number of cells (ToBeRead)
  Nx = atoi(argv[1]);          // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 1
  cout << "Discretizing (" << xleft << "," << xright << ") with " << Nx << " cells" << endl;
  double     dx = abs(xright-xleft)/(double(Nx)); // regular grid space step
  // -- Trajectory
  double     T = 10.0;         // Final Time (TBR)
  T = atoi(argv[2]);           // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 2
  cout << "Final time: " << T << endl;
  int    K0 = 1;               // Number of shocks (TBR)  
  K0 = atoi(argv[3]);          // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 3
  cout << "Number of initial extremas: " << T << endl;
  // -- Source
  double sigma = 1.;           // Amplitude (TBR)
  sigma = atof(argv[4]);       // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 4
  cout << "White noise amplitude: " << sigma << endl;
  int    Kx = min(Nx,1);       // Number of Fourier modes (TBR)
  Kx = atoi(argv[5]);          // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 5
  cout << "Using " << Kx << " Fourier modes" << endl;
//   cout << "Using " << Mc << " Monte-Carlo realizations" << endl;
  int    Mc = 1;               // Monte-Carlo realizations -- or MPI nodes ??
//   Mc = atoi(argv[6]);          // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 6
//   cout << "Using " << Mc << " Monte-Carlo realizations" << endl;
  MTRand mt(1);                // Monte-Carlo Random number generator
  double u1 = mt();
  double u2 = mt();
  // -- Solver
  int  mysolver = 0;           // TBR ?
  double   cfl = 1./2;         // (solver dependent: TBR ?)
  const char * solvername;
  switch(mysolver){
    case 0:
      solvername = "exact";
    case 1:
      solvername = "Rusanov";
  }
  cout << "Using " << solvername << " Riemann solver with CFL = " << cfl << endl;
  // -- State initialization -------------------------------------------------------------
  double state[Nx+2][Mc];      // current state with 2 neighbouring ghost cells
  for(int i=0; i<Nx; i++){
    for(int m=0; m<Mc; m++){
      state[i+1][m]=sin(2.*K0*M_PI*(double(i)/Nx)); // TBR ?
    }
  }
  double fluxes[Nx+1][Mc];    // current flux (size equals number of interfaces)
  double cintfc[Nx+1][Mc];    // current flux parameter
  double fourier[Nx][Mc];     // current state Fourier coefficients
  double mass[Mc];        // current state "mass"
  double normel2[Mc];         // current state "energy"
  double normelinf[Mc];       // current state "max"
  // -- Time discretization
  double     t = 0.0;         // current time (initialized here)
  int     nsav = 0;
  double dtsav = .1;          // output memory saving time
  double tsav = dtsav;        // next saving time
  double dtmax = min(dtsav,cfl*dx); // maximal time step -- Rmk: here we know 1 = maximal speed at t=0 !!
  int       Nt = floor(T/dtmax); // minimal number of time steps not larger than dt
  int Nstepmax = max(int(1e9),Nt); // to stop the code in case of adaptive time step // stopping criterium
  double    dt = dtmax;
  int Nstep = 100;             // number of steps before copying to a file (balance between storage capacity and IO
  cout << "Saving snapshots on (" << t << "," << T << ") every " << dtsav << endl;
  // -- Boundary conditions
  for(int m=0; m<Mc; m++){
    state[0][m] = state[Nx][m];   // periodic
    state[Nx+1][m] = state[1][m]; // periodic
  }
  // -- Output
  string mystring; stringstream mystringstream;
  mystring.assign("res");
  mystringstream << Nx; mystring.append(mystringstream.str());
  mystringstream.str(string()); mystringstream.clear();
  mystring.append("T");
  mystringstream << T; mystring.append(mystringstream.str());
  mystringstream.str(string()); mystringstream.clear();
  mystring.append("K0");
  mystringstream << K0; mystring.append(mystringstream.str());
  mystringstream.str(string()); mystringstream.clear();
  mystring.append("sigma");
  mystringstream << sigma; mystring.append(mystringstream.str());
  mystringstream.str(string()); mystringstream.clear();
  mystring.append("Kx");
  mystringstream << Kx; mystring.append(mystringstream.str());
  mystringstream.str(string()); mystringstream.clear();
//   mystring.append("Mc");
//   mystringstream << Mc; mystring.append(mystringstream.str());
  string mystring1 = mystring;
  string mystring2 = mystring;
  mystring.append(".txt"); 
  cout << "Outputfilename is " << mystring << endl;
  ofstream outputfile; outputfile.open(mystring.c_str());
  mystring1.append(".log");
  cout << "Logfilename is " << mystring1 << endl;
  ofstream outputfile1; outputfile1.open(mystring1.c_str());
  mystring2.append(".dat"); 
  cout << "Postprocfilename is " << mystring2 << endl;
  ofstream outputfile2; outputfile2.open(mystring2.c_str());
  // Iterations: explicit FV scheme --------------------------------------------------
  int nstep = 0;
  int nreal = 0;
  while( (t<T) ) // & (nstep<Nstepmax) ) // stopping criterium
  {
    nstep++;
 //dt = dtmax; // upper-bound dtmax updated at every time step depending on updated dtsav
 dt = T-t; // other choice is useful only to force discrete times on a grid a priori
    // 1 -- Godunov approach: numerical fluxes at interfaces (explicit here)
    for(int m=0; m<Mc; m++){ // use MPI here ??
    for(int i=0; i<Nx+1; i++){
      switch(mysolver)
      {
	case 0: // Lax-Oleinik formula (explicit for Burgers flux) exact Riemann flux
	{
	  if (state[i+1][m]>state[i][m]) // relaxation wave
	  {
	    if (state[i+1][m]<=0)
	    {
	      fluxes[i][m] = FLUX(state[i+1][m]);
	      dt = min( dt, cfl*dx/abs(state[i+1][m]) );
	    }else{
	      if (state[i][m]>=0)
	      {
		fluxes[i][m] = FLUX(state[i][m]);
		dt = min( dt, cfl*dx/abs(state[i][m]) );
	      }else{ // state[i]<0 and state[i+1]>0
		fluxes[i][m] = 0.;
	      }
	    }
	  }else{ // shock wave
	    if (state[i][m]+state[i+1][m]<=0)
	    {
	      fluxes[i][m] = FLUX(state[i+1][m]);
	      dt = min( dt, cfl*dx/abs(state[i+1][m]) );
	    }else{
	      fluxes[i][m] = FLUX(state[i][m]);
	      dt = min( dt, cfl*dx/abs(state[i][m]) );
	    }
	  }
	  break;
	}
	case 1: // Rusanov approximate Riemann solver
	{
	  // parameter
	  cintfc[i][m] = max(abs(state[i][m]),abs(state[i+1][m]));
	  // conservative flux (right-going equal to left-going)
	  fluxes[i][m] = FLUX(state[i][m]) - cintfc[i][m]*( state[i+1][m] - state[i][m] );
	  dt = min( dt, cfl*dx/cintfc[i][m] );
	  break;
	}
      }
    }
    }
    // 2 -- CFL condition: new time step
   dt = min(dt,T-t); // maximal time step authorized // contrary to implementation above, should not necessarily be shared between all MC realizations !! (pathwise construction)
// if(dtsav>=1.e-4&&dt<1.e-4){ cout << "WARNING: fixed dt violates cfl" << endl; }
// dt = 1.e-4; // TO IMPOSE A FIXED TIME STEP : CHECK FOR EACH NEW TESTCASE !! <<<<<<<<<<<<<<<<<<<<<<<<<<
    // 3 -- Explicit FV: state update with numerical fluxes
    for(int m=0; m<Mc; m++){
    for(int i=0; i<Nx; i++){
      state[i+1][m] += (dt/dx)*fluxes[i][m];
      state[i+1][m] -= (dt/dx)*fluxes[i+1][m];
    }
    }
    // 4 -- Boundary conditions update
    for(int m=0; m<Mc; m++){
      state[0][m] = state[Nx][m];   // periodic
      state[Nx+1][m] = state[1][m]; // periodic
    }
    // 5 -- Source term (operator splitting)
    for(int m=0; m<Mc; m++){
    for(int i=0; i<=Nx+1; i++){
      for(int k=1; k<=Kx; k++){
        u1 = mt(); u2 = mt();
        state[i][m] += sqrt(dt/dx)*sqrt(-2.*log(u1))*cos(2.*M_PI*u2)*sigma*cos(2*k*M_PI*dx*i); //*sqrt(2/Kx)
        //u1 = mt(); u2 = mt();
        state[i][m] += sqrt(dt/dx)*sqrt(-2.*log(u1))*sin(2.*M_PI*u2)*sigma*sin(2*k*M_PI*dx*i); //*sqrt(2/Kx)
// state[i][m] += dt*sigma*cos(2*k*M_PI*dx*i); // ** a time-independent source term to test determinstic long-time convergence
      }
    }
    }
    // 6 -- Output: *** Fourier coeffs at successive times -- not necessarily on a grid if we use one single long Markov Chain !! (option) ***
    t += dt;
    for(int m=0; m<Mc; m++){
    for(int k=0; k<Nx/2; k++){
      for(int i=0; i<Nx; i++){
        fourier[k][m]=state[i][m]*cos(2.*k*M_PI*(double(i)/Nx))/sqrt(2.);
      }
      for(int i=0; i<Nx; i++){
        fourier[Nx/2+k][m]=state[i][m]*sin(2.*k*M_PI*(double(i)/Nx))/sqrt(2.);
      }
    }
    }    
// if(t>=tsav)
    {
    for(int m=0; m<Mc; m++){
      for(int i=1; i<Nx+1; i++){
	mass[m] += state[i][m];
	normel2[m] += state[i][m]*state[i][m];
	normelinf[m] = max(normelinf[m],state[i][m]);
      }
      normel2[m] = sqrt(normel2[m]);
    }
// nsav++; tsav += dtsav; // increment the counter for next output
    outputfile1 << dt << "  " << t << "  " << mass[0] << "  " << normel2[0] << "  " << normelinf[0] << "\n";
    for(int m=0; m<Mc; m++){
      for(int i=1; i<Nx+1; i++){ outputfile << state[i][m] << "  "; } outputfile << "\n"; // one period only: [0] and [Nx+1] excluded
//       for(int i=0; i<Nx; i++){ outputfile2 << fourier[i][m] << "  "; } outputfile2 << "\n"; // to construct fourier coefficients
      for(int i=0; i<Nx; i++){ outputfile2 << fourier[i][m] << "  "; } outputfile2 << "\n"; // to construct fourier coefficients
      cout << "Iteration " << nstep << " computed at time " << t << endl; //<< " and saved as " << nsav << endl;
    }
    }
  }
  outputfile.close();
  outputfile1.close();
  //return 0;
  return EXIT_SUCCESS;

}

