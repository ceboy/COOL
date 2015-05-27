#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include <algorithm> // Pour sort()
#include <vector>  // pour la classe vector
#include <ctime>  // Calcul du time d execution
#include "mtrand.cpp"
#include <mpi.h>

using namespace std;

#ifndef M_PI
   #define M_PI 3.14159265358979323846
#endif

#define MPI_MASTER 0

int main(int argc, char** argv)
{

  int rank, size;
  MPI_Init(&argc,&argv); //Initialize MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Current Process ID
  MPI_Comm_size(MPI_COMM_WORLD, &size); // Total Number of Processes
  MPI_Status status;
  int tag = 11 ;
 // int tag2= 22 ;
  //int nb_realisation = 0;
  
  double time_initial= clock();
  double time_final  = 0.0 ;
  double time_CPU    = 0.0 ;
//   double t1=0.0 ;
//   double t2=0.0 ;
//   double t_rea=0.0;
  string sep(38, '=');
  // Jeu de donnees diverses
  int i = 0;
  int j = 0;
  int k = 0;

  // Generateur de nombre aleatoire
  MTRand mt(1);

  double u1 = mt();
  double u2 = mt();
  double sigma = 0.1;

     int Effect=0;
  
 // Discretisation d'espace
  double left  = 0.;
  double right = 1.;
  int       Nx = 500;      // on prend Nx points
  int        l = Nx/10;    // moitié de nombre de mode de Fourier
  
 // Discretisation en time
  double dx = abs(right-left)/(Nx);
  double cfl = 0.5;
  //double c=cfl;
  double t = 0.05   ;     // instant final
  double dt = cfl*dx; // pas de time
  int Nt = floor(t/dt); // nb d'instants
  vector <double> v1; // realisations de l'instant final
  vector <double> v2; // realisations de l'instant final
  vector <double> v3; // realisations de l'instant final
  double u_norm;
  //double t_temp=0.;

  
  // Generation des realisations
  int Nreal= 1024;
  int rea = 0;
  
  // Fichiers de sortie
  char file1[100]; sprintf(file1, "BurgerSolBruit_Godunov_MPI-T50.dat");         ofstream parBR(file1);
  char file2[100]; sprintf(file2, "BurgerSolSansBruit_Godunov_MPI-T50.dat");     ofstream parSB(file2);
  char file3[100]; sprintf(file3, "BurgerSolMoyen_Godunov_MPI-T50.dat");          ofstream parMY(file3);
  char file11[100];sprintf(file11,"BurgerVariance_MPIT50.dat");                  ofstream parVR(file11);
  //####################Fichier pour R
  char file4[100]; sprintf(file4, "Burger_Godunov-Point1_MPI-T50.txt");          ofstream parFip1(file4);
  char file5[100]; sprintf(file5, "Burger_Godunov-Point2_MPI-T50.txt");          ofstream parFip2(file5);
  char file6[100]; sprintf(file6, "Burger_Godunov-Point3_MPI-T50.txt");           ofstream parFip3(file6);

  char file8[100]; sprintf(file8, "Frepar_Godunov_MPI-1_T50.dat");              ofstream parFP1(file8);
  char file9[100]; sprintf(file9, "Frepar_Godunov_MPI-2_T50.dat");              ofstream parFP2(file9);
  char file7[100]; sprintf(file7, "Frepar_Godunov_MPI-3_T50.dat");              ofstream parFP3(file7);
  cout << " Allocation memoire..." << endl;
  double ** w = new double* [Nt];
  double ** U = new double* [Nt+1];
  double ** c = new double* [Nt+1];
  double ** F = new double* [Nt+1];
  double ** U_GD = new double* [Nt+1];
  double ** c_GD = new double* [Nt+1];    // Fluxes numériques
  double ** F_GD = new double* [Nt+1];    // Fluxes numériques
  
  double *  U_moy  = new double [Nx];
  double *  Vari_U = new double [Nx];   // Sauvegarder les variance pour tous les point d'espace
  double *  valtmp = new double [2*l+2];
  double *  Stock1 = new double[Nreal];
  double *  Stock2 = new double[Nreal];
  double *  Stock3 = new double[Nreal];

   for (i=0; i<Nt; i++)
  {
      w[i] = new double [Nx];  
  }
      

  for (i=0; i<Nt+1; i++)
  {
      U[i]  =  new double [Nx+2];
      c[i]  =  new double [Nx+2];
      F[i]  =  new double [Nx+2];
    U_GD[i] =  new double [Nx+2];  // Solution numerique d'équation transport
    c_GD[i] =  new double [Nx+2];          // Maybe here you don't need that much
    F_GD[i] =  new double [Nx+2];          // Maybe here you don't need that much
   }
  
  
    cout << " Nombre de pas de time : " << Nt << endl;


//Condition initial

  for (j=1; j<Nx+1; j++)
    {   
      U[0][j] =// (sin(2.*M_PI*(left+(j)*dx))-sin(2.*M_PI*(left+(j-1)*dx)))/(2.*M_PI*dx);
             (cos(2.*M_PI*(left+(j-1)*dx))-cos(2.*M_PI*(left+j*dx)))/(2.*M_PI*dx);
      U_GD[0][j]  = U[0][j];  
    }
    
     U[0][0] = U[0][Nx];
     U[0][Nx+1]=U[0][1];
    
     U_GD[0][0] =U_GD[0][Nx];
     U_GD[0][Nx+1]=U_GD[0][1];
    
    
// Solution sans bruit
  for (i=1; i<Nt+1; i++)
    {  
      
      for (j=0; j<Nx+1; j++)
       { 
	 c_GD[i-1][j]=0.5*(U_GD[i-1][j+1]+U_GD[i-1][j]);
	 
	 if (U_GD[i-1][j]<=U_GD[i-1][j+1]) {
	    if (U_GD[i-1][j]*U_GD[i-1][j+1]<=0) {
	       F_GD[i-1][j]=  0;
	    }
	    else {
	      if (U_GD[i-1][j+1]<0)
	       F_GD[i-1][j]=  pow(U_GD[i-1][j+1],2)*0.5;
	      else if (U_GD[i-1][j]>0)
	       F_GD[i-1][j]=  pow(U_GD[i-1][j],2)*0.5;
	    }
	  }
	  
	 else {
	   if (c_GD[i-1][j]>=0) 
	       F_GD[i-1][j]=  pow(U_GD[i-1][j],2)*0.5;
	   else
	       F_GD[i-1][j]=  pow(U_GD[i-1][j+1],2)*0.5;
	  } 	 
       }
       
      for (j=1; j<Nx+1; j++)
       {
        U_GD[i][j]  =  U_GD[i-1][j] - (dt/dx)*(F_GD[i-1][j]-F_GD[i-1][j-1]); // le Schéma
       }
      
      U_GD[i][0]    = U_GD[i][Nx];               // conditions aux bords périodiques
      U_GD[i][Nx+1] = U_GD[i][1];
      
    }
  
//########################Solution avec bruit############################################
for ( rea=rank*Nreal/size; rea<(rank+1)*Nreal/size; rea++ )
{
  
   // t1 = clock();  
  
    cout << "This is the REALISATION No. " << (rea+1) << " From Processor "<<rank<<endl;
   //#######################Engendrer du bruit####################################
    mt.seed(time(NULL)+rank+rea);
    for (i=0;i<Nt;i++)
     {
       
	
      	 for (k=1; k<=2*l+2; k++)
           {
             u1 = mt();                               // Generation VA de loi U(0,1)
             u2 = mt();                               // Generation VA de loi U(0,1)
            valtmp[k-1] = sqrt(-2.*log(u1))*cos(2.*M_PI*u2);   //engendre des V-A qui suivent une loi Gaussianne
    	  }
// 	    
	     for (j=0;j<Nx;j++){
	       	w[i][j]=0;	                                 
 	      for (  k=1; k<=l+1;  k++)
                { w[i][j] =w[i][j] + valtmp[k-1]*sin(2.*M_PI*(k-0.5*l-1)*(j+1.0-0.5)*dx);}
              for (k=l+2; k<=2*l+2; k++)
	        { w[i][j] =w[i][j] + valtmp[k-1]*cos(2.*M_PI*(k-l-0.5*l-2)*(j+1.0-0.5)*dx);}
 	     }
     }
    
 // ######################RESOLUTION#########################

  //##################################   
   for (i=1; i<Nt+1; i++)
    {       
          for (j=0; j<Nx+1; j++)
             { 
        	 c[i-1][j]=0.5*(U[i-1][j+1]+U[i-1][j]);
	 
	         if (U[i-1][j]<=U[i-1][j+1]) {
	           if (U[i-1][j]*U[i-1][j+1]<=0) {
	                F[i-1][j]=  0;
	            }
	           else {
	               if (U[i-1][j+1]<0)
	                F[i-1][j]=  pow(U[i-1][j+1],2)*0.5;
	               else if (U[i-1][j]>0)
	                F[i-1][j]=  pow(U[i-1][j],2)*0.5;
	             }
	           }
	  
	       else {
	           if (c[i-1][j]>=0) 
	              F[i-1][j]=  pow(U[i-1][j],2)*0.5;
	          else
	              F[i-1][j]=  pow(U[i-1][j+1],2)*0.5;
	         } 	 
             }//Calculer des fluxes numériques
      //##################################   
      for (j=1; j<Nx+1; j++)
       {
        U[i][j]  =  U[i-1][j] - (dt/dx)*(F[i-1][j]-F[i-1][j-1])+sigma*sqrt(dt)*w[i-1][j-1]; // Schema
       }
      
      U[i][0]    = U[i][Nx];
      U[i][Nx+1] = U[i][1];      
   }
//###################FIN De RÉSOLUTION##########################

   MPI_Barrier(MPI_COMM_WORLD);
    if (rank == MPI_MASTER){   
           u_norm=0.0;
             for (j=1;j<Nx+1;j++){
	      u_norm+=pow(U[Nt][j]-U_GD[Nt][j],2);
	    }
      if (pow(u_norm,0.5)<1000)                         // Assurer que la solution n'avait pas explosé
         {
           parFip1<< U[Nt][111]-U_GD[Nt][111]<<endl;
           parFip2<< U[Nt][289]-U_GD[Nt][289]<<endl; 
           parFip3<< U[Nt][426]-U_GD[Nt][426]<<endl; 
     
	   Stock1[Effect]=U[Nt][111]-U_GD[Nt][111];
           Stock2[Effect]=U[Nt][289]-U_GD[Nt][289];
           Stock3[Effect]=U[Nt][426]-U_GD[Nt][426];
	    Effect=Effect+1;  
      
            for (j=1; j<Nx+1; j++)
           {
	      U_moy[j-1] += U[Nt][j];            //ici, c'est plutôt la somme, on calcule le vrai moyen à la ligne #291
	      Vari_U[j-1] += pow(U[Nt][j],2);
           }
          }
         for(int idproc=1;idproc<size;idproc++)
	 {
	  MPI_Recv(&(U[Nt][0]),Nx,MPI_DOUBLE,idproc,tag,MPI_COMM_WORLD,&status);//  processor 0 recevoit U[Nt] d'autres processsors 
	             u_norm=0.0;
                  for (j=1;j<Nx+1;j++){
		     u_norm+=pow(U[Nt][j]-U_GD[Nt][j],2);
		   }
	  if (pow(u_norm,0.5)<1000)
            {
             parFip1<< U[Nt][111]-U_GD[Nt][111]<<endl;
             parFip2<< U[Nt][289]-U_GD[Nt][289]<<endl; 
             parFip3<< U[Nt][426]-U_GD[Nt][426]<<endl; 
	  
	     Stock1[Effect]=U[Nt][111]-U_GD[Nt][111];
             Stock2[Effect]=U[Nt][289]-U_GD[Nt][289];
             Stock3[Effect]=U[Nt][426]-U_GD[Nt][426];
             Effect=Effect+1;
	     for (j=1; j<Nx+1; j++)
	      {
	        U_moy[j-1] += U[Nt][j];              
	       Vari_U[j-1] += pow(U[Nt][j],2);	     
             }
	    } 
	   }
    }//End for the "rank == MPI_MASTER"
     else{
         MPI_Send(&(U[Nt][0]),Nx, MPI_DOUBLE, MPI_MASTER, tag, MPI_COMM_WORLD);   // d'autres processsors envoient U[Nt] au processor 0    
     }  
  
//     t2 = clock();
//     t_rea=(t2-t1)/CLOCKS_PER_SEC ;
//     
//     {
//     cout << "+" << sep << "+" <<endl;
//     cout << "|";
//     cout << setw(18) << " Temps d'execution :";
//     cout << setw(3) << setprecision(0) << floor(t_rea/3600.) << " H ";
//     cout << setw(2) << floor(( t_rea - floor(t_rea/3600.)*3600. )/60.) << " MIN ";
//     cout << setw(2) << floor(t_rea - floor(( t_rea - floor(t_rea/3600.)*3600. )/60.)*60.) << " S |"<<endl;
//     cout << "|";
//     cout << std::left << setw(19) << " Temps en S" << ":";
//     cout << std::right << setw(15) << t_rea << " S |" << endl;
//     cout << "+" << sep << "+" <<endl;
//     }
}//End for one realisation,voir ligne 182

cout<<"Effect is "<<Effect<<" on processor "<<rank<<endl;

//##################Fonction répartition##########################
  if(rank==MPI_MASTER){
   for (i=0; i<Effect; i++)
   {
     v1.push_back(Stock1[i]);
     v2.push_back(Stock2[i]);
     v3.push_back(Stock3[i]);
   }
  
    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());
    sort(v3.begin(), v3.end());

 // Fonction de repartition au time final
   for (i=0; i<Effect; i++)
    {
     parFP1 << v1[i] << " " << (double)i/Effect << endl;
     parFP2 << v2[i] << " " << (double)i/Effect << endl;
     parFP3 << v3[i] << " " << (double)i/Effect << endl;
    }
    
    
      for (j=1; j<Nx+1; j++)
      {
	Vari_U[j-1] = pow(Effect,2)*Vari_U[j-1]-Effect*pow(U_moy[j-1],2);
      }
     for (j=1;j<Nx+1;j++) 
      {
        U_moy[j-1]  = U_moy[j-1]/double(Effect);
        Vari_U[j-1] = Vari_U[j-1]/double((Effect-1)*pow(Effect,2));	
      }
  
  //##################Fonction répartition##########################

// ######################RÉSOLUTION#########################   

    
        for (j=1; j<Nx+1; j++)
   {
     parBR << left + (j-1)*dx+dx/2. << " " << U[Nt][j] << endl;
     parSB << left + (j-1)*dx+dx/2. << " " << U_GD[Nt][j] << endl;
     parMY << left + (j-1)*dx+dx/2. << " " << U_moy[j-1] << endl;
     parVR << left + (j-1)*dx+dx/2. << " " << Vari_U[j-1] << endl;     
    }
  }  
    
    
    
    

  
  
 // Cleaning 
  parBR.close();
  parSB.close();
  parFip1.close();
  parMY.close();
  parVR.close();
  parFip2.close();
  parFip3.close();

  parFP1.close();
  parFP2.close();
  parFP3.close();
  
  v1.clear();
  v2.clear();
  v3.clear();
  
 for (i=0; i<Nt+1; i++)
  {
    delete [] U[i];
    delete [] U_GD[i];
  }
    for (i=0; i<Nt; i++)
  {
    delete [] w[i];
  }


  delete [] w;
  delete [] U;
  delete [] U_GD;
  delete valtmp;
  delete U_moy;
  delete Stock1;
  delete Stock2;
  delete Stock3;
  
  
  // Affichage du time d'execution
  time_final = clock ();
  time_CPU = (time_final - time_initial) / CLOCKS_PER_SEC ; // secondes
  cout << "+" << sep << "+" <<endl;
  cout << "|";
  cout << setw(18) << " Temps d'execution :";
  cout << setw(3) << setprecision(0) << floor(time_CPU/3600.) << " H ";
  cout << setw(2) << floor(( time_CPU - floor(time_CPU/3600.)*3600. )/60.) << " MIN ";
  cout << setw(2) << floor(time_CPU - floor(( time_CPU - floor(time_CPU/3600.)*3600. )/60.)*60.) << " S |"<<endl;
  cout << "|";
  cout << std::left << setw(19) << " Temps en S" << ":";
  cout << std::right << setw(15) << time_CPU << " S |" << endl;
  cout << "+" << sep << "+" <<endl;
  
  MPI_Finalize();
  return 0;
}