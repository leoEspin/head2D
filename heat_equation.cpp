/****************************************************************************
   CODE:        heat_equation.cpp
   AUTHOR:      Leonardo Espin
                lespines@umn.edu

   DESCRIPTION: ADI code for solving the 2D heat equation
        h_t = c^2 (h_xx + h_yy)
   subject to Dirichlet conditions in rectangle (0,2)X(0,1)
   Discretization parameters are read from text file 'input'

   Date written  01/11/2016
   revision date 01/14/2016  Wrote a test code for thomas algorithm
                 01/15/2016  Code validated checking for a steady state. code
                             was tested and debugged. All of the bugs were
                             caused by the index-shifting which occur when
                             inforcing the c++ index standard. because of
                             this it might be more efficient to write
                             subroutines for each individual task and pass
                             pointers as parameters. inside of each subr,
                             the task is easily done using the i=0;i<size
                             convention.
			     *Also I tried using parallelization, but it is
                              not trivial doing so, because of the pointers.
                              A quick search revealed that declaring dynamic
                              arrays inside of a parallel region is better
                              because openMP automatically takes care of
                              making copies of the variables. 

****************************************************************************/
#include <fstream>
#include <iostream>
#include <string>
#include <cmath> //used for sqrt in frobenius
using namespace std;

void getInput(int& nx,int& ny,double& dt,double& c,double& tfinal);
double** initialize(int,int);
void annihilate(double**, int);
void initialCond(double**,int,int,double);
void thomas(int size,double a,double b,double c,double* rhs);
void storeOutput(double** sol,double time, int indx,int indy);
double frobenius(double** matA,double** matB,int sx, int sy);

int main(){
  int nx,ny,tot_x,tot_y,i,j,counter=0;
  double t=0.,c,dt,tfinal,L=2.,P=1.;
  double dx,dy,ax,bx,cx,ay,by,cy;
  double** h;
  double** hf;

  getInput(nx,ny,dt,c,tfinal);
  dx=1./nx;
  dy=1./ny;
  ax=-c*c/(dx*dx);
  bx=2.*(1./dt - ax);
  cx=ax;
  ay=-c*c/(dy*dy);
  by=2.*(1./dt - ay);
  cy=ay;
  tot_x=static_cast<int>(nx*L)+2;//include extra boundary points
  tot_y=static_cast<int>(ny*P)+2;
  h = initialize(tot_x,tot_y);
  hf= initialize(tot_x,tot_y);
  initialCond(h,tot_x,tot_y,1.);
  initialCond(hf,tot_x,tot_y,0.);
  double* rhsx = new double [tot_x -2];//linear system only for interior points
  double* rhsy = new double [tot_y -2];
  while(t < tfinal){//solving system only for interior points
#pragma omp parallel for default(shared) private(rhsx) num_threads(4)
    for(j=1;j< tot_y -1;j++){
      for(i=1;i< tot_x -1;i++){ //j fixed. solve system in x direction
	rhsx[i-1]=-ay*h[i][j+1]+2.*(1./dt + ay)*h[i][j] -ay*h[i][j-1]; //rhsx index shifted by -1
      }
      rhsx[tot_x -3]=rhsx[tot_x -3]-ax;//size of rhsx is tot_x -2!!
      thomas(tot_x -2,ax,bx,cx,rhsx);
      for(i=1;i< tot_x -1;i++){ //j fixed.
	hf[i][j]=rhsx[i-1];     //rhsx index shifted by -1
      }
      hf[0][j]=0.;
      hf[tot_x -1][j]=1.;
    }
    // for(i=0;i< tot_x;i++){//is this BC loop necessary? not really (check linear system)
    //   hf[i][0]=0.;
    //   hf[i][tot_y -1]=1.;
    // }
#pragma omp parallel for default(shared) private(rhsy) num_threads(4)
    for(i=1;i< tot_x -1;i++){
      for(j=1;j< tot_y -1;j++){ //i fixed. solve system in y direction
	rhsy[j-1]=-ax*hf[i+1][j]+2.*(1./dt + ax)*hf[i][j] -ax*hf[i-1][j]; //rhsx index shifted by -1
      }
      rhsy[tot_y -3]=rhsy[tot_y -3]-ay;//size of rhsy is tot_y -2!!
      thomas(tot_y -2,ay,by,cy,rhsy);
      for(j=1;j< tot_y -1;j++){ //i fixed.
	h[i][j]=rhsy[j-1];     //rhsx index shifted by -1
      }
      h[i][0]=0.;
      h[i][tot_y -1]=1.;
    }
    for(j=0;j< tot_y;j++){
      h[0][j]=0.;
      h[tot_x -1][j]=1.;
    }
    t+=dt;
    counter++;
    if(counter%10000 == 0){
      cout << "time = " << t << ", two iterations differ by: "<< frobenius(h,hf,tot_x,tot_y)<<endl;
      //I'm assuming that h and hf converge to the same steady state. which
      //it is true. so hf is indeed h at t+1/2dt ?
    }
  }
  storeOutput(h,t,tot_x,tot_y);
  annihilate(h,tot_x);
  annihilate(hf,tot_x);
  delete [] rhsx;
  delete [] rhsy;
  return 0;
}

void getInput(int& nx,int& ny,double& dt,double& c,double& tfinal){
  ifstream file;
  string caca;
  file.open("input");
  file >> nx >> caca >> ny >> caca >> dt
       >> caca >> c >> caca >> tfinal;
  file.close();
}
void storeOutput(double** sol,double time, int indx,int indy){
  int i,j;
  ofstream file;
  file.open("solution.dat");
  for(i=0;i< indx;i++){
    for(j=0;j< indy;j++){
      file << sol[i][j] << ",";
    }
    file << endl;
  }
  file.close();
}
void initialCond(double** mat,int size_x,int size_y,double value){
  for(int i=0;i<size_x;i++){
    for(int j=0;j<size_y;j++){
      mat[i][j]=value;
    }
  }
}
double** initialize(int size_x,int size_y){//size_x rows, size_y columns
  double** mat;
  //cout << size_x << "|" << size_y << endl;
  mat = new double* [size_x];
  for(int i=0;i<size_x;i++){
    mat[i]=new double [size_y];
  }
  return mat;
}
void annihilate(double** mat, int size){//size = number of rows
  for(int i=0;i<size;i++){
    delete [] mat[i];
  }
  delete [] mat;
}
double frobenius(double** matA,double** matB,int sx, int sy){
  //summing only interior points
  double tmp=0.;
  for(int i=1;i<sx-1;i++){
    for(int j=1;j<sy-1;j++){
      tmp+=pow(matA[i][j] -matB[i][j],2);
    }
  }
  return sqrt(tmp);
}
void thomas(int size,double a,double b,double c,double* rhs){
  //solution is returned in rhs vector
  double* p = new double [size+1];
  double* q = new double [size+1];
  int i;
  double denom;
  p[0]=0.;//auxiliar component
  q[0]=0.;
  for(i=0;i< size;i++){
    denom=b -a*p[i];
    p[i+1]=c/denom;
    q[i+1]=(rhs[i] -a*q[i])/denom;
  }
  rhs[size-1]=q[size];
  for(i=size-2;i >=0;i--){
    rhs[i]=q[i+1] -p[i+1]*rhs[i+1];
  }
  delete [] p;
  delete [] q;
}
