#include "MultiTricubicSpline2.h"

template<int N> void 
MultiTricubicSpline<N>::Update()
{
  for (int n=0; n<N; n++) {
    // Do dF/dx
    UpdateX(0, 1, n);
    // Do dF/dy
    UpdateY(0, 2, n);
    // Do dF/dy
    UpdateZ(0, 3, n);
    // Do d2F/dxdy
    UpdateY(1, 4, n);
    // Do d2F/dxdz
    UpdateZ(1, 5, n);
    // Do d2F/dydz
    UpdateZ(2, 6, n);
    // Do d3F/dxdydz
    UpdateZ(4, 7, n);
  }
  UpToDate=true;
}

template<int N> void 
MultiTricubicSpline<N>::UpdateX(int source, int dest, int n)
{
  Grid &x = *Xgrid;
  Array<double,1> mu(Nx);
  
  // Loop over all y and z
  for (int iy=0; iy<Ny; iy++)
    for (int iz=0; iz<Nz; iz++) {
      // Set up tridiagonal set of equations
      // Initialize RHS of equations
      F(0,iy,iz)(dest,n) = 
	1.5*(F(1,iy,iz)(source,n)-F(0,iy,iz)(source,n))/(x(1)-x(0));
      F(Nx-1,iy,iz)(dest,n) = 
	1.5*(F(Nx-1,iy,iz)(source,n)-F(Nx-2,iy,iz)(source,n))
	/(x(Nx-1)-x(Nx-2));
      mu(0) = 0.5;

      // Solve tri-diagonal set of equations.  First eliminate lower
      // elements.
      for (int j=1; j<(Nx-1); j++) {
	double lambda = 0.5*(x(j+1)-x(j))/(x(j+1)-x(j-1));
	mu(j) = 0.5 - lambda;
	F(j,iy,iz)(dest,n) = 
	  3.0*(lambda*(F(j,iy,iz)(source,n)-
		       F(j-1,iy,iz)(source,n))/(x(j)-x(j-1))+ 
	       mu(j)*(F(j+1,iy,iz)(source,n)-F(j,iy,iz)(source,n))
	       /(x(j+1)-x(j)));
	double cj = 1.0 - lambda * mu(j-1);
	F(j,iy,iz)(dest,n) -= lambda * F(j-1,iy,iz)(dest,n);
	mu(j) /= cj;
	F(j,iy,iz)(dest,n) /= cj;
      }
      
      // Last element
      int j = Nx-1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      F(j,iy,iz)(dest,n) = 
	3.0*(lambda*(F(j,iy,iz)(source,n)-F(j-1,iy,iz)(source,n))/(x(j)-x(j-1)));
      double cj = 1.0 - lambda * mu(j-1);
      F(j,iy,iz)(dest,n) -= lambda * F(j-1,iy,iz)(dest,n);
      mu(j) /= cj;
      F(j,iy,iz)(dest,n) /= cj;
      
      // Now last d/dx is correct.  We proceed upward, back substituting.
      for (j=Nx-2; j>=0; j--) {
	F(j,iy,iz)(dest,n) -= mu(j) * F(j+1,iy,iz)(dest,n);
      }
    }      
}


template<int N> void 
MultiTricubicSpline<N>::UpdateY(int source, int dest, int n)
{
  Grid &y = *Ygrid;
  Array<double,1> mu(Ny);
  
  // Loop over all x and z
  for (int ix=0; ix<Nx; ix++)
    for (int iz=0; iz<Nz; iz++) {
      // Set up tridiagonal set of equations
      // Initialize RHS of equations
      F(ix,0,iz)(dest,n) = 1.5*(F(ix,1,iz)(source,n)-F(ix,0,iz)(source,n))/(y(1)-y(0));
      F(ix,Ny-1,iz)(dest,n) = 1.5*(F(ix,Ny-1,iz)(source,n)-F(ix,Ny-2,iz)(source,n))
	/(y(Ny-1)-y(Ny-2));
      mu(0) = 0.5;

      // Solve tri-diagonal set of equations.  First eliminate lower
      // elements.
      for (int j=1; j<(Ny-1); j++) {
	double lambda = 0.5*(y(j+1)-y(j))/(y(j+1)-y(j-1));
	mu(j) = 0.5 - lambda;
	F(ix,j,iz)(dest,n) = 
	  3.0*(lambda*(F(ix,j,iz)(source,n)-F(ix,j-1,iz)(source,n))/(y(j)-y(j-1))+ 
	       mu(j) *(F(ix,j+1,iz)(source,n)-F(ix,j,iz)(source,n))/(y(j+1)-y(j)));
	double cj = 1.0 - lambda * mu(j-1);
	F(ix,j,iz)(dest,n) -= lambda * F(ix,j-1,iz)(dest,n);
	mu(j) /= cj;
	F(ix,j,iz)(dest,n) /= cj;
      }
      
      // Last element
      int j = Ny-1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      F(ix,j,iz)(dest,n) = 
	3.0*(lambda*(F(ix,j,iz)(source,n)-F(ix,j-1,iz)(source,n))/(y(j)-y(j-1)));
      double cj = 1.0 - lambda * mu(j-1);
      F(ix,j,iz)(dest,n) -= lambda * F(ix,j-1,iz)(dest,n);
      mu(j) /= cj;
      F(ix,j,iz)(dest,n) /= cj;
      
      // Now last d/dx is correct.  We proceed upward, back substituting.
      for (j=Ny-2; j>=0; j--) {
	F(ix,j,iz)(dest,n) -= mu(j) * F(ix,j+1,iz)(dest,n);
      }
    }      
}


template<int N> void 
MultiTricubicSpline<N>::UpdateZ(int source, int dest, int n)
{
  Grid &z = *Zgrid;
  Array<double,1> mu(Nz);
  
  // Loop over all x and y
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++) {
      // Set up tridiagonal set of equations
      // Initialize RHS of equations
      F(ix,iy,0)(dest,n) = 1.5*(F(ix,iy,1)(source,n)-F(ix,iy,0)(source,n))/(z(1)-z(0));
      F(ix,iy,Nz-1)(dest,n) = 1.5*(F(ix,iy,Nz-1)(source,n)-F(ix,iy,Nz-2)(source,n))
	/(z(Nz-1)-z(Nz-2));
      mu(0) = 0.5;

      // Solve tri-diagonal set of equations.  First eliminate lower
      // elements.
      for (int j=1; j<(Nz-1); j++) {
	double lambda = 0.5*(z(j+1)-z(j))/(z(j+1)-z(j-1));
	mu(j) = 0.5 - lambda;
	F(ix,iy,j)(dest,n) = 
	  3.0*(lambda*(F(ix,iy,j)(source,n)-F(ix,iy,j-1)(source,n))/(z(j)-z(j-1))+ 
	       mu(j) *(F(ix,iy,j+1)(source,n)-F(ix,iy,j)(source,n))/(z(j+1)-z(j)));
	double cj = 1.0 - lambda * mu(j-1);
	F(ix,iy,j)(dest,n) -= lambda * F(ix,iy,j-1)(dest,n);
	mu(j) /= cj;
	F(ix,iy,j)(dest,n) /= cj;
      }
      
      // Last element
      int j = Nz-1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      F(ix,iy,j)(dest,n) = 
	3.0*(lambda*(F(ix,iy,j)(source,n)-F(ix,iy,j-1)(source,n))/(z(j)-z(j-1)));
      double cj = 1.0 - lambda * mu(j-1);
      F(ix,iy,j)(dest,n) -= lambda * F(ix,iy,j-1)(dest,n);
      mu(j) /= cj;
      F(ix,iy,j)(dest,n) /= cj;
      
      // Now last d/dx is correct.  We proceed upward, back substituting.
      for (j=Nz-2; j>=0; j--) {
	F(ix,iy,j)(dest,n) -= mu(j) * F(ix,iy,j+1)(dest,n);
      }
    }      
}


