#include "ComplexMultiTricubicSpline.h"

void 
ComplexMultiTricubicSpline::Update()
{
  if (!Periodic)
    for (int i=0; i<N; i++) {
      ///////////////         ////////////////////
      // Real part //         // Imaginary part //
      ///////////////         ////////////////////
      // Do dF/dx
      UpdateX(0, 1, i);       UpdateX(8,   9, i);
      // Do dF/dy
      UpdateY(0, 2, i);       UpdateY(8,  10, i);
      // Do dF/dy
      UpdateZ(0, 3, i);       UpdateZ(8,  11, i);
      // Do d2F/dxdy
      UpdateY(1, 4, i);       UpdateY(9,  12, i);
      // Do d2F/dxdz
      UpdateZ(1, 5, i);       UpdateZ(9,  13, i);
      // Do d2F/dydz
      UpdateZ(2, 6, i);       UpdateZ(10, 14, i);
      // Do d3F/dxdydz
      UpdateZ(4, 7, i);       UpdateZ(12, 15, i);
    }
  else
    for (int i=0; i<N; i++) {
      ///////////////                ////////////////////
      // Real part //                // Imaginary part //
      ///////////////                ////////////////////
      // Do dF/dx
      UpdateXPeriodic(0, 1, i);      UpdateXPeriodic(8,   9, i);
      // Do dF/dy
      UpdateYPeriodic(0, 2, i);      UpdateYPeriodic(8,  10, i);
      // Do dF/dy
      UpdateZPeriodic(0, 3, i);      UpdateZPeriodic(8,  11, i);
      // Do d2F/dxdy
      UpdateYPeriodic(1, 4, i);      UpdateYPeriodic(9,  12, i);
      // Do d2F/dxdz
      UpdateZPeriodic(1, 5, i);      UpdateZPeriodic(9,  13, i);
      // Do d2F/dydz
      UpdateZPeriodic(2, 6, i);      UpdateZPeriodic(10, 14, i);
      // Do d3F/dxdydz
      UpdateZPeriodic(4, 7, i);      UpdateZPeriodic(12, 15, i);
    }
  UpToDate=true;
}


void 
ComplexMultiTricubicSpline::UpdateX(int source, int dest, int i)
{
  Grid &x = *Xgrid;
  Array<double,1> mu(Nx);
  ///////////////
  // Real part //
  ///////////////
  // Loop over all y and z
  for (int iy=0; iy<Ny; iy++)
    for (int iz=0; iz<Nz; iz++) {
      // Set up tridiagonal set of equations
      // Initialize RHS of equations
      F(0,iy,iz,i)[dest] = 1.5*(F(1,iy,iz,i)[source]-F(0,iy,iz,i)[source])/(x(1)-x(0));
      F(Nx-1,iy,iz,i)[dest] = 1.5*(F(Nx-1,iy,iz,i)[source]-F(Nx-2,iy,iz,i)[source])
	/(x(Nx-1)-x(Nx-2));
      mu(0) = 0.5;

      // Solve tri-diagonal set of equations.  First eliminate lower
      // elements.
      for (int j=1; j<(Nx-1); j++) {
	double lambda = 0.5*(x(j+1)-x(j))/(x(j+1)-x(j-1));
	mu(j) = 0.5 - lambda;
	F(j,iy,iz,i)[dest] = 
	  3.0*(lambda*(F(j,iy,iz,i)[source]-F(j-1,iy,iz,i)[source])/(x(j)-x(j-1))+ 
	       mu(j) *(F(j+1,iy,iz,i)[source]-F(j,iy,iz,i)[source])/(x(j+1)-x(j)));
	double cj = 1.0 - lambda * mu(j-1);
	F(j,iy,iz,i)[dest] -= lambda * F(j-1,iy,iz,i)[dest];
	mu(j) /= cj;
	F(j,iy,iz,i)[dest] /= cj;
      }
      
      // Last element
      int j = Nx-1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      F(j,iy,iz,i)[dest] = 
	3.0*(lambda*(F(j,iy,iz,i)[source]-F(j-1,iy,iz,i)[source])/(x(j)-x(j-1)));
      double cj = 1.0 - lambda * mu(j-1);
      F(j,iy,iz,i)[dest] -= lambda * F(j-1,iy,iz,i)[dest];
      mu(j) /= cj;
      F(j,iy,iz,i)[dest] /= cj;
      
      // Now last d/dx is correct.  We proceed upward, back substituting.
      for (j=Nx-2; j>=0; j--) {
	F(j,iy,iz,i)[dest] -= mu(j) * F(j+1,iy,iz,i)[dest];
      }
    }
}


void
ComplexMultiTricubicSpline::UpdateXPeriodic(int source, int dest, int i)
{
  Grid &x = *Xgrid;
  int M = x.NumPoints - 1;
  // Check to make sure first and last point are the same.
  Array<double,1> lambda(M), mu(M), d(M), gamma(M);

  ///////////////
  // Real part //
  ///////////////
  for (int iy=0; iy<Ny; iy++) 
    for (int iz=0; iz<Nz; iz++) {
      assert (F(0,iy,iz,i)[source] == F(M,iy,iz,i)[source]);
      // Setup lambdas, mus, and d's
      lambda(0) = 0.5*(x(1)-x(0))/(x(M)-x(M-2));
      mu(0)     = 0.5-lambda(0);
      d(0)      = 3.0*(lambda(0)*(F(0,iy,iz,i)[source]-F(M-1,iy,iz,i)[source])/(x(M)-x(M-1)) +
		       mu(0)    *(F(1,iy,iz,i)[source]-F(0,iy,iz,i)[source])  /(x(1)-x(0)));
      lambda(M-1) = 0.5*(x(M)-x(M-1))/(x(M)-x(M-2));
      mu(M-1)     = 0.5*lambda(M-1);
      d(M-1)      = 3.0*(lambda(M-1)*(F(M-1,iy,iz,i)[source]-F(M-2,iy,iz,i)[source])/(x(M-1)-x(M-2)) +
			 mu(M-1)    *(F(M,iy,iz,i)[source]  -F(M-1,iy,iz,i)[source])/(x(M)  -x(M-1)));
      gamma(M-1)  = 1.0;
      gamma(0) = lambda(0);
      
      for (int ix=1; ix<M; ix++) {
	lambda(ix) = 0.5*(x(ix+1)-x(ix))/(x(ix+1)-x(ix-1));
	mu(ix)     = 0.5-lambda(ix);
	d(ix)      = 3.0*(lambda(ix)*(F(ix,iy,iz,i)[source]-F(ix-1,iy,iz,i)[source])/(x(ix)-x(ix-1)) +
			 mu(ix)   *(F(ix+1,iy,iz,i)[source]-F(ix,iy,iz,i)[source])/(x(ix+1)-x(ix)));
      }
      // Solve down lower triangular part
      for (int ix=1; ix<(M-1); ix++) {
	double diag = 1.0-lambda(ix)*mu(ix-1);
	gamma(ix) = -gamma(ix-1)*lambda(ix);
	d(ix) -= lambda(ix) * d(ix-1);
	double diagInv = 1.0/diag;
	gamma(ix) *= diagInv;
	d(ix)     *= diagInv;
	mu(ix)    *= diagInv;
	// last row
	d(M-1) -= mu(M-1) * d(ix-1);
	gamma(M-1) -= mu(M-1) * gamma(ix-1);
	mu(M-1) = -mu(M-1)*mu(ix-1);
      }
      // last row
      // mu is really on top of lambda in the last row
      lambda(M-1) += mu(M-1);
      d(M-1) -= lambda(M-1) * d(M-2);
      gamma(M-1) -= lambda(M-1) * (gamma(M-2)+ mu(M-2));
      // Compute last derivative;
      F(M-1,iy,iz,i)[dest] = d(M-1)/gamma(M-1);
      
      // Now proceed up upper diagonal, backsubstituting
      for (int ix=M-2; ix>=0; ix--)
	F(ix,iy,iz,i)[dest] = d(ix) - mu(ix)*F(ix+1,iy,iz,i)[dest] - gamma(ix)*F(M-1,iy,iz,i)[dest];
      
      // Finally, assign repeated last element for PBC
      F(M,iy,iz,i)[dest] = F(0,iy,iz,i)[dest];
    }
}





void 
ComplexMultiTricubicSpline::UpdateY(int source, int dest, int i)
{
  Grid &y = *Ygrid;
  Array<double,1> mu(Ny);

  ///////////////
  // Real part //
  /////////////// 
  // Loop over all x and z
  for (int ix=0; ix<Nx; ix++)
    for (int iz=0; iz<Nz; iz++) {
      // Set up tridiagonal set of equations
      // Initialize RHS of equations
      F(ix,0,iz,i)[dest] = 1.5*(F(ix,1,iz,i)[source]-F(ix,0,iz,i)[source])/(y(1)-y(0));
      F(ix,Ny-1,iz,i)[dest] = 1.5*(F(ix,Ny-1,iz,i)[source]-F(ix,Ny-2,iz,i)[source])
	/(y(Ny-1)-y(Ny-2));
      mu(0) = 0.5;

      // Solve tri-diagonal set of equations.  First eliminate lower
      // elements.
      for (int j=1; j<(Ny-1); j++) {
	double lambda = 0.5*(y(j+1)-y(j))/(y(j+1)-y(j-1));
	mu(j) = 0.5 - lambda;
	F(ix,j,iz,i)[dest] = 
	  3.0*(lambda*(F(ix,j,iz,i)[source]-F(ix,j-1,iz,i)[source])/(y(j)-y(j-1))+ 
	       mu(j) *(F(ix,j+1,iz,i)[source]-F(ix,j,iz,i)[source])/(y(j+1)-y(j)));
	double cj = 1.0 - lambda * mu(j-1);
	F(ix,j,iz,i)[dest] -= lambda * F(ix,j-1,iz,i)[dest];
	mu(j) /= cj;
	F(ix,j,iz,i)[dest] /= cj;
      }
      
      // Last element
      int j = Ny-1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      F(ix,j,iz,i)[dest] = 
	3.0*(lambda*(F(ix,j,iz,i)[source]-F(ix,j-1,iz,i)[source])/(y(j)-y(j-1)));
      double cj = 1.0 - lambda * mu(j-1);
      F(ix,j,iz,i)[dest] -= lambda * F(ix,j-1,iz,i)[dest];
      mu(j) /= cj;
      F(ix,j,iz,i)[dest] /= cj;
      
      // Now last d/dx is correct.  We proceed upward, back substituting.
      for (j=Ny-2; j>=0; j--) {
	F(ix,j,iz,i)[dest] -= mu(j) * F(ix,j+1,iz,i)[dest];
      }
    }
}



void
ComplexMultiTricubicSpline::UpdateYPeriodic(int source, int dest, int i)
{
  Grid &y = *Ygrid;
  int M = y.NumPoints - 1;
  // Check to make sure first and last point are the same.
  Array<double,1> lambda(M), mu(M), d(M), gamma(M);

  ///////////////
  // Real part //
  ///////////////  
  for (int ix=0; ix<Nx; ix++) 
    for (int iz=0; iz<Nz; iz++) {
      assert (F(ix,0,iz,i)[source] == F(ix,M,iz,i)[source]);
  
      // Setup lambdas, mus, and d's
      lambda(0) = 0.5*(y(1)-y(0))/(y(M)-y(M-2));
      mu(0)     = 0.5-lambda(0);
      d(0)      = 3.0*(lambda(0)*(F(ix,0,iz,i)[source]-F(ix,M-1,iz,i)[source])/(y(M)-y(M-1)) +
		       mu(0)    *(F(ix,1,iz,i)[source]-F(ix,0,iz,i)[source])  /(y(1)-y(0)));
      lambda(M-1) = 0.5*(y(M)-y(M-1))/(y(M)-y(M-2));
      mu(M-1)     = 0.5*lambda(M-1);
      d(M-1)      = 3.0*(lambda(M-1)*(F(ix,M-1,iz,i)[source]-F(ix,M-2,iz,i)[source])/(y(M-1)-y(M-2)) +
			 mu(M-1)    *(F(ix,M,iz,i)[source]  -F(ix,M-1,iz,i)[source])/(y(M)  -y(M-1)));
      gamma(M-1)  = 1.0;
      gamma(0) = lambda(0);
      
      for (int iy=1; iy<M; iy++) {
	lambda(iy) = 0.5*(y(iy+1)-y(iy))/(y(iy+1)-y(iy-1));
	mu(iy)     = 0.5-lambda(iy);
	d(iy)      = 3.0*(lambda(iy)*(F(ix,iy,iz,i)[source]-F(ix,iy-1,iz,i)[source])/(y(iy)-y(iy-1)) +
			 mu(iy)   *(F(ix,iy+1,iz,i)[source]-F(ix,iy,iz,i)[source])/(y(iy+1)-y(iy)));
      }
      // Solve down lower triangular part
      for (int iy=1; iy<(M-1); iy++) {
	double diag = 1.0-lambda(iy)*mu(iy-1);
	gamma(iy) = -gamma(iy-1)*lambda(iy);
	d(iy) -= lambda(iy) * d(iy-1);
	double diagInv = 1.0/diag;
	gamma(iy) *= diagInv;
	d(iy)     *= diagInv;
	mu(iy)    *= diagInv;
	// last row
	d(M-1) -= mu(M-1) * d(iy-1);
	gamma(M-1) -= mu(M-1) * gamma(iy-1);
	mu(M-1) = -mu(M-1)*mu(iy-1);
      }
      // last row
      // mu is really on top of lambda in the last row
      lambda(M-1) += mu(M-1);
      d(M-1) -= lambda(M-1) * d(M-2);
      gamma(M-1) -= lambda(M-1) * (gamma(M-2)+ mu(M-2));
      // Compute last derivative;
      F(ix,M-1,iz,i)[dest] = d(M-1)/gamma(M-1);
      
      // Now proceed up upper diagonal, backsubstituting
      for (int iy=M-2; iy>=0; iy--)
	F(ix,iy,iz,i)[dest] = d(iy) - mu(iy)*F(ix,iy+1,iz,i)[dest] - gamma(iy)*F(ix,M-1,iz,i)[dest];
      
      // Finally, assign repeated last element for PBC
      F(ix,M,iz,i)[dest] = F(ix,0,iz,i)[dest];
    }
}



void ComplexMultiTricubicSpline::UpdateZ(int source, int dest, int i)
{
  Grid &z = *Zgrid;
  Array<double,1> mu(Nz);
  
  ///////////////
  // Real part //
  ///////////////
  // Loop over all x and y
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++) {
      // Set up tridiagonal set of equations
      // Initialize RHS of equations
      F(ix,iy,0,i)[dest] = 1.5*(F(ix,iy,1,i)[source]-F(ix,iy,0,i)[source])/(z(1)-z(0));
      F(ix,iy,Nz-1,i)[dest] = 1.5*(F(ix,iy,Nz-1,i)[source]-F(ix,iy,Nz-2,i)[source])
	/(z(Nz-1)-z(Nz-2));
      mu(0) = 0.5;

      // Solve tri-diagonal set of equations.  First eliminate lower
      // elements.
      for (int j=1; j<(Nz-1); j++) {
	double lambda = 0.5*(z(j+1)-z(j))/(z(j+1)-z(j-1));
	mu(j) = 0.5 - lambda;
	F(ix,iy,j,i)[dest] = 
	  3.0*(lambda*(F(ix,iy,j,i)[source]-F(ix,iy,j-1,i)[source])/(z(j)-z(j-1))+ 
	       mu(j) *(F(ix,iy,j+1,i)[source]-F(ix,iy,j,i)[source])/(z(j+1)-z(j)));
	double cj = 1.0 - lambda * mu(j-1);
	F(ix,iy,j,i)[dest] -= lambda * F(ix,iy,j-1,i)[dest];
	mu(j) /= cj;
	F(ix,iy,j,i)[dest] /= cj;
      }
      
      // Last element
      int j = Nz-1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      F(ix,iy,j,i)[dest] = 
	3.0*(lambda*(F(ix,iy,j,i)[source]-F(ix,iy,j-1,i)[source])/(z(j)-z(j-1)));
      double cj = 1.0 - lambda * mu(j-1);
      F(ix,iy,j,i)[dest] -= lambda * F(ix,iy,j-1,i)[dest];
      mu(j) /= cj;
      F(ix,iy,j,i)[dest] /= cj;
      
      // Now last d/dx is correct.  We proceed upward, back substituting.
      for (j=Nz-2; j>=0; j--) {
	F(ix,iy,j,i)[dest] -= mu(j) * F(ix,iy,j+1,i)[dest];
      }
    }
}


void
ComplexMultiTricubicSpline::UpdateZPeriodic(int source, int dest, int i)
{
  Grid &z = *Zgrid;
  int M = z.NumPoints - 1;
  // Check to make sure first and last point are the same.
  Array<double,1> lambda(M), mu(M), d(M), gamma(M);
  
  ///////////////
  // Real part //
  ///////////////
  for (int ix=0; ix<Nx; ix++) 
    for (int iy=0; iy<Nz; iy++) {
      assert (F(ix,iy,0,i)[source] == F(ix,iy,M,i)[source]);
  
      // Setup lambdas, mus, and d's
      lambda(0) = 0.5*(z(1)-z(0))/(z(M)-z(M-2));
      mu(0)     = 0.5-lambda(0);
      d(0)      = 3.0*(lambda(0)*(F(ix,iy,0,i)[source]-F(ix,iy,M-1,i)[source])/(z(M)-z(M-1)) +
		       mu(0)    *(F(ix,iy,1,i)[source]-F(ix,iy,0,i)[source])  /(z(1)-z(0)));
      lambda(M-1) = 0.5*(z(M)-z(M-1))/(z(M)-z(M-2));
      mu(M-1)     = 0.5*lambda(M-1);
      d(M-1)      = 3.0*(lambda(M-1)*(F(ix,iy,M-1,i)[source]-F(ix,iy,M-2,i)[source])/(z(M-1)-z(M-2)) +
			 mu(M-1)    *(F(ix,iy,M,i)[source]  -F(ix,iy,M-1,i)[source])/(z(M)  -z(M-1)));
      gamma(M-1)  = 1.0;
      gamma(0) = lambda(0);
      
      for (int iz=1; iz<M; iz++) {
	lambda(iz) = 0.5*(z(iz+1)-z(iz))/(z(iz+1)-z(iz-1));
	mu(iz)     = 0.5-lambda(iz);
	d(iz)      = 3.0*(lambda(iz)*(F(ix,iy,iz,i)[source]-F(ix,iy,iz-1,i)[source])/(z(iz)-z(iz-1)) +
			 mu(iz)   *(F(ix,iy,iz+1,i)[source]-F(ix,iy,iz,i)[source])/(z(iz+1)-z(iz)));
      }
      // Solve down lower triangular part
      for (int iz=1; iz<(M-1); iz++) {
	double diag = 1.0-lambda(iz)*mu(iz-1);
	gamma(iz) = -gamma(iz-1)*lambda(iz);
	d(iz) -= lambda(iz) * d(iz-1);
	double diagInv = 1.0/diag;
	gamma(iz) *= diagInv;
	d(iz)     *= diagInv;
	mu(iz)    *= diagInv;
	// last row
	d(M-1) -= mu(M-1) * d(iz-1);
	gamma(M-1) -= mu(M-1) * gamma(iz-1);
	mu(M-1) = -mu(M-1)*mu(iz-1);
      }
      // last row
      // mu is really on top of lambda in the last row
      lambda(M-1) += mu(M-1);
      d(M-1) -= lambda(M-1) * d(M-2);
      gamma(M-1) -= lambda(M-1) * (gamma(M-2)+ mu(M-2));
      // Compute last derivative;
      F(ix,iy,M-1,i)[dest] = d(M-1)/gamma(M-1);
      
      // Now proceed up upper diagonal, backsubstituting
      for (int iz=M-2; iz>=0; iz--)
	F(ix,iy,iz,i)[dest] = d(iz) - mu(iz)*F(ix,iy,iz+1,i)[dest] - gamma(iz)*F(ix,iy,M-1,i)[dest];
      
      // Finally, assign repeated last element for PBC
      F(ix,iy,M,i)[dest] = F(ix,iy,0,i)[dest];
    }
}


void 
MakePeriodic (Array<complex<double>,4> &A)
{
  int Nx = A.extent(0);
  int Ny = A.extent(1);
  int Nz = A.extent(2);
  int M  = A.extent(3);
  
  for (int i=0; i<M; i++) {
    // X face
    for (int iy=0; iy<Ny-1; iy++)
      for (int iz=0; iz<Nz-1; iz++)
	A(Nx-1,iy,iz,i) = A(0,iy,iz,i);
    // Y face
    for (int ix=0; ix<Nx-1; ix++)
      for (int iz=0; iz<Nz-1; iz++)
	A(ix,Ny-1,iz,i) = A(ix,0,iz,i);
    // Z face
    for (int ix=0; ix<Nx-1; ix++)
      for (int iy=0; iy<Ny-1; iy++)
	A(ix,iy,Nz-1,i) = A(ix,iy,0,i);
    // XY edge
    for (int iz=0; iz<Nz-1; iz++)
      A(Nx-1,Ny-1,iz,i) = A(0,0,iz,i);
    // XZ edge
    for (int iy=0; iy<Ny-1; iy++)
      A(Nx-1,iy,Nz-1,i) = A(0,iy,0,i);
    // YZ edge
    for (int ix=0; ix<Nx-1; ix++)
      A(ix,Ny-1,Nz-1,i) = A(ix,0,0,i);
    /// Corner
    A(Nx-1,Ny-1,Nz-1,i) = A(0,0,0,i);
  }
}
