/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "BicubicSpline.h"


void BicubicSpline::XUpdate(int iy)
{
  Grid &x = *Xgrid;
  Array<double,1> mu(Nx);
  
  // Set up tridiagonal set of equations
  // Initialize RHS of equations
  F(0,iy).dzdx = 1.5*(F(1,iy).z-F(0,iy).z)/(x(1)-x(0));
  F(Nx-1,iy).dzdx = 1.5*(F(Nx-1,iy).z-F(Nx-2,iy).z)/(x(Nx-1)-x(Nx-2));
  mu(0) = 0.5;

  // Solve tri-diagonal set of equations.  First eliminate lower
  // elements.
  for (int j=1; j<(Nx-1); j++) {
    double lambda = 0.5*(x(j+1)-x(j))/(x(j+1)-x(j-1));
    mu(j) = 0.5 - lambda;
    F(j,iy).dzdx = 3.0*(lambda*(F(j,iy).z-F(j-1,iy).z)/(x(j)-x(j-1)) + 
			mu(j) *(F(j+1,iy).z-F(j,iy).z)/(x(j+1)-x(j)));
    double cj = 1.0 - lambda * mu(j-1);
    F(j,iy).dzdx -= lambda * F(j-1,iy).dzdx;
    mu(j) /= cj;
    F(j,iy).dzdx /= cj;
  }
  // Last element
  int j = Nx-1;
  double lambda = 0.5;
  mu(j) = 0.5 - lambda;
  F(j,iy).dzdx = 3.0*(lambda*(F(j,iy).z-F(j-1,iy).z)/(x(j)-x(j-1)));
  double cj = 1.0 - lambda * mu(j-1);
  F(j,iy).dzdx -= lambda * F(j-1,iy).dzdx;
  mu(j) /= cj;
  F(j,iy).dzdx /= cj;

  // Now last dzdx is correct.  We proceed upward, back substituting.
  for (j=Nx-2; j>=0; j--) {
    F(j,iy).dzdx -= mu(j) * F(j+1,iy).dzdx;
  }
  XUpToDate(iy) = true;
}





void BicubicSpline::YUpdate(int ix)
{
  Grid &y = *Ygrid;
  Array<double,1> mu(Ny);
  
  // Set up tridiagonal set of equations
  // Initialize RHS of equations
  F(ix,0).dzdy = 1.5*(F(ix,1).z-F(ix,0).z)/(y(1)-y(0));
  F(ix,Ny-1).dzdy = 1.5*(F(ix,Ny-1).z-F(ix,Ny-2).z)/(y(Ny-1)-y(Ny-2));
  mu(0) = 0.5;

  // Solve tri-diagonal set of equations.  First eliminate lower
  // elements.
  for (int j=1; j<(Ny-1); j++) {
    double lambda = 0.5*(y(j+1)-y(j))/(y(j+1)-y(j-1));
    mu(j) = 0.5 - lambda;
    F(ix,j).dzdy = 3.0*(lambda*(F(ix,j).z-F(ix,j-1).z)/(y(j)-y(j-1)) + 
			mu(j) *(F(ix,j+1).z-F(ix,j).z)/(y(j+1)-y(j)));
    double cj = 1.0 - lambda * mu(j-1);
    F(ix,j).dzdy -= lambda * F(ix,j-1).dzdy;
    mu(j) /= cj;
    F(ix,j).dzdy /= cj;
  }
  // Last element
  int j = Ny-1;
  double lambda = 0.5;
  mu(j) = 0.5 - lambda;
  F(ix,j).dzdy = 3.0*(lambda*(F(ix,j).z-F(ix,j-1).z)/(y(j)-y(j-1)));
  double cj = 1.0 - lambda * mu(j-1);
  F(ix,j).dzdy -= lambda * F(ix,j-1).dzdy;
  mu(j) /= cj;
  F(ix,j).dzdy /= cj;

  // Now last dzdy is correct.  We proceed upward, back substituting.
  for (j=Ny-2; j>=0; j--) {
    F(ix,j).dzdy -= mu(j) * F(ix,j+1).dzdy;
  }
  YUpToDate(ix) = true;
}




void BicubicSpline::BiUpdate()
{
  // First, update X and Y derivatives
  for (int i=0; i<Ny; i++)
    if (!XUpToDate(i)) XUpdate(i);
  for (int i=0; i<Nx; i++)
    if (!YUpToDate(i)) YUpdate(i);

  Grid &x = *Xgrid;
  Array<double,1> mu(Nx);
  for (int iy=0; iy<Ny; iy++) {
    // Set up tridiagonal set of equations
    // Initialize RHS of equations
    F(0,iy).d2zdxdy = 1.5*(F(1,iy).dzdy-F(0,iy).dzdy)/(x(1)-x(0));
    F(Nx-1,iy).d2zdxdy = 
      1.5*(F(Nx-1,iy).dzdy-F(Nx-2,iy).dzdy)/(x(Nx-1)-x(Nx-2));
    mu(0) = 0.5;
    
    // Solve tri-diagonal set of equations.  First eliminate lower
    // elements.
    for (int j=1; j<(Nx-1); j++) {
      double lambda = 0.5*(x(j+1)-x(j))/(x(j+1)-x(j-1));
      mu(j) = 0.5 - lambda;
      F(j,iy).d2zdxdy = 
	3.0*(lambda*(F(j,iy).dzdy-F(j-1,iy).dzdy)/(x(j)-x(j-1))
	     +mu(j) *(F(j+1,iy).dzdy-F(j,iy).dzdy)/(x(j+1)-x(j)));
      double cj = 1.0 - lambda * mu(j-1);
      F(j,iy).d2zdxdy -= lambda * F(j-1,iy).d2zdxdy;
      mu(j) /= cj;
      F(j,iy).d2zdxdy /= cj;
    }
    // Last element
    int j = Nx-1;
    double lambda = 0.5;
    mu(j) = 0.5 - lambda;
    F(j,iy).d2zdxdy = 3.0*(lambda*(F(j,iy).dzdy-F(j-1,iy).dzdy)/(x(j)-x(j-1)));
    double cj = 1.0 - lambda * mu(j-1);
    F(j,iy).d2zdxdy -= lambda * F(j-1,iy).d2zdxdy;
    mu(j) /= cj;
    F(j,iy).d2zdxdy /= cj;
    
    // Now last d2zdxdy is correct.  We proceed upward, back substituting.
    for (j=Nx-2; j>=0; j--) {
      F(j,iy).d2zdxdy -= mu(j) * F(j+1,iy).d2zdxdy;
    }
  }
  BiUpToDate = true;
}





void SymmBicubicSpline::XUpdate(int iy)
{
  Grid &x = *Xgrid;
  Array<double,1> mu(Nx);
  
  // Set up tridiagonal set of equations
  // Initialize RHS of equations
  dx(0,iy) = 1.5*(z(1,iy)-z(0,iy))/(x(1)-x(0));
  dx(Nx-1,iy) = 1.5*(z(Nx-1,iy)-z(Nx-2,iy))/(x(Nx-1)-x(Nx-2));
  mu(0) = 0.5;

  // Solve tri-diagonal set of equations.  First eliminate lower
  // elements.
  for (int j=1; j<(Nx-1); j++) {
    double lambda = 0.5*(x(j+1)-x(j))/(x(j+1)-x(j-1));
    mu(j) = 0.5 - lambda;
    dx(j,iy) = 3.0*(lambda*(z(j,iy)-z(j-1,iy))/(x(j)-x(j-1)) + 
			mu(j) *(z(j+1,iy)-z(j,iy))/(x(j+1)-x(j)));
    double cj = 1.0 - lambda * mu(j-1);
    dx(j,iy) -= lambda * dx(j-1,iy);
    mu(j) /= cj;
    dx(j,iy) /= cj;
  }
  // Last element
  int j = Nx-1;
  double lambda = 0.5;
  mu(j) = 0.5 - lambda;
  dx(j,iy) = 3.0*(lambda*(z(j,iy)-z(j-1,iy))/(x(j)-x(j-1)));
  double cj = 1.0 - lambda * mu(j-1);
  dx(j,iy) -= lambda * dx(j-1,iy);
  mu(j) /= cj;
  dx(j,iy) /= cj;

  // Now last dzdx is correct.  We proceed upward, back substituting.
  for (j=Nx-2; j>=0; j--) {
    dx(j,iy) -= mu(j) * dx(j+1,iy);
  }
  XUpToDate(iy) = true;
}





void SymmBicubicSpline::YUpdate(int ix)
{
  Grid &y = *Ygrid;
  Array<double,1> mu(Ny);
  
  // Set up tridiagonal set of equations
  // Initialize RHS of equations
  dy(ix,0) = 1.5*(z(ix,1)-z(ix,0))/(y(1)-y(0));
  dy(ix,Ny-1) = 1.5*(z(ix,Ny-1)-z(ix,Ny-2))/(y(Ny-1)-y(Ny-2));
  mu(0) = 0.5;

  // Solve tri-diagonal set of equations.  First eliminate lower
  // elements.
  for (int j=1; j<(Ny-1); j++) {
    double lambda = 0.5*(y(j+1)-y(j))/(y(j+1)-y(j-1));
    mu(j) = 0.5 - lambda;
    dy(ix,j) = 3.0*(lambda*(z(ix,j)-z(ix,j-1))/(y(j)-y(j-1)) + 
			mu(j) *(z(ix,j+1)-z(ix,j))/(y(j+1)-y(j)));
    double cj = 1.0 - lambda * mu(j-1);
    dy(ix,j) -= lambda * dy(ix,j-1);
    mu(j) /= cj;
    dy(ix,j) /= cj;
  }
  // Last element
  int j = Ny-1;
  double lambda = 0.5;
  mu(j) = 0.5 - lambda;
  dy(ix,j) = 3.0*(lambda*(z(ix,j)-z(ix,j-1))/(y(j)-y(j-1)));
  double cj = 1.0 - lambda * mu(j-1);
  dy(ix,j) -= lambda * dy(ix,j-1);
  mu(j) /= cj;
  dy(ix,j) /= cj;

  // Now last dzdy is correct.  We proceed upward, back substituting.
  for (j=Ny-2; j>=0; j--) {
    dy(ix,j) -= mu(j) * dy(ix,j+1);
  }
  YUpToDate(ix) = true;
}




void SymmBicubicSpline::BiUpdate()
{
  // First, update X and Y derivatives
  for (int i=0; i<Ny; i++)
    if (!XUpToDate(i)) XUpdate(i);
  for (int i=0; i<Nx; i++)
    if (!YUpToDate(i)) YUpdate(i);

  Grid &x = *Xgrid;
  Array<double,1> mu(Nx);
  for (int iy=0; iy<Ny; iy++) {
    // Set up tridiagonal set of equations
    // Initialize RHS of equations
    dxdy(0,iy) = 1.5*(dy(1,iy)-dy(0,iy))/(x(1)-x(0));
    dxdy(Nx-1,iy) = 
      1.5*(dy(Nx-1,iy)-dy(Nx-2,iy))/(x(Nx-1)-x(Nx-2));
    mu(0) = 0.5;
    
    // Solve tri-diagonal set of equations.  First eliminate lower
    // elements.
    for (int j=1; j<(Nx-1); j++) {
      double lambda = 0.5*(x(j+1)-x(j))/(x(j+1)-x(j-1));
      mu(j) = 0.5 - lambda;
      dxdy(j,iy) = 
	3.0*(lambda*(dy(j,iy)-dy(j-1,iy))/(x(j)-x(j-1))
	     +mu(j) *(dy(j+1,iy)-dy(j,iy))/(x(j+1)-x(j)));
      double cj = 1.0 - lambda * mu(j-1);
      dxdy(j,iy) -= lambda * dxdy(j-1,iy);
      mu(j) /= cj;
      dxdy(j,iy) /= cj;
    }
    // Last element
    int j = Nx-1;
    double lambda = 0.5;
    mu(j) = 0.5 - lambda;
    dxdy(j,iy) = 3.0*(lambda*(dy(j,iy)-dy(j-1,iy))/(x(j)-x(j-1)));
    double cj = 1.0 - lambda * mu(j-1);
    dxdy(j,iy) -= lambda * dxdy(j-1,iy);
    mu(j) /= cj;
    dxdy(j,iy) /= cj;
    
    // Now last d2zdxdy is correct.  We proceed upward, back substituting.
    for (j=Nx-2; j>=0; j--) {
      dxdy(j,iy) -= mu(j) * dxdy(j+1,iy);
    }
  }
  BiUpToDate = true;
}





void MultiBicubicSpline::XUpdate(int iy)
{
  Grid &x = *Xgrid;
  Array<double,1> mu(Nx);
  
  for (int iz=0; iz<Nz; iz++) {
    // Set up tridiagonal set of equations
    // Initialize RHS of equations
    F(0,iy,iz).dzdx = 1.5*(F(1,iy,iz).z-F(0,iy,iz).z)/(x(1)-x(0));
    F(Nx-1,iy,iz).dzdx = 1.5*(F(Nx-1,iy,iz).z-F(Nx-2,iy,iz).z)/(x(Nx-1)-x(Nx-2));
    mu(0) = 0.5;
    
    // Solve tri-diagonal set of equations.  First eliminate lower
    // elements.
    for (int j=1; j<(Nx-1); j++) {
      double lambda = 0.5*(x(j+1)-x(j))/(x(j+1)-x(j-1));
      mu(j) = 0.5 - lambda;
      F(j,iy,iz).dzdx = 
	3.0*(lambda*(F(j,iy,iz).z-F(j-1,iy,iz).z)/(x(j)-x(j-1)) + 
	     mu(j) *(F(j+1,iy,iz).z-F(j,iy,iz).z)/(x(j+1)-x(j)));
      double cj = 1.0 - lambda * mu(j-1);
      F(j,iy,iz).dzdx -= lambda * F(j-1,iy,iz).dzdx;
      mu(j) /= cj;
      F(j,iy,iz).dzdx /= cj;
    }
    // Last element
    int j = Nx-1;
    double lambda = 0.5;
    mu(j) = 0.5 - lambda;
    F(j,iy,iz).dzdx = 3.0*(lambda*(F(j,iy,iz).z-F(j-1,iy,iz).z)/(x(j)-x(j-1)));
    double cj = 1.0 - lambda * mu(j-1);
    F(j,iy,iz).dzdx -= lambda * F(j-1,iy,iz).dzdx;
    mu(j) /= cj;
    F(j,iy,iz).dzdx /= cj;
    
    // Now last dzdx is correct.  We proceed upward, back substituting.
    for (j=Nx-2; j>=0; j--) {
      F(j,iy,iz).dzdx -= mu(j) * F(j+1,iy,iz).dzdx;
    }
  }
  XUpToDate(iy) = true;
}





void MultiBicubicSpline::YUpdate(int ix)
{
  Grid &y = *Ygrid;
  Array<double,1> mu(Ny);

  for (int iz=0; iz<Nz; iz++) {
    // Set up tridiagonal set of equations
    // Initialize RHS of equations
    F(ix,0,iz).dzdy = 1.5*(F(ix,1,iz).z-F(ix,0,iz).z)/(y(1)-y(0));
    F(ix,Ny-1,iz).dzdy = 
      1.5*(F(ix,Ny-1,iz).z-F(ix,Ny-2,iz).z)/(y(Ny-1)-y(Ny-2));
    mu(0) = 0.5;
    
    // Solve tri-diagonal set of equations.  First eliminate lower
    // elements.
    for (int j=1; j<(Ny-1); j++) {
      double lambda = 0.5*(y(j+1)-y(j))/(y(j+1)-y(j-1));
      mu(j) = 0.5 - lambda;
      F(ix,j,iz).dzdy = 
	3.0*(lambda*(F(ix,j,iz).z-F(ix,j-1,iz).z)/(y(j)-y(j-1)) + 
	     mu(j) *(F(ix,j+1,iz).z-F(ix,j,iz).z)/(y(j+1)-y(j)));
      double cj = 1.0 - lambda * mu(j-1);
      F(ix,j,iz).dzdy -= lambda * F(ix,j-1,iz).dzdy;
      mu(j) /= cj;
      F(ix,j,iz).dzdy /= cj;
    }
    // Last element
    int j = Ny-1;
    double lambda = 0.5;
    mu(j) = 0.5 - lambda;
    F(ix,j,iz).dzdy = 3.0*(lambda*(F(ix,j,iz).z-F(ix,j-1,iz).z)/(y(j)-y(j-1)));
    double cj = 1.0 - lambda * mu(j-1);
    F(ix,j,iz).dzdy -= lambda * F(ix,j-1,iz).dzdy;
    mu(j) /= cj;
    F(ix,j,iz).dzdy /= cj;
    
    // Now last dzdy is correct.  We proceed upward, back substituting.
    for (j=Ny-2; j>=0; j--) {
      F(ix,j,iz).dzdy -= mu(j) * F(ix,j+1,iz).dzdy;
    }
  }
  YUpToDate(ix) = true;
}




void MultiBicubicSpline::BiUpdate()
{
  // First, update X and Y derivatives
  for (int i=0; i<Ny; i++)
    if (!XUpToDate(i)) XUpdate(i);
  for (int i=0; i<Nx; i++)
    if (!YUpToDate(i)) YUpdate(i);

  Grid &x = *Xgrid;
  Array<double,1> mu(Nx);
  for (int iz=0; iz<Nz; iz++) {
    for (int iy=0; iy<Ny; iy++) {
      // Set up tridiagonal set of equations
      // Initialize RHS of equations
      F(0,iy,iz).d2zdxdy = 1.5*(F(1,iy,iz).dzdy-F(0,iy,iz).dzdy)/(x(1)-x(0));
      F(Nx-1,iy,iz).d2zdxdy = 
	1.5*(F(Nx-1,iy,iz).dzdy-F(Nx-2,iy,iz).dzdy)/(x(Nx-1)-x(Nx-2));
      mu(0) = 0.5;
      
      // Solve tri-diagonal set of equations.  First eliminate lower
      // elements.
      for (int j=1; j<(Nx-1); j++) {
	double lambda = 0.5*(x(j+1)-x(j))/(x(j+1)-x(j-1));
	mu(j) = 0.5 - lambda;
	F(j,iy,iz).d2zdxdy = 
	  3.0*(lambda*(F(j,iy,iz).dzdy-F(j-1,iy,iz).dzdy)/(x(j)-x(j-1))
	       +mu(j) *(F(j+1,iy,iz).dzdy-F(j,iy,iz).dzdy)/(x(j+1)-x(j)));
	double cj = 1.0 - lambda * mu(j-1);
	F(j,iy,iz).d2zdxdy -= lambda * F(j-1,iy,iz).d2zdxdy;
	mu(j) /= cj;
	F(j,iy,iz).d2zdxdy /= cj;
      }
      // Last element
      int j = Nx-1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      F(j,iy,iz).d2zdxdy = 
	3.0*(lambda*(F(j,iy,iz).dzdy-F(j-1,iy,iz).dzdy)/(x(j)-x(j-1)));
      double cj = 1.0 - lambda * mu(j-1);
      F(j,iy,iz).d2zdxdy -= lambda * F(j-1,iy,iz).d2zdxdy;
      mu(j) /= cj;
      F(j,iy,iz).d2zdxdy /= cj;
      
      // Now last d2zdxdy is correct.  We proceed upward, back substituting.
      for (j=Nx-2; j>=0; j--) {
	F(j,iy,iz).d2zdxdy -= mu(j) * F(j+1,iy,iz).d2zdxdy;
      }
    }
  }
  BiUpToDate = true;
}

  



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



/*void
BicubicSpline::XUpdate(int i)
{
  Grid &x = *Xgrid;
  int N = x.NumPoints;

  Array<double,1> U(N);

  if (isnan(XStartDeriv(i))) {
    d2Fdx2(0,i)=0.0; // Use "Natural" boundary conditions--ie d^2F/dx^2 = 0 
    U(0) = 0.0;
  }
  else {
    d2Fdx2(0,i) = -0.5;
    U(0) = (3.0/(x(1)-x(0)))*((F(1,i)-F(0,i))/(x(1)-x(0))-XStartDeriv(i));
  }

  d2Fdx2(Nx-1,i) = 0.0;
  
  // This part translated from Numerical Recipes.
  for (int j=1; j<Nx-1; j++) {
    double sig = (x(j) - x(j-1)) / (x(j+1)-x(j-1));
    double p = sig *d2Fdx2(j-1,i)+2.0;
    d2Fdx2(j,i) = (sig-1.0)/p;
    U(j) = (6.0*((F(j+1,i)-F(j,i))/(x(j+1)-x(j))-
		 (F(j,i)-F(j-1,i))/(x(j)-x(j-1))) /
	    (x(j+1)-x(j-1))-sig*U(j-1))/p;   
  }

  double Un, Qn;
  if (XEndDeriv(i) > 0.99e30)
    {
      Un = 0.0;
      Qn = 0.0;
    }
  else
    {
      Qn = 0.5;
      Un = (3.0/(x(N-1)-x(N-2)))*
	(XEndDeriv(i)-(F(N-1,i)-F(N-2,i))/(x(N-1)-x(N-2)));
    }
  
  d2Fdx2(N-1,i) =
    (Un-Qn*U(N-2))/(Qn*d2Fdx2(N-2,i)+1.0);
  
  for (int k=N-2; k>=0; k--)
    d2Fdx2(k,i) = d2Fdx2(k,i)*d2Fdx2(k+1,i) + U(k);

  XUpToDate(i) = true;
}
*/


/*
void
BicubicSpline::YUpdate(int i)
{
  Grid &y = *Ygrid;
  int N = y.NumPoints;

  Array<double,1> U(N);

  if (isnan(YStartDeriv(i))) {
    d2Fdy2(i,0)=0.0; // Use "Natural" boundary conditions--ie d^2F/dy^2 = 0 
    U(0) = 0.0;
  }
  else {
    d2Fdy2(i,0) = -0.5;
    U(0) = (3.0/(y(1)-y(0)))*((F(i,1)-F(i,0))/(y(1)-y(0))-YStartDeriv(i));
  }

  d2Fdy2(i,Ny-1) = 0.0;
  
  // This part translated from Numerical Recipes.
  for (int j=1; j<Ny-1; j++) {
    double sig = (y(j) - y(j-1)) / (y(j+1)-y(j-1));
    double p = sig *d2Fdy2(i,j-1)+2.0;
    d2Fdy2(i,j) = (sig-1.0)/p;
    U(j) = (6.0*((F(i,j+1)-F(i,j))/(y(j+1)-y(j))-
		 (F(i,j)-F(i,j-1))/(y(j)-y(j-1))) /
	    (y(j+1)-y(j-1))-sig*U(j-1))/p;   
  }

  double Un, Qn;
  if (YEndDeriv(i) > 0.99e30)
    {
      Un = 0.0;
      Qn = 0.0;
    }
  else
    {
      Qn = 0.5;
      Un = (3.0/(y(N-1)-y(N-2)))*
	(YEndDeriv(i)-(F(i,N-1)-F(i,N-2))/(y(N-1)-y(N-2)));
    }
  
  d2Fdy2(i,N-1) =
    (Un-Qn*U(N-2))/(Qn*d2Fdy2(i,N-2)+1.0);
  
  for (int k=N-2; k>=0; k--)
    d2Fdy2(i,k) = d2Fdy2(i,k)*d2Fdy2(i,k+1) + U(k);

  YUpToDate(i) = true;
}

*/
