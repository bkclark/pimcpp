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

#include "FastCubicSpline.h"

void
FastCubicSpline::Init (double xStart, double xEnd, Array<double,1> NewYs,
		       double startDeriv, double endDeriv)
{
  int nx = NewYs.size();
  F.resize(nx);
  for (int i=0; i<F.size(); i++)
    F(i)[0] = NewYs(i);
  Xstart = xStart;
  Xend   = xEnd;
  dx = (xEnd-xStart)/(double)(nx-1);
  dxInv = 1.0/dx;
  Fixed    = true;
  Periodic = false;
  StartDeriv = startDeriv;
  EndDeriv = endDeriv;
  UpdateFixed();
}


void
FastCubicSpline::Init (double xStart, double xEnd, Array<double,1> NewYs,
		       bool isPeriodic)
{
  int nx = NewYs.size();
  F.resize(nx);
  for (int i=0; i<F.size(); i++)
    F(i)[0] = NewYs(i);
  Xstart = xStart;
  Xend   = xEnd;
  dx = (xEnd-xStart)/(double)(nx-1);
  dxInv = 1.0/dx;
  Fixed    = false;
  Periodic = isPeriodic;
  if (Periodic)
    UpdatePeriodic();
  else
    UpdateNatural();
}

 
void
FastCubicSpline::UpdatePeriodic()
{
  int M = F.extent(0)-1;
  Array<double,1> lambda(M), mu(M), gamma(M);
  Array<double,1> d(M);
  assert (F(0)[0] == F(M)[0]);
  // Setup lambdas, mus, and d's
  lambda(0) = 0.25;
  mu(0)     = 0.5-lambda(0);
  d(0)      = 3.0*(lambda(0)*(F(0)[0]-F(M-1)[0])*dxInv +
		   mu(0)    *(F(1)[0]-F(0)[0]) * dxInv);
  lambda(M-1) = 0.25;
  mu(M-1)     = 0.5*lambda(M-1);
  d(M-1)      = 3.0*(lambda(M-1)*(F(M-1)[0]-F(M-2)[0])*dxInv +
		     mu(M-1)    *(F(M)[0]  -F(M-1)[0])*dxInv);
  gamma(M-1)  = 1.0;
  gamma(0) = lambda(0);
  
  for (int ix=1; ix<M; ix++) {
    lambda(ix) = 0.25;
    mu(ix)     = 0.5-lambda(ix);
    d(ix)      = 3.0*(lambda(ix)*(F(ix)[0]-F(ix-1)[0])*dxInv +
		      mu(ix)   *(F(ix+1)[0]-F(ix)[0])*dxInv);
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
  F(M-1)[1] = d(M-1)/gamma(M-1);
  
  // Now proceed up upper diagonal, backsubstituting
  for (int ix=M-2; ix>=0; ix--)
    F(ix)[1] = d(ix) - mu(ix)*F(ix+1)[1] - gamma(ix)*F(M-1)[1];
  
  // Finally, assign repeated last element for PBC
  F(M)[1] = F(0)[1];
}

void
FastCubicSpline::UpdateFixed()
{

}

void
FastCubicSpline::UpdateNatural()
{
  if (Periodic) {
    UpdatePeriodic();
    return;
  }
  Array<double,1> mu(F.size());
  int Nx = F.size();
  // Set up tridiagonal set of equations
  // Initialize RHS of equations
  F(0)[1] = 1.5*(F(1)[0]-F(0)[0])*dxInv;
  F(Nx-1)[1] = 1.5*(F(Nx-1)[0]-F(Nx-2)[0])
    *dxInv;
  mu(0) = 0.5;
  
  // Solve tri-diagonal set of equations.  First eliminate lower
  // elements.
  for (int j=1; j<(Nx-1); j++) {
    double lambda = 0.25;
    mu(j) = 0.5 - lambda;
    F(j)[1] = 3.0*(lambda*(F(j)[0]-F(j-1)[0])*dxInv+ 
		   mu(j) *(F(j+1)[0]-F(j)[0])*dxInv);
    double cj = 1.0 - lambda * mu(j-1);
    F(j)[1] -= lambda * F(j-1)[1];
    mu(j) /= cj;
    F(j)[1] /= cj;
  }
  
  // Last element
  int j = Nx-1;
  double lambda = 0.5;
  mu(j) = 0.5 - lambda;
  F(j)[1] = 
    3.0*(lambda*(F(j)[0]-F(j-1)[0])*dxInv);
  double cj = 1.0 - lambda * mu(j-1);
  F(j)[1] -= lambda * F(j-1)[1];
  mu(j) /= cj;
  F(j)[1] /= cj;
  
  // Now last d/dx is correct.  We proceed upward, back substituting.
  for (j=Nx-2; j>=0; j--) {
    F(j)[1] -= mu(j) * F(j+1)[1];
  }
  UpToDate = true;
}

