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

#include "PeriodicSpline.h"

// We assume that the last y value is the same as the first.
void
PeriodicSpline::Update()
{
  Grid &x = *grid;
  int N = grid->NumPoints - 1;
  // Check to make sure first and last point are the same.
  double L = x.End-x.Start;
  assert (F(0)[0] == F(N)[0]);
  double lambda[N], mu[N], d[N], gamma[N];
  
  // Setup lambdas, mus, and d's
  lambda[0] = 0.5*(x(1)-x(0))/(x(N)-x(N-2));
  mu[0]     = 0.5-lambda[0];
  d[0]      = 3.0*(lambda[0]*(F(0)[0]-F(N-1)[0])/(x(N)-x(N-1)) +
		   mu[0]    *(F(1)[0]-F(0)[0])  /(x(1)-x(0)));
  lambda[N-1] = 0.5*(x(N)-x(N-1))/(x(N)-x(N-2));
  mu[N-1]     = 0.5*lambda[N-1];
  d[N-1]      = 3.0*(lambda[N-1]*(F(N-1)[0]-F(N-2)[0])/(x(N-1)-x(N-2)) +
		     mu[N-1]    *(F(N)[0]  -F(N-1)[0])/(x(N)  -x(N-1)));
  gamma[N-1]  = 1.0;
  gamma[0] = lambda[0];
  
  for (int i=1; i<N; i++) {
    lambda[i] = 0.5*(x(i+1)-x(i))/(x(i+1)-x(i-1));
    mu[i]     = 0.5-lambda[i];
    d[i]      = 3.0*(lambda[i]*(F(i)[0]-F(i-1)[0])/(x(i)-x(i-1)) +
		     mu[i]   *(F(i+1)[0]-F(i)[0])/(x(i+1)-x(i)));
  }
  // Solve down lower triangular part
  for (int i=1; i<(N-1); i++) {
    double diag = 1.0-lambda[i]*mu[i-1];
    gamma[i] = -gamma[i-1]*lambda[i];
    d[i] -= lambda[i] * d[i-1];
    double diagInv = 1.0/diag;
    gamma[i] *= diagInv;
    d[i]     *= diagInv;
    mu[i]    *= diagInv;
    // last row
    d[N-1] -= mu[N-1] * d[i-1];
    gamma[N-1] -= mu[N-1] * gamma[i-1];
    mu[N-1] = -mu[N-1]*mu[i-1];
  }
  // last row
  // mu is really on top of lambda in the last row
  lambda[N-1] += mu[N-1];
  d[N-1] -= lambda[N-1] * d[N-2];
  gamma[N-1] -= lambda[N-1] * (gamma[N-2]+ mu[N-2]);
  // Compute last derivative;
  F(N-1)[1] = d[N-1]/gamma[N-1];
  
  // Now proceed up upper diagonal, backsubstituting
  for (int i=N-2; i>=0; i--)
    F(i)[1] = d[i] - mu[i]*F(i+1)[1] - gamma[i]*F(N-1)[1];

  // Finally, assign repeated last element for PBC
  F(N) = F(0);
  
  IsUp2Date = true;
}

void 
PeriodicSpline::Init (Grid *newGrid, const Array<double,1> &data)
{
  grid = newGrid;
  assert (grid->NumPoints == data.size());
  F.resize(data.size());
  for (int i=0; i<data.size(); i++)
    F(i)[0] = data(i);
  Update();
}
