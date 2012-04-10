#include "CubicBspline.h"
#include "BsplineHelper.h"
#include <iostream>
#include <cstdlib>

// Solve tridiagonal linear system in periodic boundary conditions
void
CubicBspline::SolvePeriodicInterp(Array<double,1> &data)
{
  double ratio = 0.25;
  int N = data.size();

  Array<double,1> d(N), gamma(N), mu(N);
  d = 1.5*data;
  P.resize(Range(0,N+2));
  // First, eliminate leading coefficients
  gamma (0) = ratio;
  mu(0) = ratio;
  mu(N-1) = ratio;
  gamma(N-1) = 1.0;
  for (int row=1; row <(N-1); row++) {
    double diag = 1.0- mu(row-1)*ratio;
    double diagInv = 1.0/diag;
    gamma(row) = -ratio*gamma(row-1)*diagInv;
    mu(row) = diagInv*ratio;
    d(row)  = diagInv*(d(row)-ratio*d(row-1));
    // Last row
    d(N-1) -= mu(N-1) * d(row-1);
    gamma(N-1) -= mu(N-1)*gamma(row-1);
    mu(N-1) = -mu(N-1)*mu(row-1);
  }
  // Last row:  gamma(N-1) hold diagonal element
  mu(N-1) += ratio;
  gamma(N-1) -= mu(N-1)*(mu(N-2)+gamma(N-2));
  d(N-1) -= mu(N-1) * d(N-2);
  P(N) = d(N-1)/gamma(N-1);
 
  // Now go back upward, back substituting
  for (int row=N-2; row>=0; row--) 
    P(row+1) = d(row) - mu(row)*P(row+2) - gamma(row)*P(N);
}





void
CubicBspline::Init(double start, double end, Array<double,1> &data,
		   bool interpolating, 
		   BoundaryCondition<double> startBC,
		   BoundaryCondition<double> endBC)
{
  Interpolating = interpolating;
  GridStart = start;
  GridEnd   = end;
  L = GridEnd - GridStart;
  Linv = 1.0/L;
  int N = data.size();
  
  if (startBC.GetType() == PERIODIC) {
    if (endBC.GetType() != PERIODIC) {
      cerr << "If startBC is periodic, endBC must also be.\n";
      abort();
    }
    Periodic = true;
    GridDelta = (end-start)/(double)(data.size());
    GridDeltaInv = 1.0/GridDelta;
    P.resize(data.size()+3);

    if (interpolating) 
      SolvePeriodicInterp (data);
    else 
      P(Range(1,N)) = data;
    // Finally, assign periodic elements
    P(0)   = P(N);
    P(N+1) = P(1);
    P(N+2) = P(2);
  }
  else {
    Periodic = false;
    if (endBC.GetType() == PERIODIC) {
     cerr << "If endBC is periodic, startBC must also be.\n";
      abort();
    }
    GridDelta = (end-start)/(double)(data.size()-1);
    GridDeltaInv = 1.0/GridDelta;
    P.resize(data.size()+2);
//     Array<double,1> d(P.size());
//     d(Range(1,N)) = data;
//     d(0) = 0.0;
//     d(N+1) = 0.0;
    if (interpolating) {
      TinyVector<double,4> rBC, lBC;
      if (startBC.GetType() == FIXED_FIRST) 
	lBC = -3.0, 0.0, 3.0, startBC.GetVal()*GridDelta;
      else if (startBC.GetType() == FIXED_SECOND) 
	lBC = 6.0, -12.0, 6.0, startBC.GetVal()*GridDelta*GridDelta;
      else if (startBC.GetType() == FLAT) 
	lBC = -3.0, 0.0, 3.0, 0.0;
      else if (startBC.GetType() == NATURAL)
	lBC = 6.0, -12.0, 6.0, 0.0;

      if (endBC.GetType() == FIXED_FIRST) 
	rBC = -3.0, 0.0, 3.0, endBC.GetVal()*GridDelta;
      else if (endBC.GetType() == FIXED_SECOND) 
	rBC = 6.0, -12.0, 6.0, endBC.GetVal()*GridDelta*GridDelta;
      else if (endBC.GetType() == FLAT) 
	rBC = -3.0, 0.0, 3.0, 0.0;
      else if (endBC.GetType() == NATURAL)
	rBC = 6.0, -12.0, 6.0, 0.0;
      

//       rBC = -3.0, 0.0, 3.0, GridDelta;
//       lBC = -3.0, 0.0, 3.0, GridDelta;
//       rBC = 6.0, -12.0, 6.0, GridDelta*GridDelta;
//       lBC = 6.0, -12.0, 6.0, GridDelta*GridDelta;
      SolveDerivInterp1D (data, P, lBC, rBC);
    }
    else {
      cerr << "Don't know how to do noninterpolating nonperiodic.\n";
      abort();
    }
  }
}


CubicBspline::CubicBspline() 
{
  A(0,0) = -1.0/6.0; A(0,1) =  3.0/6.0; A(0,2) = -3.0/6.0; A(0,3) = 1.0/6.0;
  A(1,0) =  3.0/6.0; A(1,1) = -6.0/6.0; A(1,2) =  3.0/6.0; A(1,3) = 0.0/6.0;
  A(2,0) = -3.0/6.0; A(2,1) =  0.0/6.0; A(2,2) =  3.0/6.0; A(2,3) = 0.0/6.0;
  A(3,0) =  1.0/6.0; A(3,1) =  4.0/6.0; A(3,2) =  1.0/6.0; A(3,3) = 0.0/6.0;

  dA(0,0)= 0.0; dA(0,1)= 0.0; dA(0,2)= 0.0; dA(0,3)= 0.0;
  dA(1,0)=-0.5; dA(1,1)= 1.5; dA(1,2)=-1.5; dA(1,3)= 0.5;
  dA(2,0)= 1.0; dA(2,1)=-2.0; dA(2,2)= 1.0; dA(2,3)= 0.0;
  dA(3,0)=-0.5; dA(3,1)= 0.0; dA(3,2)= 0.5; dA(3,3)= 0.0;

  d2A(0,0)= 0.0; d2A(0,1)= 0.0; d2A(0,2)= 0.0; d2A(0,3)= 0.0;
  d2A(1,0)= 0.0; d2A(1,1)= 0.0; d2A(1,2)= 0.0; d2A(1,3)= 0.0;
  d2A(2,0)=-1.0; d2A(2,1)= 3.0; d2A(2,2)=-3.0; d2A(2,3)= 1.0;
  d2A(3,0)= 1.0; d2A(3,1)=-2.0; d2A(3,2)= 1.0; d2A(3,3)= 0.0;

  d3A(0,0)= 0.0; d3A(0,1)= 0.0; d3A(0,2)= 0.0; d3A(0,3)= 0.0;
  d3A(1,0)= 0.0; d3A(1,1)= 0.0; d3A(1,2)= 0.0; d3A(1,3)= 0.0;
  d3A(2,0)= 0.0; d3A(2,1)= 0.0; d3A(1,2)= 2.0; d3A(2,3)= 0.0;
  d3A(3,0)=-1.0; d3A(3,1)= 3.0; d3A(3,2)=-3.0; d3A(3,3)= 1.0;
}
