#ifndef CUBIC_B_SPLINE_H
#define CUBIC_B_SPLINE_H

#include <blitz/array.h>
#include <blitz/tinymat.h>
#include "BsplineHelper.h"

using namespace blitz;

class CubicBspline
{
private:
  TinyMatrix<double,4,4> A, dA, d2A, d3A;
  double GridStart, GridEnd, GridDelta, GridDeltaInv, L, Linv;
  bool Interpolating, Periodic;

  // The control points
  Array<double,1> P;

  // Interpolating solvers:
  void SolvePeriodicInterp(Array<double,1> &data);  

  // Data related to locating grid points
  mutable int i0, i1, i2, i3;
  mutable double t;
  mutable double tp[4];
  void Find (double x) const;

public:
  inline double GetControlPoint (int i);
  void Init(double start, double end, Array<double,1> &data, 
	    bool interpolating=true, 
	    BoundaryCondition<double> startBC=BoundaryCondition<double>(PERIODIC),
	    BoundaryCondition<double>   endBC=BoundaryCondition<double>(PERIODIC));
  inline double operator()(double x) const;
  inline double Deriv     (double x) const;
  inline double Deriv2    (double x) const;
  inline double Deriv3    (double x) const;

  CubicBspline();
};

inline void
CubicBspline::Find(double x) const
{
  double delta = x - GridStart;
  // Enforce PBC
  if (Periodic) 
    delta -= nearbyint(delta*Linv)*L;
    
  double fi = delta*GridDeltaInv;
  double ipart;
  t = modf (fi, &ipart);
  int i = (int) ipart;
  i0 = i;
  i1 = i+1;
  i2 = i+2;
  i3 = i+3;
  
  tp[0] = t*t*t;
  tp[1] = t*t;
  tp[2] = t;
  tp[3] = 1.0;
}

inline double 
CubicBspline::operator()(double x) const
{
  Find(x);
//   return (P(i0)*(A(0,0)*tp[0]+A(1,0)*tp[1]+A(2,0)*tp[2]+A(3,0)*tp[3])+
//   	     P(i  )*(A(0,1)*tp[0]+A(1,1)*tp[1]+A(2,1)*tp[2]+A(3,1)*tp[3])+
// 	     P(ip1)*(A(0,2)*tp[0]+A(1,2)*tp[1]+A(2,2)*tp[2]+A(3,2)*tp[3])+
// 	     P(ip2)*(A(0,3)*tp[0]+A(1,3)*tp[1]+A(2,3)*tp[2]+A(3,3)*tp[3]));
  return 
    (tp[0]*(A(0,0)*P(i0)+A(0,1)*P(i1)+A(0,2)*P(i2)+A(0,3)*P(i3))+
     tp[1]*(A(1,0)*P(i0)+A(1,1)*P(i1)+A(1,2)*P(i2)+A(1,3)*P(i3))+
     tp[2]*(A(2,0)*P(i0)+A(2,1)*P(i1)+A(2,2)*P(i2)+A(2,3)*P(i3))+
     tp[3]*(A(3,0)*P(i0)+A(3,1)*P(i1)+A(3,2)*P(i2)+A(3,3)*P(i3)));
}


inline double 
CubicBspline::Deriv(double x) const
{
  Find(x);
  return GridDeltaInv *
    (tp[1]*(dA(1,0)*P(i0)+dA(1,1)*P(i1)+dA(1,2)*P(i2)+dA(1,3)*P(i3))+
     tp[2]*(dA(2,0)*P(i0)+dA(2,1)*P(i1)+dA(2,2)*P(i2)+dA(2,3)*P(i3))+
     tp[3]*(dA(3,0)*P(i0)+dA(3,1)*P(i1)+dA(3,2)*P(i2)+dA(3,3)*P(i3)));
}


inline double 
CubicBspline::Deriv2(double x) const
{
  Find(x);
  return GridDeltaInv * GridDeltaInv*
    (tp[2]*(d2A(2,0)*P(i0)+d2A(2,1)*P(i1)+d2A(2,2)*P(i2)+d2A(2,3)*P(i3))+
     tp[3]*(d2A(3,0)*P(i0)+d2A(3,1)*P(i1)+d2A(3,2)*P(i2)+d2A(3,3)*P(i3)));
}

inline double 
CubicBspline::Deriv3(double x) const
{
  Find(x);
  return GridDeltaInv * GridDeltaInv* GridDeltaInv*
    (tp[3]*(d2A(3,0)*P(i0)+d2A(3,1)*P(i1)+d2A(3,2)*P(i2)+d2A(3,3)*P(i3)));
}




#endif
