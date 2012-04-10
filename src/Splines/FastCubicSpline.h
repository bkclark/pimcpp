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

#ifndef FAST_CUBIC_SPLINE_H
#define FAST_CUBIC_SPLINE_H

#include <iostream>
#include "../Blitz.h"

/// The CubicSpline class is a third-order spline representation of a
/// function.  It stores a pointer to a grid and the values of the
/// function and its second derivative at the points defined by the
/// grid. 
class FastCubicSpline
{
private:
  /// This flag records whether or not the stored second derivatives
  /// are in sync with the function values.  It is used to determine
  /// whether the second derivatives need recomputation.
  bool UpToDate;
  /// The function values and first derivatives on the grid points.
  /// The values are in the [0] elements and the second derivs in [1]
  Array<Vec2, 1> F;   
  double dx, dxInv;
  double Xstart, Xend;

  /// The values of the derivative of the represented function on the
  /// boundary.  If each value is greater that 1e30, we compute
  /// bondary conditions assuming that the second derivative is zero at
  /// that boundary.
  double StartDeriv, EndDeriv;

  /// Recompute the second derivatives from the function values
  void UpdateFixed();
  void UpdateNatural();
  void UpdatePeriodic();
  /// Stores whether to use periodic boundary conditions.
  bool Periodic, Fixed;
public:
  inline int size() { return F.size(); }
  /// Returns the interpolated value.
  inline double operator()(double x);
  /// Returns the interpolated first derivative.
  inline double d_dx(double x);
  /// Returns the interpolated second derivative.
  inline double d2_dx2(double x);
  /// Returns the interpolated third derivative.
  inline double d3_dx3(double x);
  
  /// Initialize the cubic spline.  See notes about start and end
  /// deriv above.
  void Init(double xStart, double xEnd, Array<double,1> NewYs,
	    double startderiv, double endderiv);
  void Init(double xStart, double xEnd, Array<double,1> NewYs,
	    bool isPeriodic=false);
  
  /// Returns the value of the function at the ith grid point.
  inline double operator()(int i) const
  { return F(i)[0];  }

  /// Returns a reference to the value at the ith grid point.
  inline double & operator()(int i)
  { UpToDate = false;  return F(i)[0];  }

  /// Trivial constructor
  FastCubicSpline() : UpToDate (false), Periodic(false)
  {
    /// do nothing
  }
};


inline 
double FastCubicSpline::operator()(double x)
{
  if (!UpToDate) {
    if (Fixed)
      UpdateFixed();
    else if (Periodic)
      UpdatePeriodic();
    else
      UpdateNatural();
  }
       

  double xpos, t, p1, p2, q1, q2, tm1, t2;
  xpos = floor((x-Xstart)*dxInv);
  int ix = (int) xpos;
  t = (x*dxInv-xpos);
  tm1 = t-1.0;
  t2 = t*t;
  p1 = tm1*tm1*(1.0+2.0*t);
  p2 = t2*(3.0-2.0*t);
  q1 = dx*t*tm1*tm1;
  q2 = dx*t2*tm1;

  return F(ix)[0]*p1 + F(ix)[1]*q1 + F(ix+1)[0]*p2 + F(ix+1)[1]*q2;
}

#endif
