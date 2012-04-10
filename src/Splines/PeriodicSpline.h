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

#ifndef PERIODIC_SPLINE_H
#define PERIODIC_SPLINE_H

#include "Grid.h"

class PeriodicSpline
{
private:
  // The first component has y, the second has dy/dx
  Array<Vec2,1> F;
  Grid *grid;
  void Update();
  bool IsUp2Date;
  
public:
  inline double operator()(double x);
  inline double Deriv     (double x);
  
  void Init (Grid *newGrid, const Array<double,1> &data);
  PeriodicSpline() : IsUp2Date(false)
  {
    // do nothing
  }
};

inline double
PeriodicSpline::operator()(double x)
{
  if (!IsUp2Date)
    Update();
  int i = grid->ReverseMap(x);
  double h = (*grid)(i+1) - (*grid)(i);
  double t = (x - (*grid)(i))/h;
  double tm1 = t-1.0;
  double p1 = (1.0 + 2.0*t)*tm1*tm1;
  double q1 = t*tm1*tm1;
  double p2 = t*t*(3.0-2.0*t);
  double q2 = t*t*tm1;

  return (F(i)[0]*p1 + F(i+1)[0]*p2 + h*(F(i)[1]*q1 + F(i+1)[1]*q2));
}


inline double 
PeriodicSpline::Deriv(double x)
{
  if (!IsUp2Date)
    Update();
  int i = grid->ReverseMap(x);
  double h = (*grid)(i+1) - (*grid)(i);
  double hinv = 1.0/h;
  double t = (x - (*grid)(i))/h;
  double tm1 = t-1.0;
  double dp1 = 6.0*t*tm1;
  double dq1 = tm1*(3.0*t-1.0);
  double dp2 = -dp1;
  double dq2 = 3.0*t*t-2.0*t;
  return hinv*(F(i)[0]*dp1 + F(i+1)[0]*dp2) + F(i)[1]*dq1 + F(i+1)[1]*dq2;
}


#endif
