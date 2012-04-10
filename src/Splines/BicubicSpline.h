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

#ifndef BICUBIC_SPLINE_H
#define BICUBIC_SPLINE_H

#include "Grid.h"

/// This structure holds the required data for a grid point, which
/// includes the value of z at (xi,yj), dz/dx, dz/dy, and d^2/dxdy.
struct BCpoint
{
  double z;
  double dzdx;
  double dzdy;
  double d2zdxdy;
};

/// This class is used to interpolate a function z(x,y) in two
/// dimensions.  It is of cubic order in both directions.  It
/// currently uses only "natural" boundary conditions:  all second
/// derivatives are 0 at the boundary.
class BicubicSpline
{
private:
  /// Stores whether dzdx needs to be recalculated for each column
  Array<bool,1> XUpToDate;
  /// Stores whether dzdy needs to be recalculated for each row
  Array<bool,1> YUpToDate;
  /// Stores whether the d2zdxdy's need to be updated
  bool BiUpToDate;
  /// Update the iy column's dzdx's 
  void XUpdate(int iy);
  /// Update the ix row's dzdy's
  void YUpdate(int ix);
  /// Update all the d2zdxdy's
  void BiUpdate();
  /// Holds the BCpoint array containing z and its derivatives.
  Array<BCpoint,2> F;
public:
  int Nx, Ny;
  Grid *Xgrid, *Ygrid;

  /// Interpolates in y at a given row
  inline double  operator() (int ix,   double y);
  /// Interpolates in x at a given column
  inline double  operator() (double x, int iy);
  /// Bicubic interpolation in x and y
  inline double  operator() (double x, double y);
  /// Returns z(x_{ix}, y_{iy})
  inline double  operator() (int ix, int iy) const;
  /// Returns a reference to z(x_{ix}, y_{iy}) that can be used as an
  /// L-value.
  inline double& operator() (int ix, int iy);
  inline double  Deriv      (int ix,   double y);
  inline double  Deriv      (double x, int iy);
  inline double  xDeriv     (int ix, int iy);
  inline double  yDeriv     (int ix, int iy);
  inline double  Deriv2     (int ix,   double y);
  inline double  Deriv2     (double x, int iy);
  inline double  Deriv3     (int ix,   double y);
  inline double  Deriv3     (double x, int iy);
  inline double  d_dx       (double x, double y);
  inline double  d_dy       (double x, double y);
  inline double  d2_dx2     (double x, double y);
  inline double  d2_dy2     (double x, double y);
  inline double  d2_dxdy    (double x, double y);
  /// Initialize the bicubic spline with the given grids and data
  inline void Init (Grid *xgrid, Grid *ygrid, Array<double,2> &f);

  /// Copy constructor
  inline BicubicSpline (const BicubicSpline &a);

  /// Assigment operator -- necessary for array resizeAndPreserve
  inline BicubicSpline & operator= (BicubicSpline &a);
  inline BicubicSpline & operator= (BicubicSpline a);
  inline BicubicSpline() {}
};

inline BicubicSpline::BicubicSpline (const BicubicSpline &a)
{
  XUpToDate.resize(a.XUpToDate.size());
  XUpToDate = a.XUpToDate;
  YUpToDate.resize(a.YUpToDate.size());
  YUpToDate = a.YUpToDate;
  F.resize(a.F.rows(), a.F.cols());
  F=a.F;
  Nx=a.Nx; Ny=a.Ny;
  Xgrid=a.Xgrid; Ygrid=a.Ygrid;
  BiUpToDate = a.BiUpToDate;
}

inline BicubicSpline& BicubicSpline::operator=(BicubicSpline &a)
{
  XUpToDate.resize(a.XUpToDate.size());
  XUpToDate = a.XUpToDate;
  YUpToDate.resize(a.YUpToDate.size());
  YUpToDate = a.YUpToDate;
  F.resize(a.F.rows(), a.F.cols());
  F=a.F;
  Nx=a.Nx; Ny=a.Ny;
  Xgrid=a.Xgrid; Ygrid=a.Ygrid;
  BiUpToDate = a.BiUpToDate;
  return *this;
}

inline BicubicSpline& BicubicSpline::operator=(BicubicSpline a)
{
  XUpToDate.resize(a.XUpToDate.size());
  XUpToDate = a.XUpToDate;
  YUpToDate.resize(a.YUpToDate.size());
  YUpToDate = a.YUpToDate;
  F.resize(a.F.rows(), a.F.cols());
  F=a.F;
  Nx=a.Nx; Ny=a.Ny;
  Xgrid=a.Xgrid; Ygrid=a.Ygrid;
  BiUpToDate = a.BiUpToDate;
  return *this;
}


inline void BicubicSpline::Init(Grid *xgrid, Grid *ygrid, Array<double,2> &f)
{
  Nx = xgrid->NumPoints;
  Ny = ygrid->NumPoints;
  Xgrid = xgrid; Ygrid = ygrid;
  assert (f.rows() == Nx);
  assert (f.cols() == Ny);
  XUpToDate.resize(Ny);
  YUpToDate.resize(Nx);
  F.resize(f.rows(),f.cols());
  for (int i=0; i<Nx; i++)
    for (int j=0; j<Ny; j++)
      F(i,j).z = f(i,j);
  for (int i=0; i<Ny; i++)
    XUpdate(i);
  for (int i=0; i<Nx; i++)
    YUpdate(i);

  BiUpdate();
}

inline double BicubicSpline::operator() (int ix, int iy) const
{ return (F(ix,iy).z);}

inline double& BicubicSpline::operator() (int ix, int iy)
{ 
  XUpToDate(iy) = false;
  YUpToDate(ix) = false;
  BiUpToDate = false;
  return (F(ix,iy).z);
}

inline double BicubicSpline::operator() (double x, int iy)
{
  if (!XUpToDate(iy))
    XUpdate(iy);
  
  int ix = Xgrid->ReverseMap(x);
  ix = max(0,ix);
  ix = min(ix, Nx-2);
  
  double t = (x - (*Xgrid)(ix))/((*Xgrid)(ix+1) - (*Xgrid)(ix));
  double tm1 = t - 1.0;
  double p1 = tm1*tm1*(1.0+2.0*t);
  double q1 = t*tm1*tm1;
  double p2 = t*t*(3.0-2.0*t);
  double q2 = t*t*tm1;
  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);

  return (F(ix,iy).z*p1 + F(ix+1,iy).z*p2 + 
	  h*(F(ix,iy).dzdx*q1 + F(ix+1,iy).dzdx*q2));
}


inline double BicubicSpline::operator() (int ix, double y)
{
  if (!YUpToDate(ix))
    YUpdate(ix);
  
  int iy = Ygrid->ReverseMap(y);
  iy = max(0,iy);
  iy = min(iy, Ny-2);
  
  double t = (y - (*Ygrid)(iy))/((*Ygrid)(iy+1) - (*Ygrid)(iy));
  double tm1 = t - 1.0;
  double p1 = tm1*tm1*(1.0+2.0*t);
  double q1 = t*tm1*tm1;
  double p2 = t*t*(3.0-2.0*t);
  double q2 = t*t*tm1;
  double h = (*Ygrid)(iy+1) - (*Ygrid)(iy);

  return (F(ix,iy).z*p1 + F(ix,iy+1).z*p2 + 
	  h*(F(ix,iy).dzdy*q1 + F(ix,iy+1).dzdy*q2));
}


inline double BicubicSpline::Deriv (int ix, double y)
{
  if (!YUpToDate(ix))
    YUpdate(ix);
  
  int iy = Ygrid->ReverseMap(y);
  iy = max(0,iy);
  iy = min(iy, Ny-2);
  
  double t = (y - (*Ygrid)(iy))/((*Ygrid)(iy+1) - (*Ygrid)(iy));
  double tm1 = t - 1.0;
  double h = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double hinv = 1.0/h;
  double dp1 = 6.0*t*tm1;  
  double dq1 = tm1*(3.0*t-1.0);
  double dp2 = -dp1;  
  double dq2 = 3.0*t*t - 2.0*t;

  return (hinv*(F(ix,iy).z*dp1 + F(ix,iy+1).z*dp2) + 
	  (F(ix,iy).dzdy*dq1 + F(ix,iy+1).dzdy*dq2));
}

inline double BicubicSpline::Deriv2 (int ix, double y)
{
  if (!YUpToDate(ix))
    YUpdate(ix);
  
  int iy = Ygrid->ReverseMap(y);
  iy = max(0,iy);
  iy = min(iy, Ny-2);
  double t = (y - (*Ygrid)(iy))/((*Ygrid)(iy+1) - (*Ygrid)(iy));
  double h = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double hinv = 1.0/h;

  return (hinv*((12.0*t-6.0)*hinv*(F(ix,iy).z -F(ix,iy+1).z)
		+(6.0*t-4.0)*F(ix,iy).dzdy + (6.0*t-2.0)*F(ix,iy+1).dzdy));
}

inline double BicubicSpline::Deriv3 (int ix, double y)
{

  if (!YUpToDate(ix))
    YUpdate(ix);
  
  int iy = Ygrid->ReverseMap(y);
  iy = max(0,iy);
  iy = min(iy, Ny-2);
  double h = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double hinv = 1.0/h;

  return (hinv*hinv*(12.0*hinv*(F(ix,iy).z -F(ix,iy+1).z)
		+6.0*(F(ix,iy).dzdy + F(ix,iy+1).dzdy)));
 }


inline double p1(double t)
{ return ((t-1.0)*(t-1.0)*(1.0+2.0*t)); }

inline double p2(double t)
{ return (t*t*(3.0-2.0*t)); }

inline double q1(double t)
{ return (t*(t-1.0)*(t-1.0)); }

inline double q2(double t)
{ return (t*t*(t-1.0)); }

inline double dp1(double t)
{ return (6.0*t*(t-1.0)); }

inline double dq1(double t)
{ return ((t-1.0)*(3.0*t-1.0)); }

inline double dp2(double t)
{ return (-dp1(t)); }

inline double dq2 (double t)
{ return (3.0*t*t - 2.0*t); }

inline double d2p1(double t)
{ return (12.0*t-6.0); }

inline double d2q1 (double t)
{ return (6.0*t - 4.0); }

inline double d2p2 (double t)
{ return (-d2p1(t)); }

inline double d2q2 (double t)
{ return (6.0*t - 2.0); } 

// inline double BicubicSpline::operator() (double x, double y)
// {
//   if (!BiUpToDate)
//     BiUpdate();
//   TinyMatrix<double,4,4> Z;
//   TinyVector<double,4> a, b;

//   int ix = Xgrid->ReverseMap(x);  
//   int iy = Ygrid->ReverseMap(y);

//   ix = max(0,ix); ix = min(ix, Nx-2);
//   iy = max(0,iy); iy = min(iy, Ny-2);

//   double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
//   double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
//   double u = (x - (*Xgrid)(ix))/h;
//   double v = (y - (*Ygrid)(iy))/k;
//   a(0) = p1(u);
//   a(1) = p2(u);
//   a(2) = h*q1(u);
//   a(3) = h*q2(u);

//   b(0) = p1(v);
//   b(1) = p2(v);
//   b(2) = k*q1(v);
//   b(3) = k*q2(v);
  
//   Z(0,0) = F(ix,iy).z;
//   Z(0,1) = F(ix,iy+1).z;
//   Z(0,2) = F(ix,iy).dzdy;
//   Z(0,3) = F(ix,iy+1).dzdy;
//   Z(1,0) = F(ix+1,iy).z;
//   Z(1,1) = F(ix+1,iy+1).z;
//   Z(1,2) = F(ix+1,iy).dzdy;
//   Z(1,3) = F(ix+1,iy+1).dzdy;
//   Z(2,0) = F(ix,iy).dzdx;
//   Z(2,1) = F(ix,iy+1).dzdx;
//   Z(2,2) = F(ix,iy).d2zdxdy;
//   Z(2,3) = F(ix,iy+1).d2zdxdy;
//   Z(3,0) = F(ix+1,iy).dzdx;
//   Z(3,1) = F(ix+1,iy+1).dzdx;
//   Z(3,2) = F(ix+1,iy).d2zdxdy;
//   Z(3,3) = F(ix+1,iy+1).d2zdxdy;
  
//   double val = 0.0;
//   for (int m=0; m<4; m++) {
//     double Zb_m = 0.0;
//     for (int n=0; n<4; n++)
//       Zb_m += Z(m,n) * b(n);
//     val += Zb_m * a(m);
//   }
//   return (val);
// }

inline double BicubicSpline::operator() (double x, double y)
{
  if (!BiUpToDate)
    BiUpdate();
  //  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

//   if (x > (Xgrid->End*1.00001)) 
//     cerr << "x too large in BiCubicSpline:  " << x << endl;
//   if (y > (Ygrid->End*1.0001)) 
//     cerr << "y too large in BiCubicSpline:  " << y << endl;


  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  a(0) = p1(u);
  a(1) = p2(u);
  a(2) = h*q1(u);
  a(3) = h*q2(u);

  b(0) = p1(v);
  b(1) = p2(v);
  b(2) = k*q1(v);
  b(3) = k*q2(v);
  
  double &Z00 = F(ix,iy).z;
  double &Z01 = F(ix,iy+1).z;
  double &Z02 = F(ix,iy).dzdy;
  double &Z03 = F(ix,iy+1).dzdy;
  double &Z10 = F(ix+1,iy).z;
  double &Z11 = F(ix+1,iy+1).z;
  double &Z12 = F(ix+1,iy).dzdy;
  double &Z13 = F(ix+1,iy+1).dzdy;
  double &Z20 = F(ix,iy).dzdx;
  double &Z21 = F(ix,iy+1).dzdx;
  double &Z22 = F(ix,iy).d2zdxdy;
  double &Z23 = F(ix,iy+1).d2zdxdy;
  double &Z30 = F(ix+1,iy).dzdx;
  double &Z31 = F(ix+1,iy+1).dzdx;
  double &Z32 = F(ix+1,iy).d2zdxdy;
  double &Z33 = F(ix+1,iy+1).d2zdxdy;
  
  double val = 
      a(0)*(Z00*b(0)+Z01*b(1)+Z02*b(2)+ Z03*b(3)) +
      a(1)*(Z10*b(0)+Z11*b(1)+Z12*b(2)+ Z13*b(3)) +
      a(2)*(Z20*b(0)+Z21*b(1)+Z22*b(2)+ Z23*b(3)) +
      a(3)*(Z30*b(0)+Z31*b(1)+Z32*b(2)+ Z33*b(3));

  return (val);
}



inline double BicubicSpline::d_dx (double x, double y)
{
  if (!BiUpToDate)
    BiUpdate();
  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double u = (x - (*Xgrid)(ix))*hinv;
  double v = (y - (*Ygrid)(iy))/k;
  a(0) = hinv*dp1(u);
  a(1) = hinv*dp2(u);
  a(2) = dq1(u);
  a(3) = dq2(u);

  b(0) = p1(v);
  b(1) = p2(v);
  b(2) = k*q1(v);
  b(3) = k*q2(v);
  
  Z(0,0) = F(ix,iy).z;
  Z(0,1) = F(ix,iy+1).z;
  Z(0,2) = F(ix,iy).dzdy;
  Z(0,3) = F(ix,iy+1).dzdy;
  Z(1,0) = F(ix+1,iy).z;
  Z(1,1) = F(ix+1,iy+1).z;
  Z(1,2) = F(ix+1,iy).dzdy;
  Z(1,3) = F(ix+1,iy+1).dzdy;
  Z(2,0) = F(ix,iy).dzdx;
  Z(2,1) = F(ix,iy+1).dzdx;
  Z(2,2) = F(ix,iy).d2zdxdy;
  Z(2,3) = F(ix,iy+1).d2zdxdy;
  Z(3,0) = F(ix+1,iy).dzdx;
  Z(3,1) = F(ix+1,iy+1).dzdx;
  Z(3,2) = F(ix+1,iy).d2zdxdy;
  Z(3,3) = F(ix+1,iy+1).d2zdxdy;
  
  double val = 0.0;
  for (int m=0; m<4; m++) {
    double Zb_m = 0.0;
    for (int n=0; n<4; n++)
      Zb_m += Z(m,n) * b(n);
    val += Zb_m * a(m);
  }
  return (val);
}

inline double BicubicSpline::d_dy (double x, double y)
{
  if (!BiUpToDate)
    BiUpdate();
  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double hinv = 1.0/h;
  double kinv = 1.0/k;
  double u = (x - (*Xgrid)(ix))*hinv;
  double v = (y - (*Ygrid)(iy))*kinv;
  a(0) = p1(u);
  a(1) = p2(u);
  a(2) = h*q1(u);
  a(3) = h*q2(u);

  b(0) = kinv*dp1(v);
  b(1) = kinv*dp2(v);
  b(2) = dq1(v);
  b(3) = dq2(v);
  
  Z(0,0) = F(ix,iy).z;
  Z(0,1) = F(ix,iy+1).z;
  Z(0,2) = F(ix,iy).dzdy;
  Z(0,3) = F(ix,iy+1).dzdy;
  Z(1,0) = F(ix+1,iy).z;
  Z(1,1) = F(ix+1,iy+1).z;
  Z(1,2) = F(ix+1,iy).dzdy;
  Z(1,3) = F(ix+1,iy+1).dzdy;
  Z(2,0) = F(ix,iy).dzdx;
  Z(2,1) = F(ix,iy+1).dzdx;
  Z(2,2) = F(ix,iy).d2zdxdy;
  Z(2,3) = F(ix,iy+1).d2zdxdy;
  Z(3,0) = F(ix+1,iy).dzdx;
  Z(3,1) = F(ix+1,iy+1).dzdx;
  Z(3,2) = F(ix+1,iy).d2zdxdy;
  Z(3,3) = F(ix+1,iy+1).d2zdxdy;
  
  double val = 0.0;
  for (int m=0; m<4; m++) {
    double Zb_m = 0.0;
    for (int n=0; n<4; n++)
      Zb_m += Z(m,n) * b(n);
    val += Zb_m * a(m);
  }
  return (val);
}

inline double BicubicSpline::d2_dxdy (double x, double y)
{
  if (!BiUpToDate)
    BiUpdate();
  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double hinv = 1.0/h;
  double kinv = 1.0/k;
  double u = (x - (*Xgrid)(ix))*hinv;
  double v = (y - (*Ygrid)(iy))*kinv;
  a(0) = hinv*dp1(u);
  a(1) = hinv*dp2(u);
  a(2) = dq1(u);
  a(3) = dq2(u);

  b(0) = kinv*dp1(v);
  b(1) = kinv*dp2(v);
  b(2) = dq1(v);
  b(3) = dq2(v);
  
  Z(0,0) = F(ix,iy).z;
  Z(0,1) = F(ix,iy+1).z;
  Z(0,2) = F(ix,iy).dzdy;
  Z(0,3) = F(ix,iy+1).dzdy;
  Z(1,0) = F(ix+1,iy).z;
  Z(1,1) = F(ix+1,iy+1).z;
  Z(1,2) = F(ix+1,iy).dzdy;
  Z(1,3) = F(ix+1,iy+1).dzdy;
  Z(2,0) = F(ix,iy).dzdx;
  Z(2,1) = F(ix,iy+1).dzdx;
  Z(2,2) = F(ix,iy).d2zdxdy;
  Z(2,3) = F(ix,iy+1).d2zdxdy;
  Z(3,0) = F(ix+1,iy).dzdx;
  Z(3,1) = F(ix+1,iy+1).dzdx;
  Z(3,2) = F(ix+1,iy).d2zdxdy;
  Z(3,3) = F(ix+1,iy+1).d2zdxdy;
  
  double val = 0.0;
  for (int m=0; m<4; m++) {
    double Zb_m = 0.0;
    for (int n=0; n<4; n++)
      Zb_m += Z(m,n) * b(n);
    val += Zb_m * a(m);
  }
  return (val);
}

inline double BicubicSpline::d2_dx2 (double x, double y)
{
  if (!BiUpToDate)
    BiUpdate();
  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double u = (x - (*Xgrid)(ix))*hinv;
  double v = (y - (*Ygrid)(iy))/k;
  a(0) = hinv*hinv*d2p1(u);
  a(1) = hinv*hinv*d2p2(u);
  a(2) = hinv*d2q1(u);
  a(3) = hinv*d2q2(u);

  b(0) = p1(v);
  b(1) = p2(v);
  b(2) = k*q1(v);
  b(3) = k*q2(v);
  
  Z(0,0) = F(ix,iy).z;
  Z(0,1) = F(ix,iy+1).z;
  Z(0,2) = F(ix,iy).dzdy;
  Z(0,3) = F(ix,iy+1).dzdy;
  Z(1,0) = F(ix+1,iy).z;
  Z(1,1) = F(ix+1,iy+1).z;
  Z(1,2) = F(ix+1,iy).dzdy;
  Z(1,3) = F(ix+1,iy+1).dzdy;
  Z(2,0) = F(ix,iy).dzdx;
  Z(2,1) = F(ix,iy+1).dzdx;
  Z(2,2) = F(ix,iy).d2zdxdy;
  Z(2,3) = F(ix,iy+1).d2zdxdy;
  Z(3,0) = F(ix+1,iy).dzdx;
  Z(3,1) = F(ix+1,iy+1).dzdx;
  Z(3,2) = F(ix+1,iy).d2zdxdy;
  Z(3,3) = F(ix+1,iy+1).d2zdxdy;
  
  double val = 0.0;
  for (int m=0; m<4; m++) {
    double Zb_m = 0.0;
    for (int n=0; n<4; n++)
      Zb_m += Z(m,n) * b(n);
    val += Zb_m * a(m);
  }
  return (val);
}




inline double BicubicSpline::d2_dy2 (double x, double y)
{
  if (!BiUpToDate)
    BiUpdate();
  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double u = (x - (*Xgrid)(ix))*hinv;
  double v = (y - (*Ygrid)(iy))*kinv;
  a(0) = p1(u);
  a(1) = p2(u);
  a(2) = h*q1(u);
  a(3) = h*q2(u);

  b(0) = kinv*kinv*d2p1(v);
  b(1) = kinv*kinv*d2p2(v);
  b(2) = kinv*d2q1(v);
  b(3) = kinv*d2q2(v);
  
  Z(0,0) = F(ix,iy).z;
  Z(0,1) = F(ix,iy+1).z;
  Z(0,2) = F(ix,iy).dzdy;
  Z(0,3) = F(ix,iy+1).dzdy;
  Z(1,0) = F(ix+1,iy).z;
  Z(1,1) = F(ix+1,iy+1).z;
  Z(1,2) = F(ix+1,iy).dzdy;
  Z(1,3) = F(ix+1,iy+1).dzdy;
  Z(2,0) = F(ix,iy).dzdx;
  Z(2,1) = F(ix,iy+1).dzdx;
  Z(2,2) = F(ix,iy).d2zdxdy;
  Z(2,3) = F(ix,iy+1).d2zdxdy;
  Z(3,0) = F(ix+1,iy).dzdx;
  Z(3,1) = F(ix+1,iy+1).dzdx;
  Z(3,2) = F(ix+1,iy).d2zdxdy;
  Z(3,3) = F(ix+1,iy+1).d2zdxdy;
  
  double val = 0.0;
  for (int m=0; m<4; m++) {
    double Zb_m = 0.0;
    for (int n=0; n<4; n++)
      Zb_m += Z(m,n) * b(n);
    val += Zb_m * a(m);
  }
  return (val);
}



//////////////////////////////////////////////////////////////////
// THIS IS HIGHLY BROKEN!!!!!!!!!!!!!!!!!!                      //
// DO NOT USE UNDER ANY CIRCUMSTANCES!!!!!!!!!!!!!!!!!!!!!      //
//////////////////////////////////////////////////////////////////
/// This class is used to interpolate a function z(x,y) in two
/// dimensions.  It is of cubic order in both directions.  It
/// currently uses only "natural" boundary conditions:  all second
/// derivatives are 0 at the boundary.
class SymmBicubicSpline
{
private:
  /// Stores whether dzdx needs to be recalculated for each column
  Array<bool,1> XUpToDate;
  /// Stores whether dzdy needs to be recalculated for each row
  Array<bool,1> YUpToDate;
  /// Stores whether the d2zdxdy's need to be updated
  bool BiUpToDate;
  /// Update the iy column's dzdx's 
  void XUpdate(int iy);
  /// Update the ix row's dzdy's
  void YUpdate(int ix);
  /// Update all the d2zdxdy's
  void BiUpdate();
  /// Holds the BCpoint array containing z and its derivatives.
  SymmArray<BCpoint> F;
public:
  int Nx, Ny;
  Grid *Xgrid, *Ygrid;

  inline double z(int ix, int iy) const
  { return (F(ix,iy).z); }
  inline double &z(int ix, int iy)
  { return (F(ix,iy).z); }
  inline double dx(int ix, int iy) const
  { return ((iy>ix) ? F(iy,ix).dzdy : F(ix,iy).dzdx); }
  inline double& dx(int ix, int iy) 
  { return ((iy>ix) ? F(iy,ix).dzdy : F(ix,iy).dzdx); }
  inline double dy(int ix, int iy) const
  { return ((iy>ix) ? F(iy,ix).dzdx : F(ix,iy).dzdy); }
  inline double &dy(int ix, int iy) 
  { return ((iy>ix) ? F(iy,ix).dzdx : F(ix,iy).dzdy); }
  inline double dxdy (int ix, int iy) const
  { return (F(ix,iy).d2zdxdy); }
  inline double& dxdy (int ix, int iy)
  { return (F(ix,iy).d2zdxdy); }

  /// Interpolates in y at a given row
  inline double  operator() (int ix,   double y);
  /// Interpolates in x at a given column
  inline double  operator() (double x, int iy);
  /// Bicubic interpolation in x and y
  inline double  operator() (double x, double y);
  /// Returns z(x_{ix}, y_{iy})
  inline double  operator() (int ix, int iy) const;
  /// Returns a reference to z(x_{ix}, y_{iy}) that can be used as an
  /// L-value.
  inline double& operator() (int ix, int iy);
  inline double  Deriv      (int ix,   double y);
  inline double  Deriv      (double x, int iy);
  inline double  xDeriv     (int ix, int iy);
  inline double  yDeriv     (int ix, int iy);
  inline double  Deriv2     (int ix,   double y);
  inline double  Deriv2     (double x, int iy);
  inline double  Deriv3     (int ix,   double y);
  inline double  Deriv3     (double x, int iy);
  inline double  d_dx       (double x, double y);
  inline double  d_dy       (double x, double y);
  inline double  d2_dx2     (double x, double y);
  inline double  d2_dy2     (double x, double y);
  inline double  d2_dxdy    (double x, double y);
  /// Initialize the bicubic spline with the given grids and data
  inline void Init (Grid *xgrid, Grid *ygrid, Array<double,2> &f);

  /// Copy constructor
  inline SymmBicubicSpline (const SymmBicubicSpline &a);

  /// Assigment operator -- necessary for array resizeAndPreserve
  inline SymmBicubicSpline & operator= (SymmBicubicSpline &a);
  inline SymmBicubicSpline & operator= (SymmBicubicSpline a);
  inline SymmBicubicSpline() {}
};

inline SymmBicubicSpline::SymmBicubicSpline (const SymmBicubicSpline &a)
{
  XUpToDate.resize(a.XUpToDate.size());
  XUpToDate = a.XUpToDate;
  YUpToDate.resize(a.YUpToDate.size());
  YUpToDate = a.YUpToDate;
  F.resize(a.F.rows());
  F=a.F;
  Nx=a.Nx; Ny=a.Ny;
  Xgrid=a.Xgrid; Ygrid=a.Ygrid;
  BiUpToDate = a.BiUpToDate;
}

inline SymmBicubicSpline& SymmBicubicSpline::operator=(SymmBicubicSpline &a)
{
  XUpToDate.resize(a.XUpToDate.size());
  XUpToDate = a.XUpToDate;
  YUpToDate.resize(a.YUpToDate.size());
  YUpToDate = a.YUpToDate;
  F.resize(a.F.rows());
  F=a.F;
  Nx=a.Nx; Ny=a.Ny;
  Xgrid=a.Xgrid; Ygrid=a.Ygrid;
  BiUpToDate = a.BiUpToDate;
  return *this;
}

inline SymmBicubicSpline& SymmBicubicSpline::operator=(SymmBicubicSpline a)
{
  XUpToDate.resize(a.XUpToDate.size());
  XUpToDate = a.XUpToDate;
  YUpToDate.resize(a.YUpToDate.size());
  YUpToDate = a.YUpToDate;
  F.resize(a.F.rows());
  F=a.F;
  Nx=a.Nx; Ny=a.Ny;
  Xgrid=a.Xgrid; Ygrid=a.Ygrid;
  BiUpToDate = a.BiUpToDate;
  return *this;
}


inline void SymmBicubicSpline::Init(Grid *xgrid, Grid *ygrid, Array<double,2> &f)
{
  Nx = xgrid->NumPoints;
  Ny = ygrid->NumPoints;
  Xgrid = xgrid; Ygrid = ygrid;
  assert (f.rows() == Nx);
  assert (f.cols() == Ny);
  XUpToDate.resize(Ny);
  YUpToDate.resize(Nx);
  F.resize(f.rows());
  for (int i=0; i<Nx; i++)
    for (int j=0; j<Ny; j++)
      z(i,j) = f(i,j);
  for (int i=0; i<Ny; i++)
    XUpdate(i);
  for (int i=0; i<Nx; i++)
    YUpdate(i);

  BiUpdate();
}

inline double SymmBicubicSpline::operator() (int ix, int iy) const
{ return (z(ix,iy));}

inline double& SymmBicubicSpline::operator() (int ix, int iy)
{ 
  XUpToDate(iy) = false;
  YUpToDate(ix) = false;
  BiUpToDate = false;
  return (z(ix,iy));
}

inline double SymmBicubicSpline::operator() (double x, int iy)
{
  if (!XUpToDate(iy))
    XUpdate(iy);
  
  int ix = Xgrid->ReverseMap(x);
  ix = max(0,ix);
  ix = min(ix, Nx-2);
  
  double t = (x - (*Xgrid)(ix))/((*Xgrid)(ix+1) - (*Xgrid)(ix));
  double tm1 = t - 1.0;
  double p1 = tm1*tm1*(1.0+2.0*t);
  double q1 = t*tm1*tm1;
  double p2 = t*t*(3.0-2.0*t);
  double q2 = t*t*tm1;
  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);

  return (z(ix,iy)*p1 + z(ix+1,iy)*p2 + 
	  h*(dx(ix,iy)*q1 + dx(ix+1,iy)*q2));
}


inline double SymmBicubicSpline::operator() (int ix, double y)
{
  if (!YUpToDate(ix))
    YUpdate(ix);
  
  int iy = Ygrid->ReverseMap(y);
  iy = max(0,iy);
  iy = min(iy, Ny-2);
  
  double t = (y - (*Ygrid)(iy))/((*Ygrid)(iy+1) - (*Ygrid)(iy));
  double tm1 = t - 1.0;
  double p1 = tm1*tm1*(1.0+2.0*t);
  double q1 = t*tm1*tm1;
  double p2 = t*t*(3.0-2.0*t);
  double q2 = t*t*tm1;
  double h = (*Ygrid)(iy+1) - (*Ygrid)(iy);

  return (z(ix,iy)*p1 + z(ix,iy+1)*p2 + 
	  h*(dy(ix,iy)*q1 + dy(ix,iy+1)*q2));
}


inline double SymmBicubicSpline::Deriv (int ix, double y)
{
  if (!YUpToDate(ix))
    YUpdate(ix);
  
  int iy = Ygrid->ReverseMap(y);
  iy = max(0,iy);
  iy = min(iy, Ny-2);
  
  double t = (y - (*Ygrid)(iy))/((*Ygrid)(iy+1) - (*Ygrid)(iy));
  double tm1 = t - 1.0;
  double h = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double hinv = 1.0/h;
  double dp1 = 6.0*t*tm1;  
  double dq1 = tm1*(3.0*t-1.0);
  double dp2 = -dp1;  
  double dq2 = 3.0*t*t - 2.0*t;

  return (hinv*(z(ix,iy)*dp1 + z(ix,iy+1)*dp2) + 
	  (dy(ix,iy)*dq1 + dy(ix,iy+1)*dq2));
}

inline double SymmBicubicSpline::Deriv2 (int ix, double y)
{
  if (!YUpToDate(ix))
    YUpdate(ix);
  
  int iy = Ygrid->ReverseMap(y);
  iy = max(0,iy);
  iy = min(iy, Ny-2);
  double t = (y - (*Ygrid)(iy))/((*Ygrid)(iy+1) - (*Ygrid)(iy));
  double h = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double hinv = 1.0/h;

  return (hinv*((12.0*t-6.0)*hinv*(z(ix,iy) -z(ix,iy+1))
		+(6.0*t-4.0)*dy(ix,iy) + (6.0*t-2.0)*dy(ix,iy+1)));
}

inline double SymmBicubicSpline::Deriv3 (int ix, double y)
{

  if (!YUpToDate(ix))
    YUpdate(ix);
  
  int iy = Ygrid->ReverseMap(y);
  iy = max(0,iy);
  iy = min(iy, Ny-2);
  double h = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double hinv = 1.0/h;

  return (hinv*hinv*(12.0*hinv*(z(ix,iy) -z(ix,iy+1))
		+6.0*(dy(ix,iy) + dy(ix,iy+1))));
}



inline double SymmBicubicSpline::operator() (double x, double y)
{
  if (y < x) {
    double t = y;
    y = x; 
    x = t;
  }

  if (!BiUpToDate)
    BiUpdate();
  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  a(0) = p1(u);
  a(1) = p2(u);
  a(2) = h*q1(u);
  a(3) = h*q2(u);

  b(0) = p1(v);
  b(1) = p2(v);
  b(2) = k*q1(v);
  b(3) = k*q2(v);
  
  Z(0,0) = z(ix,iy);
  Z(0,1) = z(ix,iy+1);
  Z(0,2) = dy(ix,iy);
  Z(0,3) = dy(ix,iy+1);
  Z(1,0) = z(ix+1,iy);
  Z(1,1) = z(ix+1,iy+1);
  Z(1,2) = dy(ix+1,iy);
  Z(1,3) = dy(ix+1,iy+1);
  Z(2,0) = dx(ix,iy);
  Z(2,1) = dx(ix,iy+1);
  Z(2,2) = dxdy(ix,iy);
  Z(2,3) = dxdy(ix,iy+1);
  Z(3,0) = dx(ix+1,iy);
  Z(3,1) = dx(ix+1,iy+1);
  Z(3,2) = dxdy(ix+1,iy);
  Z(3,3) = dxdy(ix+1,iy+1);
  
  double val = 0.0;
  for (int m=0; m<4; m++) {
    double Zb_m = 0.0;
    for (int n=0; n<4; n++)
      Zb_m += Z(m,n) * b(n);
    val += Zb_m * a(m);
  }
  return (val);
}

inline double SymmBicubicSpline::d_dx (double x, double y)
{
  if (!BiUpToDate)
    BiUpdate();
  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double u = (x - (*Xgrid)(ix))*hinv;
  double v = (y - (*Ygrid)(iy))/k;
  a(0) = hinv*dp1(u);
  a(1) = hinv*dp2(u);
  a(2) = dq1(u);
  a(3) = dq2(u);

  b(0) = p1(v);
  b(1) = p2(v);
  b(2) = k*q1(v);
  b(3) = k*q2(v);
  
  Z(0,0) = z(ix,iy);
  Z(0,1) = z(ix,iy+1);
  Z(0,2) = dy(ix,iy);
  Z(0,3) = dy(ix,iy+1);
  Z(1,0) = z(ix+1,iy);
  Z(1,1) = z(ix+1,iy+1);
  Z(1,2) = dy(ix+1,iy);
  Z(1,3) = dy(ix+1,iy+1);
  Z(2,0) = dx(ix,iy);
  Z(2,1) = dx(ix,iy+1);
  Z(2,2) = dxdy(ix,iy);
  Z(2,3) = dxdy(ix,iy+1);
  Z(3,0) = dx(ix+1,iy);
  Z(3,1) = dx(ix+1,iy+1);
  Z(3,2) = dxdy(ix+1,iy);
  Z(3,3) = dxdy(ix+1,iy+1);
  
  double val = 0.0;
  for (int m=0; m<4; m++) {
    double Zb_m = 0.0;
    for (int n=0; n<4; n++)
      Zb_m += Z(m,n) * b(n);
    val += Zb_m * a(m);
  }
  return (val);
}

inline double SymmBicubicSpline::d_dy (double x, double y)
{
  if (!BiUpToDate)
    BiUpdate();
  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double hinv = 1.0/h;
  double kinv = 1.0/k;
  double u = (x - (*Xgrid)(ix))*hinv;
  double v = (y - (*Ygrid)(iy))*kinv;
  a(0) = p1(u);
  a(1) = p2(u);
  a(2) = h*q1(u);
  a(3) = h*q2(u);

  b(0) = kinv*dp1(v);
  b(1) = kinv*dp2(v);
  b(2) = dq1(v);
  b(3) = dq2(v);
  
  Z(0,0) = z(ix,iy);
  Z(0,1) = z(ix,iy+1);
  Z(0,2) = dy(ix,iy);
  Z(0,3) = dy(ix,iy+1);
  Z(1,0) = z(ix+1,iy);
  Z(1,1) = z(ix+1,iy+1);
  Z(1,2) = dy(ix+1,iy);
  Z(1,3) = dy(ix+1,iy+1);
  Z(2,0) = dx(ix,iy);
  Z(2,1) = dx(ix,iy+1);
  Z(2,2) = dxdy(ix,iy);
  Z(2,3) = dxdy(ix,iy+1);
  Z(3,0) = dx(ix+1,iy);
  Z(3,1) = dx(ix+1,iy+1);
  Z(3,2) = dxdy(ix+1,iy);
  Z(3,3) = dxdy(ix+1,iy+1);
  
  double val = 0.0;
  for (int m=0; m<4; m++) {
    double Zb_m = 0.0;
    for (int n=0; n<4; n++)
      Zb_m += Z(m,n) * b(n);
    val += Zb_m * a(m);
  }
  return (val);
}

inline double SymmBicubicSpline::d2_dxdy (double x, double y)
{
  if (!BiUpToDate)
    BiUpdate();
  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double hinv = 1.0/h;
  double kinv = 1.0/k;
  double u = (x - (*Xgrid)(ix))*hinv;
  double v = (y - (*Ygrid)(iy))*kinv;
  a(0) = hinv*dp1(u);
  a(1) = hinv*dp2(u);
  a(2) = dq1(u);
  a(3) = dq2(u);

  b(0) = kinv*dp1(v);
  b(1) = kinv*dp2(v);
  b(2) = dq1(v);
  b(3) = dq2(v);
  
  Z(0,0) = z(ix,iy);
  Z(0,1) = z(ix,iy+1);
  Z(0,2) = dy(ix,iy);
  Z(0,3) = dy(ix,iy+1);
  Z(1,0) = z(ix+1,iy);
  Z(1,1) = z(ix+1,iy+1);
  Z(1,2) = dy(ix+1,iy);
  Z(1,3) = dy(ix+1,iy+1);
  Z(2,0) = dx(ix,iy);
  Z(2,1) = dx(ix,iy+1);
  Z(2,2) = dxdy(ix,iy);
  Z(2,3) = dxdy(ix,iy+1);
  Z(3,0) = dx(ix+1,iy);
  Z(3,1) = dx(ix+1,iy+1);
  Z(3,2) = dxdy(ix+1,iy);
  Z(3,3) = dxdy(ix+1,iy+1);
  
  double val = 0.0;
  for (int m=0; m<4; m++) {
    double Zb_m = 0.0;
    for (int n=0; n<4; n++)
      Zb_m += Z(m,n) * b(n);
    val += Zb_m * a(m);
  }
  return (val);
}

inline double SymmBicubicSpline::d2_dx2 (double x, double y)
{
  if (!BiUpToDate)
    BiUpdate();
  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double u = (x - (*Xgrid)(ix))*hinv;
  double v = (y - (*Ygrid)(iy))/k;
  a(0) = hinv*hinv*d2p1(u);
  a(1) = hinv*hinv*d2p2(u);
  a(2) = hinv*d2q1(u);
  a(3) = hinv*d2q2(u);

  b(0) = p1(v);
  b(1) = p2(v);
  b(2) = k*q1(v);
  b(3) = k*q2(v);
  
  Z(0,0) = z(ix,iy);
  Z(0,1) = z(ix,iy+1);
  Z(0,2) = dy(ix,iy);
  Z(0,3) = dy(ix,iy+1);
  Z(1,0) = z(ix+1,iy);
  Z(1,1) = z(ix+1,iy+1);
  Z(1,2) = dy(ix+1,iy);
  Z(1,3) = dy(ix+1,iy+1);
  Z(2,0) = dx(ix,iy);
  Z(2,1) = dx(ix,iy+1);
  Z(2,2) = dxdy(ix,iy);
  Z(2,3) = dxdy(ix,iy+1);
  Z(3,0) = dx(ix+1,iy);
  Z(3,1) = dx(ix+1,iy+1);
  Z(3,2) = dxdy(ix+1,iy);
  Z(3,3) = dxdy(ix+1,iy+1);
  
  double val = 0.0;
  for (int m=0; m<4; m++) {
    double Zb_m = 0.0;
    for (int n=0; n<4; n++)
      Zb_m += Z(m,n) * b(n);
    val += Zb_m * a(m);
  }
  return (val);
}




inline double SymmBicubicSpline::d2_dy2 (double x, double y)
{
  if (!BiUpToDate)
    BiUpdate();
  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double u = (x - (*Xgrid)(ix))*hinv;
  double v = (y - (*Ygrid)(iy))*kinv;
  a(0) = p1(u);
  a(1) = p2(u);
  a(2) = h*q1(u);
  a(3) = h*q2(u);

  b(0) = kinv*kinv*d2p1(v);
  b(1) = kinv*kinv*d2p2(v);
  b(2) = kinv*d2q1(v);
  b(3) = kinv*d2q2(v);
  
  Z(0,0) = z(ix,iy);
  Z(0,1) = z(ix,iy+1);
  Z(0,2) = dy(ix,iy);
  Z(0,3) = dy(ix,iy+1);
  Z(1,0) = z(ix+1,iy);
  Z(1,1) = z(ix+1,iy+1);
  Z(1,2) = dy(ix+1,iy);
  Z(1,3) = dy(ix+1,iy+1);
  Z(2,0) = dx(ix,iy);
  Z(2,1) = dx(ix,iy+1);
  Z(2,2) = dxdy(ix,iy);
  Z(2,3) = dxdy(ix,iy+1);
  Z(3,0) = dx(ix+1,iy);
  Z(3,1) = dx(ix+1,iy+1);
  Z(3,2) = dxdy(ix+1,iy);
  Z(3,3) = dxdy(ix+1,iy+1);
  
  double val = 0.0;
  for (int m=0; m<4; m++) {
    double Zb_m = 0.0;
    for (int n=0; n<4; n++)
      Zb_m += Z(m,n) * b(n);
    val += Zb_m * a(m);
  }
  return (val);
}
 



class MultiBicubicSpline
{
private:
  Array<bool,1> XUpToDate;
  Array<bool,1> YUpToDate;
  bool BiUpToDate;
  void XUpdate(int ix);
  void YUpdate(int iy);
  void BiUpdate();
  Array<BCpoint,3> F;
public:
  int Nx, Ny, Nz;
  Grid *Xgrid, *Ygrid;

  inline double  operator() (int ix,   double y, int iz);
  inline double  operator() (double x, int iy, int iz);
  inline double  operator() (double x, double y, int iz);
  inline double  operator() (int ix, int iy, int iz) const;
  inline double& operator() (int ix, int iy, int iz);
  inline void    operator() (double x, double y, Array<double,1> &z);
  inline double  Deriv      (int ix,   double y, int iz);
  inline double  Deriv      (double x, int iy, int iz);
  inline double  xDeriv     (int ix, int iy, int iz);
  inline double  yDeriv     (int ix, int iy, int iz);
  inline double  Deriv2     (int ix,   double y, int iz);
  inline double  Deriv2     (double x, int iy, int iz);
  inline double  Deriv3     (int ix,   double y, int iz);
  inline double  Deriv3     (double x, int iy, int iz);
  inline void Init (Grid *xgrid, Grid *ygrid, Array<double,3> &f);

  /// Copy constructor
  inline MultiBicubicSpline (const MultiBicubicSpline &a);

  /// Assigment operator -- necessary for array resizeAndPreserve
  inline MultiBicubicSpline & operator= (MultiBicubicSpline &a);
  inline MultiBicubicSpline & operator= (MultiBicubicSpline a);
  inline MultiBicubicSpline() {}

};

inline void 
MultiBicubicSpline::Init(Grid *xgrid, Grid *ygrid, Array<double,3> &f)
{
  Nx = xgrid->NumPoints;
  Ny = ygrid->NumPoints;
  Nz = f.extent(2);
  Xgrid = xgrid; Ygrid = ygrid;
  assert (f.rows() == Nx);
  assert (f.cols() == Ny);
  XUpToDate.resize(Ny);
  YUpToDate.resize(Nx);
  F.resize(Nx, Ny, Nz);
  for (int i=0; i<Nx; i++)
    for (int j=0; j<Ny; j++)
      for (int k=0; k<Nz; k++)
	F(i,j,k).z = f(i,j,k);
  for (int i=0; i<Ny; i++)
    XUpdate(i);
  for (int i=0; i<Nx; i++)
    YUpdate(i);

  BiUpdate();
}


inline double MultiBicubicSpline::operator() (int ix, int iy, int iz) const
{ return (F(ix,iy,iz).z);}

inline double& MultiBicubicSpline::operator() (int ix, int iy, int iz)
{ 
  XUpToDate(iy) = false;
  YUpToDate(ix) = false;
  BiUpToDate = false;
  return (F(ix,iy,iz).z);
}


inline double MultiBicubicSpline::operator() (double x, int iy, int iz)
{
  if (!XUpToDate(iy))
    XUpdate(iy);
  
  int ix = Xgrid->ReverseMap(x);
  ix = max(0,ix);
  ix = min(ix, Nx-2);
  
  double t = (x - (*Xgrid)(ix))/((*Xgrid)(ix+1) - (*Xgrid)(ix));
  double tm1 = t - 1.0;
  double p1 = tm1*tm1*(1.0+2.0*t);
  double q1 = t*tm1*tm1;
  double p2 = t*t*(3.0-2.0*t);
  double q2 = t*t*tm1;
  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);

  return (F(ix,iy,iz).z*p1 + F(ix+1,iy,iz).z*p2 + 
	  h*(F(ix,iy,iz).dzdx*q1 + F(ix+1,iy,iz).dzdx*q2));
}


inline double MultiBicubicSpline::operator() (int ix, double y, int iz)
{
  if (!YUpToDate(ix))
    YUpdate(ix);
  
  int iy = Ygrid->ReverseMap(y);
  iy = max(0,iy);
  iy = min(iy, Ny-2);
  
  double t = (y - (*Ygrid)(iy))/((*Ygrid)(iy+1) - (*Ygrid)(iy));
  double tm1 = t - 1.0;
  double p1 = tm1*tm1*(1.0+2.0*t);
  double q1 = t*tm1*tm1;
  double p2 = t*t*(3.0-2.0*t);
  double q2 = t*t*tm1;
  double h = (*Ygrid)(iy+1) - (*Ygrid)(iy);

  return (F(iy,ix,iz).z*p1 + F(iy+1,ix,iz).z*p2 + 
	  h*(F(iy,ix,iz).dzdy*q1 + F(iy+1,ix,iz).dzdy*q2));
}

inline double MultiBicubicSpline::operator() (double x, double y, int iz)
{
  if (!BiUpToDate)
    BiUpdate();
  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  a(0) = p1(u);
  a(1) = p2(u);
  a(2) = h*q1(u);
  a(3) = h*q2(u);

  b(0) = p1(v);
  b(1) = p2(v);
  b(2) = k*q1(v);
  b(3) = k*q2(v);
  
  Z(0,0) = F(ix,iy,iz).z;
  Z(0,1) = F(ix,iy+1,iz).z;
  Z(0,2) = F(ix,iy,iz).dzdy;
  Z(0,3) = F(ix,iy+1,iz).dzdy;
  Z(1,0) = F(ix+1,iy,iz).z;
  Z(1,1) = F(ix+1,iy+1,iz).z;
  Z(1,2) = F(ix+1,iy,iz).dzdy;
  Z(1,3) = F(ix+1,iy+1,iz).dzdy;
  Z(2,0) = F(ix,iy,iz).dzdx;
  Z(2,1) = F(ix,iy+1,iz).dzdx;
  Z(2,2) = F(ix,iy,iz).d2zdxdy;
  Z(2,3) = F(ix,iy+1,iz).d2zdxdy;
  Z(3,0) = F(ix+1,iy,iz).dzdx;
  Z(3,1) = F(ix+1,iy+1,iz).dzdx;
  Z(3,2) = F(ix+1,iy,iz).d2zdxdy;
  Z(3,3) = F(ix+1,iy+1,iz).d2zdxdy;
  
  double val = 0.0;
  for (int m=0; m<4; m++) {
    double Zb_m = 0.0;
    for (int n=0; n<4; n++)
      Zb_m += Z(m,n) * b(n);
    val += Zb_m * a(m);
  }
  return (val);
}



inline void MultiBicubicSpline::operator() (double x, double y, 
					    Array<double,1> &z)
{
  if (!BiUpToDate)
    BiUpdate();
  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  a(0) = p1(u);
  a(1) = p2(u);
  a(2) = h*q1(u);
  a(3) = h*q2(u);

  b(0) = p1(v);
  b(1) = p2(v);
  b(2) = k*q1(v);
  b(3) = k*q2(v);
  
  for (int iz=0; iz<Nz; iz++) {
    Z(0,0) = F(ix,iy,iz).z;
    Z(0,1) = F(ix,iy+1,iz).z;
    Z(0,2) = F(ix,iy,iz).dzdy;
    Z(0,3) = F(ix,iy+1,iz).dzdy;
    Z(1,0) = F(ix+1,iy,iz).z;
    Z(1,1) = F(ix+1,iy+1,iz).z;
    Z(1,2) = F(ix+1,iy,iz).dzdy;
    Z(1,3) = F(ix+1,iy+1,iz).dzdy;
    Z(2,0) = F(ix,iy,iz).dzdx;
    Z(2,1) = F(ix,iy+1,iz).dzdx;
    Z(2,2) = F(ix,iy,iz).d2zdxdy;
    Z(2,3) = F(ix,iy+1,iz).d2zdxdy;
    Z(3,0) = F(ix+1,iy,iz).dzdx;
    Z(3,1) = F(ix+1,iy+1,iz).dzdx;
    Z(3,2) = F(ix+1,iy,iz).d2zdxdy;
    Z(3,3) = F(ix+1,iy+1,iz).d2zdxdy;
  
    z(iz) = 0.0;
    for (int m=0; m<4; m++) {
      double Zb_m = 0.0;
      for (int n=0; n<4; n++)
	Zb_m += Z(m,n) * b(n);
      z(iz) += Zb_m * a(m);
    }
  }
}


/// Copy constructor
inline MultiBicubicSpline::MultiBicubicSpline(const MultiBicubicSpline &a)
{
  XUpToDate.resize(a.XUpToDate.size());
  XUpToDate = a.XUpToDate;
  YUpToDate.resize(a.YUpToDate.size());
  YUpToDate = a.YUpToDate;
  F.resize(a.F.extent(0), a.F.extent(1), a.F.extent(2));
  F=a.F;
  Nx=a.Nx; Ny=a.Ny; Nz=a.Nz;
  Xgrid=a.Xgrid; Ygrid=a.Ygrid;
  BiUpToDate = a.BiUpToDate;
}


inline MultiBicubicSpline& MultiBicubicSpline::operator=(MultiBicubicSpline a)
{
  XUpToDate.resize(a.XUpToDate.size());
  XUpToDate = a.XUpToDate;
  YUpToDate.resize(a.YUpToDate.size());
  YUpToDate = a.YUpToDate;
  F.resize(a.F.extent(0), a.F.extent(1), a.F.extent(2));
  F=a.F;
  Nx=a.Nx; Ny=a.Ny; Nz=a.Nz;
  Xgrid=a.Xgrid; Ygrid=a.Ygrid;
  BiUpToDate = a.BiUpToDate;
  return *this;
}

inline MultiBicubicSpline& MultiBicubicSpline::operator=(MultiBicubicSpline &a)
{
  XUpToDate.resize(a.XUpToDate.size());
  XUpToDate = a.XUpToDate;
  YUpToDate.resize(a.YUpToDate.size());
  YUpToDate = a.YUpToDate;
  F.resize(a.F.extent(0), a.F.extent(1), a.F.extent(2));
  F=a.F;
  Nx=a.Nx; Ny=a.Ny; Nz=a.Nz;
  Xgrid=a.Xgrid; Ygrid=a.Ygrid;
  BiUpToDate = a.BiUpToDate;
  return *this;
}



#endif
