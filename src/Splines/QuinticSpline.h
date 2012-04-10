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

#ifndef QUINTIC_SPLINE_H
#define QUINTIC_SPLINE_H

#include "Grid.h"
#include "../nan.h"
#include <iostream>

#define F77_FUNC(name,NAME) NAME
/// The QuinticSpline class is a third-order spline representation of a
/// function.  It stores a pointer to a grid and the values of the
/// function and its second derivative at the points defined by the
/// grid. 
class QuinticSpline
{
private:
  /// This flag records whether or not the stored second derivatives
  /// are in sync with the function values.  It is used to determine
  /// whether the second derivatives need recomputation.
  int UpToDate;
  /// The function values on the grid points.
  Array<double, 1> Y;
  /// The parameters to the fortran routine QUINAT
  Array<double, 1> FX, FY, FB, FC, FD, FE, FF;
  int offset;
  /// The start and end first and second derivatives of the function.
  /// If these are NAN's, then we assume natural boundary conditions.
  double StartDeriv, StartDeriv2, EndDeriv, EndDeriv2;
  /// The coefficients of the 5th-order polynomial
  Array<double,1>  B, C, D, E, F;
  int I, J;
  inline void GetIJ(double x);
  Grid *grid;
public:
  int NumParams;

  /// Returns the interpolated value.
  inline double operator()(double x);
  /// Returns the interpolated first derivative.
  inline double Deriv(double x);
  /// Returns the interpolated second derivative.
  inline double Deriv2(double x);
  /// Returns the interpolated third derivative.
  inline double Deriv3(double x);
  /// Returns the interpolated fourth derivative.
  inline double Deriv4(double x);
  /// Returns the number of points;
  inline int NumPoints() const;

  inline Array<double,1>& Data();
  inline const Array<double,1>& Data() const;
  /// Recompute the quintic polynomial coefficients
  void Update();
  
  /// Initialize the cubic spline.  See notes about start and end
  /// deriv above.
  inline void Init(Grid *NewGrid, Array<double,1> NewYs,
		   double startderiv, double endderiv,
		   double startderiv2, double endderiv2);

  /// Simplified form which assumes that the second derivative at both
  /// boundaries are zero.
  inline void Init (Grid *NewGrid, Array<double,1> NewYs)
  {
    Init (NewGrid, NewYs, NAN, NAN, NAN, NAN);
  }
  
  /// Simplified constructor.
  inline QuinticSpline (Grid *NewGrid, Array<double,1> NewYs)
  {
    Init (NewGrid, NewYs, NAN, NAN, NAN, NAN);
  }

  /// Full constructor.
  inline QuinticSpline (Grid *NewGrid, Array<double,1> NewYs,
			double startderiv, double startderiv2,
			double endderiv, double endderiv2)
  {
    Init (NewGrid, NewYs, startderiv, startderiv2, endderiv, endderiv2);
    Update();
  }

  /// Returns the value of the function at the ith grid point.
  inline double operator()(int i) const
  {
    return (Y(i));
  }
  /// Returns a reference to the value at the ith grid point.
  inline double & operator()(int i)
  {
    UpToDate = false;
    return (Y(i));
  }
  void Write(IOSectionClass &outSection);
  void Read(IOSectionClass &inSection);
  QuinticSpline& operator=(const QuinticSpline& spline);

  /// Trivial constructor
  QuinticSpline()
  {
    UpToDate = false;
  }
};



inline int QuinticSpline::NumPoints() const 
{
  return Y.size();
}

void QuinticSpline::Init(Grid *NewGrid, Array<double,1> NewY,
			 double startderiv, double endderiv,
			 double startderiv2, double endderiv2)
{
  StartDeriv = startderiv; StartDeriv2 = startderiv2;
  EndDeriv   = endderiv;   EndDeriv2   = endderiv2;
  grid = NewGrid;
  if (NewGrid->NumPoints != NewY.size()) {
    cerr << "Size mismatch in QuinticSpline.\n";
    cerr << "Number of grid points = " << NewGrid->NumPoints << endl;
    cerr << "Number of Y's         = " << NewY.size() << endl;
    exit(1);
  }
  Y.resize(NewY.size());
  B.resize(NewY.size());
  C.resize(NewY.size());
  D.resize(NewY.size());
  E.resize(NewY.size());
  F.resize(NewY.size());
  Y = NewY;

  NumParams = grid->NumPoints;
  if (!myIsNAN(StartDeriv) /*&& !isinf(StartDeriv)*/)
    {
      NumParams++;
      if (!myIsNAN(StartDeriv2) /*&& !isinf(StartDeriv2)*/)
	NumParams++;
    }
  if (!myIsNAN(EndDeriv) /* && !isinf(EndDeriv)*/)
    {
      NumParams++;
      if (!myIsNAN(EndDeriv2) /* && !isinf(EndDeriv2)*/)
	NumParams++;
    }
  FX.resize(NumParams);
  FY.resize(NumParams);
  FB.resize(NumParams);
  FC.resize(NumParams);
  FD.resize(NumParams);
  FE.resize(NumParams);
  FF.resize(NumParams);
  UpToDate = false;
}


inline void QuinticSpline::GetIJ(double x)
{    
  Grid &X = *grid;
#ifdef BZ_DEBUG
  if (x > X.End)
    {
      if (x < (X.End * 1.000000001))
	x = X.End;
      else
	{
	  cerr << "x outside grid in QuinticSpline.\n";
	  cerr << "x = " << x << " X.End = " << X.End << "\n";
	  exit(1);
	}
    }
#endif
  J = X.ReverseMap(x)+1;
  I = J-1;
  if (I<0)
    {
      I = 0;
      J = 1;
    }
  if (J>(X.NumPoints-1))
    {
      J = (X.NumPoints-1);
      I = J-1;
    }
}



inline double QuinticSpline::operator()(double x)
{
  Grid &X = *grid;
  if (!UpToDate)
    Update();

  GetIJ(x);
  
  double P = (x-X(I));

//   cerr << "x = " << x << "I = " << I << endl;
//   cerr << "P = " << P << endl;

  double S = ((((F(I)*P+E(I))*P+D(I))*P+C(I))*P+B(I))*P+Y(I);
  return (S);
}



inline double QuinticSpline::Deriv(double x)
{
  Grid &X = *grid;
  if (!UpToDate)
    Update();

  GetIJ(x);
  double P = x-X(I);

  double S = (((5.0*F(I)*P+4.0*E(I))*P+3.0*D(I))*P+2.0*C(I))*P+B(I);
  return (S);
}



inline double QuinticSpline::Deriv2(double x)
{
  Grid &X = *grid;
  if (!UpToDate)
    Update();

  GetIJ(x);
  double P = (x-X(I));

  double S = ((20.0*F(I)*P+12.0*E(I))*P+6.0*D(I))*P+2.0*C(I);
  return (S);
}




inline double QuinticSpline::Deriv3(double x)
{
  Grid &X = *grid;
  if (!UpToDate)
    Update();

  GetIJ(x);
  double P = (x-X(I));
  double S = (60.0*F(I)*P+24.0*E(I))*P+6.0*D(I);
  return (S);
}


inline double QuinticSpline::Deriv4(double x)
{
  Grid &X = *grid;
  if (!UpToDate)
    Update();

  GetIJ(x);
  double P = (x-X(I));
  double S = 120.0*F(I)*P+24.0*E(I);
  return (S);
}

inline Array<double,1>& QuinticSpline::Data()
{
  UpToDate = false;
  return (Y);
}

inline const Array<double,1>& QuinticSpline::Data() const
{
  return (Y);
}


#endif
