/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2013  B. Clark, K. Esler, E. Brown   //
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
// http://code.google.com/p/pimcplusplus/                  //
/////////////////////////////////////////////////////////////

#include "Integrate.h"


/* Note:  The first column of result must have the starting */
/* values already stored when called */
void RungeKutta4 (const Grid &grid, 
		  int StartPoint, int EndPoint,
		  Array<pscalar, 2> &Result,
		  Array<pscalar,1> (*DerivFunc)(pscalar x,
					       Array<pscalar,1> y,
					       void *Params),
		  void *Params)
{
  int NumVars = Result.rows();
  Array <pscalar,1> k1(NumVars), k2(NumVars), k3(NumVars), k4(NumVars);
  Array <pscalar,1> y(NumVars), yplus(NumVars);
  pscalar h, x;
  int direction;

  if (StartPoint < EndPoint)
    direction = 1;
  else
    direction = -1;
  
  const static pscalar OneSixth = 1.0/6.0;
  const static pscalar OneThird = 1.0/3.0;

  int i = StartPoint;
  while (i!=EndPoint)
    {
      x = grid(i);
      h = grid(i+direction) - x;
      y = Result (Range::all(), i);
      k1 = h * (*DerivFunc)(x,       y, Params);
      yplus = y + 0.5*k1;
      k2 = h * (*DerivFunc)(x+0.5*h, yplus, Params);
      yplus = y + 0.5*k2;
      k3 = h * (*DerivFunc)(x+0.5*h, yplus, Params);
      yplus = y + k3;
      k4 = h * (*DerivFunc)(x+h, yplus, Params);
      Result(Range::all(), i+direction) = 
	y + (OneSixth*(k1+k4) + OneThird*(k2+k3));
      int TooBig = 0;
      for (int j=0; j<Result.rows(); j++)
	if (fabs(Result(j, i+direction)) > 1.0e50)
	  TooBig=1;
      if (TooBig)
	{
	  for (int j=StartPoint; j!=(EndPoint+direction); j+=direction)
	    Result(Range::all(), j) *= 1.0e-50;
	  cerr << "Integration overflow.\n";
	}
	
      i += direction;
    }
}



void IntegrateFirstOrderNS (const Grid &grid,
			   int StartPoint, int EndPoint,
			   Array<pscalar,1> &Result,
			   pscalar (*DerivFunc)(pscalar x, pscalar y,
					     void *Params),
			   void *Params)
{
  pscalar k1, k2, k3, k4;
  pscalar y, yplus;
  pscalar h, x;
  int direction;

  if (StartPoint < EndPoint)
    direction = 1;
  else
    direction = -1;
  
  const static pscalar OneSixth = 1.0/6.0;
  const static pscalar OneThird = 1.0/3.0;

  int i = StartPoint;
  while (i!=EndPoint)
    {
      x = grid(i);
      h = grid(i+direction) - x;
      y = Result (i);
      k1 = h * (*DerivFunc)(x,       y, Params);
      yplus = y + 0.5*k1;
      k2 = h * (*DerivFunc)(x+0.5*h, yplus, Params);
      yplus = y + 0.5*k2;
      k3 = h * (*DerivFunc)(x+0.5*h, yplus, Params);
      yplus = y + k3;
      k4 = h * (*DerivFunc)(x+h, yplus, Params);
      Result(i+direction) = 
	y + (OneSixth*(k1+k4) + OneThird*(k2+k3));
      i += direction;
    }
}





void IntegrateFirstOrder (const Grid &grid,
			   int StartPoint, int EndPoint,
			   Array<pscalar,1> &Result,
			   pscalar (*DerivFunc)(pscalar x, pscalar y,
					     void *Params),
			   void *Params)
{
  pscalar k1, k2, k3, k4;
  pscalar y, yplus;
  pscalar h, x;
  int direction;

  if (StartPoint < EndPoint)
    direction = 1;
  else
    direction = -1;
  
  const static pscalar OneSixth = 1.0/6.0;
  const static pscalar OneThird = 1.0/3.0;

  int i = StartPoint;
  while (i!=EndPoint)
    {
      x = grid(i);
      h = grid(i+direction) - x;
      y = Result (i);
      k1 = h * (*DerivFunc)(x,       y, Params);
      yplus = y + 0.5*k1;
      k2 = h * (*DerivFunc)(x+0.5*h, yplus, Params);
      yplus = y + 0.5*k2;
      k3 = h * (*DerivFunc)(x+0.5*h, yplus, Params);
      yplus = y + k3;
      k4 = h * (*DerivFunc)(x+h, yplus, Params);
      Result(i+direction) = 
	y + (OneSixth*(k1+k4) + OneThird*(k2+k3));
      int TooBig = 0;
      if (fabs(Result(i+direction)) > 1.0e10)
	for (int j=StartPoint; j!=(EndPoint+direction); j+=direction)
	  Result(j) *= 1.0e-10;
      i += direction;
    }
}




void IntegrateSecondOrder (const Grid &grid,
			   int StartPoint, int EndPoint,
			   Array<Vec2,1> &Result,
			   Vec2 (*DerivFunc)(pscalar x, Vec2 y,
					     void *Params),
			   void *Params)
{
  Vec2 k1, k2, k3, k4;
  Vec2 y, yplus;
  pscalar h, x;
  int direction;

  if (StartPoint < EndPoint)
    direction = 1;
  else
    direction = -1;
  
  const static pscalar OneSixth = 1.0/6.0;
  const static pscalar OneThird = 1.0/3.0;

  int i = StartPoint;
  while (i!=EndPoint)
    {
      x = grid(i);
      h = grid(i+direction) - x;
      y = Result (i);
      k1 = h * (*DerivFunc)(x,       y, Params);
      yplus = y + 0.5*k1;
      k2 = h * (*DerivFunc)(x+0.5*h, yplus, Params);
      yplus = y + 0.5*k2;
      k3 = h * (*DerivFunc)(x+0.5*h, yplus, Params);
      yplus = y + k3;
      k4 = h * (*DerivFunc)(x+h, yplus, Params);
      Result(i+direction) = 
	y + (OneSixth*(k1+k4) + OneThird*(k2+k3));
      int TooBig = 0;
      if ( (fabs(Result(i+direction)[0]) > 1.0e10) ||
	   (fabs(Result(i+direction)[1]) > 1.0e10))
	for (int j=StartPoint; j!=(EndPoint+direction); j+=direction)
	  Result(j) *= 1.0e-10;
      i += direction;
    }
}
