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

#ifndef VARIATIONAL_PI__CLASS_H
#define VARIATIONAL_PI__CLASS_H

#include "NodalActionClass.h"
#include <Common/Splines/CubicSpline.h>


/// VartionalPIClass implements the nodal action corresponding to
/// the free fermion density matrix in periodic boundary conditions.
/// Currently, the beta-derivative of the action is computed only
/// approximately, but in a way that becomes accurate as \f$ \tau
/// \rightarrow 0 \f$.
class VariationalPIClass : public NodalActionClass
{
private:
  PathClass &Path;
  void calc_u();
  Array<double,1> u;
  Array<double,1> newCol;
  /// These splines will hold the free-particle action for
  /// periodic boundary conditions.  The array is over time-slice
  /// separation from the reference slice.
  Array<double,2> DetMatrix;
public:
  double SingleAction (int slice1, int slice2,
		       const Array<int,1> &activeParticles,
		       int level);
  double d_dBeta (int slice1, int slice2, 
		  int level);
  string GetName();
  void BuildDeterminantMatrix();
  void CheckDeterminantMatrix();
  void calc_u(const Array<int,1> &changePtcls);
  void Read (IOSectionClass &in);
  bool IsGroundState();
  bool IsPositive(int x);
  NodeType Type();
  void AcceptCopy(int slice1, int slice2);
  void RejectCopy(int slice1, int slice2);
  void WriteInfo (IOSectionClass &out);
  int ChangedColumn;
  
  VariationalPIClass (PathDataClass &pathData);
};

#endif

