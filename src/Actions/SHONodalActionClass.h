/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2013  B. Clark, K. Esler, E. Brown   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the SHO Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
// http://code.google.com/p/pimcplusplus/                  //
/////////////////////////////////////////////////////////////

#ifndef SHO_NODAL_ACTION_CLASS_H
#define SHO_NODAL_ACTION_CLASS_H

#include "NodalActionClass.h"
#include "../Splines/CubicSpline.h"

/// SHONodalActionClass implements the nodal action corresponding to
/// the free fermion density matrix in periodic boundary conditions.
/// Currently, the beta-derivative of the action is computed only
/// approximately, but in a way that becomes accurate as \f$ \tau
/// \rightarrow 0 \f$.
class SHONodalActionClass : public NodalActionClass
{
private:
  PathClass &Path;

  void SetupActions();

  Array<double,1> c1, c2, c3; // constants
  double omega;

  double GetRhoij(int slice, int sliceDiff, int refPtcl, int ptcl);
  double GetRhoij(int slice, int sliceDiff, int refPtcl, int ptcl, Array<dVec,1> &tempPath);
  void GetActionDeriv(int slice, int sliceDiff, int refPtcl, int ptcl, dVec &gradPhi, Array<double,2> &detMatrix);
  void GetActionDeriv(int slice, int sliceDiff, int refPtcl, int ptcl, dVec &gradPhi, Array<double,2> &detMatrix, Array<dVec,1> &tempPath);
public:
  void Read (IOSectionClass &in);
  NodeType Type();
  void WriteInfo (IOSectionClass &out);
  string GetName();
  SHONodalActionClass (PathDataClass &pathData, int speciesNum);
};


#endif
