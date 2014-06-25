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

#ifndef PARAMETRIZED_FREE_NODAL_ACTION_CLASS_H
#define PARAMETRIZED_FREE_NODAL_ACTION_CLASS_H

//#include <omp.h>
#include "NodalActionClass.h"
#include "../Splines/CubicSpline.h"

/// ParametrizedFreeNodalActionClass implements the nodal action corresponding to
/// the free fermion density matrix in periodic boundary conditions.
/// Currently, the beta-derivative of the action is computed only
/// approximately, but in a way that becomes accurate as \f$ \tau
/// \rightarrow 0 \f$.
class ParametrizedFreeNodalActionClass : public NodalActionClass
{
private:
  PathClass &Path;

  /// These splines will hold the free-particle action for
  /// periodic boundary conditions.  The array is over time-slice
  /// separation from the reference slice.
  Array<TinyVector<CubicSpline,NDIM>,2> ActionSplines;
  TinyVector<LinearGrid,NDIM> ActionGrids;
  double ActionImageSum (double L, double lambdaTau, double disp, int periodic);
  double ActionkSum (double L, double lambdaTau, double disp);

  double GetRhoij(int slice, int sliceDiff, int refPtcl, int ptcl);
  double GetRhoij(int slice, int sliceDiff, int refPtcl, int ptcl, Array<dVec,1> &tempPath);
  void GetActionDeriv(int slice, int sliceDiff, int refPtcl, int ptcl, dVec &gradPhi, Array<double,2> &detMatrix);
  void GetActionDeriv(int slice, int sliceDiff, int refPtcl, int ptcl, dVec &gradPhi, Array<double,2> &detMatrix, Array<dVec,1> &tempPath);

  double GetFourLambdanTauInv (double lambdanTau);
  void SetupActions();
public:
  /// Variational parameters
  int NumModels, NumParams, ModelType;
  Array<double,2> ParamList;
  int model;
  void ChangeModel(int tmpModel);
  int GetModel();
  int GetNumModels();
  void Read (IOSectionClass &in);

  //  Array<double,2> GetMatrix(int slice);
  NodeType Type();
  void WriteInfo (IOSectionClass &out);
  string GetName();
  ParametrizedFreeNodalActionClass (PathDataClass &pathData, int speciesNum);
};




#endif
