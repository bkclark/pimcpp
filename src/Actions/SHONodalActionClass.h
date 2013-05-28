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

  /// These splines will hold the free-particle action for
  /// periodic boundary conditions.  The array is over time-slice
  /// separation from the reference slice.
  Array<TinyVector<CubicSpline,NDIM>,1> ActionSplines;
  TinyVector<LinearGrid,NDIM> ActionGrids;
  void SetupSHOActions();
  double ActionImageSum (double L, double lambdaTau, double disp);
  double ActionkSum (double L, double lambdaTau, double disp);

  Array<double,2> DetMatrix, Cofactors;
  Array<dVec,1> GradVec, SavePath;
  void GradientDet (int slice, double &det, Array<dVec,1> &gradient);
  void GradientDetFD (int slice, double &det, Array<dVec,1> &gradient);
  double NodalDist (int slice);
  double HybridDist(int slice, double lambdaTau);
  /// This returns the upper bound on the distance to a node by
  /// returning the minimum distance to particle coincidence.
  double MaxDist(int slice);
  /// This calculates the distance to the node along the line in the
  /// direction of the gradient by a bisection search
  double LineSearchDist (int slice);
  /// This calculates the nodal distance by an iterated Newton-Raphson
  /// approach
  double NewtonRaphsonDist (int slice);
  int SpeciesNum;
  int NumGradDists, NumLineDists;
public:
  double SingleAction (int slice1, int slice2,
		       const Array<int,1> &activeParticles,
		       int level);
  double d_dBeta (int slice1, int slice2, 
		  int level);
  double SimpleAction (int slice1, int slice2,
		       const Array<int,1> &activeParticles,
		       int level);
  /// Returns true if the nodal restriction is satisfied for my
  /// species at timeslice slice.  If slice is the reference slice,
  /// returns true.
  bool IsPositive (int slice);
  double Det (int slice);
  //  Array<double,2> GetMatrix(int slice);
  void Read (IOSectionClass &in);
  bool IsGroundState();
  NodeType Type();
  void WriteInfo (IOSectionClass &out);
  string GetName();
  SHONodalActionClass (PathDataClass &pathData, int speciesNum);
};




#endif
