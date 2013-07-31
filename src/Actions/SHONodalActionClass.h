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

  void SetupSHOActions();

  Array<double,2> DetMatrix, Cofactors;
  Array<dVec,1> GradVec, SavePath;
  void GradientDet (int slice, double &det, Array<dVec,1> &gradient);
  void GradientDet (int slice, double &det, Array<dVec,1> &gradient, Array<dVec,1> &tempPath);
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
  int nSingular;
  bool FirstDistTime, FirstDetTime;

  Array<double,1> c1, c2, c3; // constants
  double omega;
  double GetAction(int refSlice, int slice, int sliceDiff, int refPtcl, int ptcl);

public:
  void Init();
  double GetNodeDist (int slice, double lambda, double levelTau, int SpeciesNum);
  bool UseHybridDist, UseNewtonRaphsonDist, UseLineSearchDist, UseMaxDist, UseNoDist;
  double SingleAction (int slice1, int slice2, const Array<int,1> &activeParticles, int level);
  double SimpleAction (int slice1, int slice2, const Array<int,1> &activeParticles, int level);
  double PreciseAction (int slice1, int slice2, const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  /// Returns true if the nodal restriction is satisfied for my
  /// species at timeslice slice.  If slice is the reference slice,
  /// returns true.
  bool IsPositive (int slice);
  double Det (int slice);
  double Det (int slice, Array<dVec,1> &tempPath);
  //  Array<double,2> GetMatrix(int slice);
  void Read (IOSectionClass &in);
  bool IsGroundState();
  NodeType Type();
  void WriteInfo (IOSectionClass &out);
  string GetName();
  SHONodalActionClass (PathDataClass &pathData, int speciesNum);
};


#endif
