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

#ifndef NODAL_ACTION_CLASS_H
#define NODAL_ACTION_CLASS_H

#include "ActionBase.h"
#include "../Splines/CubicSpline.h"

  class PathClass;

typedef enum { FREE_PARTICLE, GROUND_STATE, GROUND_STATE_FP } NodeType;

class NodalActionClass : public ActionBaseClass
{
public:
  virtual void ChangeModel(int tmpModel);
  virtual int GetModel();
  virtual int GetNumModels();

  virtual bool IsPositive (int slice);
  //  virtual double Det(int slice)       = 0;
  //  virtual Array<double,2> GetMatrix (int slice=0) = 0;
  void AcceptCopy (int slice1, int slice2);
  void RejectCopy (int slice1, int slice2);
  virtual void Init();
  virtual bool IsGroundState();
  virtual NodeType Type() = 0;
  virtual void Setk (Vec3 kVec);
  virtual void Update();
  virtual void Read (IOSectionClass &in);
  virtual void SetupActions();
  virtual double GetNodeDist (int slice, double lambda, double levelTau, int SpeciesNum);
  virtual double NodalDist (int slice);
  virtual double HybridDist(int slice, double lambdaTau);
  /// This returns the upper bound on the distance to a node by
  /// returning the minimum distance to particle coincidence.
  virtual double MaxDist(int slice);
  /// This calculates the distance to the node along the line in the
  /// direction of the gradient by a bisection search
  virtual double LineSearchDist (int slice);
  /// This calculates the nodal distance by an iterated Newton-Raphson approach
  virtual double NewtonRaphsonDist (int slice);
  virtual double SingleAction (int slice1, int slice2, const Array<int,1> &activeParticles, int level);
  virtual double SimpleAction (int slice1, int slice2, const Array<int,1> &activeParticles, int level);
  virtual double PreciseAction (int slice1, int slice2, const Array<int,1> &activeParticles, int level);
  virtual double NodeImportanceAction (int slice1, int slice2, const Array<int,1> &activeParticles, int level);
  virtual double d_dBeta (int slice1, int slice2, int level);

  virtual double Det (int slice);
  virtual double Det (int slice, Array<dVec,1> &tempPath);
  virtual void GradientDet (int slice, double &det, Array<dVec,1> &gradient);
  virtual void GradientDet (int slice, double &det, Array<dVec,1> &gradient, Array<dVec,1> &tempPath);
  virtual void GradientDetFD (int slice, double &det, Array<dVec,1> &gradient);

  virtual double GetRhoij(int slice, int sliceDiff, int refPtcl, int ptcl);
  virtual double GetRhoij(int slice, int sliceDiff, int refPtcl, int ptcl, Array<dVec,1> &tempPath);
  virtual void GetActionDeriv(int slice, int sliceDiff, int refPtcl, int ptcl, dVec &gradPhi, Array<double,2> &detMatrix);
  virtual void GetActionDeriv(int slice, int sliceDiff, int refPtcl, int ptcl, dVec &gradPhi, Array<double,2> &detMatrix, Array<dVec,1> &tempPath);
  virtual void WriteInfo(IOSectionClass &out);

  bool UseHybridDist, UseNewtonRaphsonDist, UseLineSearchDist, UseMaxDist, UseNoDist;
  bool FirstDistTime, FirstDetTime;
  int SpeciesNum;
  int nSingular;
  int NumGradDists, NumLineDists;
  Array<double,2> DetMatrix, Cofactors;
  Array<dVec,1> GradVec, SavePath;

  NodalActionClass (PathDataClass &pathData) :
    ActionBaseClass (pathData)
  {
    NumLineDists = NumGradDists = 0;
    nSingular = 0;
  }
};


#endif
