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

// #include "../MPI/Communication.h"
#include "ParametrizedFreeNodalActionClass.h"
#include "../PathDataClass.h"
#include "../MatrixOps/MatrixOps.h"


ParametrizedFreeNodalActionClass::ParametrizedFreeNodalActionClass (PathDataClass &pathData, int speciesNum) :
  NodalActionClass (pathData), Path (pathData.Path)
{
  SpeciesNum = speciesNum;
  int N = Path.Species(speciesNum).LastPtcl - Path.Species(speciesNum).FirstPtcl+1;
  DetMatrix.resize(N,N);
  Cofactors.resize(N,N);
  GradVec.resize(N);
  SavePath.resize(N);
}


void ParametrizedFreeNodalActionClass::Read (IOSectionClass &in)
{
  assert (in.ReadVar ("ModelType", ModelType));
  assert (in.ReadVar ("NumModels", NumModels));
  assert (in.ReadVar ("NumParams", NumParams));
  assert (in.ReadVar ("ParamList", ParamList));

  // First model will be setup by base class
  model = 0;
  NodalActionClass::Read(in);

  // Do first model again, just in case (probably could do this better)
  for (int i=0; i<NumModels; i++) {
    model = i;
    SetupActions();
  }

}

double ParametrizedFreeNodalActionClass::GetFourLambdanTauInv (double lambdanTau)
{
  if (ModelType == 0) {
    lambdanTau = lambdanTau*ParamList(model,0);
  } else if (ModelType == 1) {
    double lambda = Path.Species(SpeciesNum).lambda;
    double nTau = lambdanTau/lambda;
    double lambdaStar = lambda * (1. + ParamList(model,0)*pow(nTau,ParamList(model,2)))
                        /(1. + ParamList(model,1)*pow(nTau,ParamList(model,2)));
    lambdanTau = lambdaStar*nTau;
  }

  return (lambdanTau!=0.0) ? (1.0/(4.0*lambdanTau)) : 0.0;
}


double ParametrizedFreeNodalActionClass::ActionImageSum (double L, double lambdaBeta,  double disp)
{
  int numImages = 10;
  double sum = 0.0;
  //double fourLambdaBetaInv = (lambdaBeta!=0.0) ? ((ParamList(model,0) * 1.0/(4.0*lambdaBeta))) : 0.0;
  double fourLambdaBetaInv = GetFourLambdanTauInv(lambdaBeta);
  // If the images won't contributed anything, let's not worry
  // about image sums.
  if ((disp*disp*fourLambdaBetaInv) > 50.0)
    return (disp*disp*fourLambdaBetaInv);
  for (int image=-numImages; image<=numImages; image++) {
    double x = disp + (double)image*L;
    sum += exp (-(x*x)*fourLambdaBetaInv);
  }
  return (-log(sum));
}


double ParametrizedFreeNodalActionClass::ActionkSum (double L, double lambdaBeta, double disp)
{
  double kmax = sqrt (50.0/lambdaBeta);
  double kInc = 2.0*M_PI/L;

  if (lambdaBeta == 0.0)
    return (0.0);

  double sum = 0.5;
  for (double k=kInc; k<=kmax; k+=kInc)
    sum += cos(k*disp)*exp(-lambdaBeta*k*k);
  return (-log(sum));
}


void ParametrizedFreeNodalActionClass::SetupActions()
{
  const int nPoints = 1000;
  // Setup grids
  for (int i=0; i<NDIM; i++)
    ActionGrids[i].Init (-0.5*Path.GetBox()[i], 0.5*Path.GetBox()[i], nPoints);

  Array<double,1> actionData(nPoints);
  int nSplines = Path.TotalNumSlices/2 + (Path.TotalNumSlices%2)+1;
  double lambdaTau = Path.tau * Path.Species(SpeciesNum).lambda;

  // Now, setup up actions
  ActionSplines.resize(NumModels,nSplines);
  for (int spline=0; spline<nSplines; spline++) {
    double lambdaBeta = lambdaTau * (double)spline;
    for (int dim=0; dim<NDIM; dim++) {
      double L = Path.GetBox()[dim];
      for (int i=0; i<nPoints; i++) {
        double disp = ActionGrids[dim](i);
        actionData(i) = ActionImageSum (L, lambdaBeta, disp) - ActionImageSum(L, lambdaBeta, 0.0);
      }
      // Since the action is periodic, the slope should be zero
      // at the boundaries
      ActionSplines(model,spline)[dim].Init (&ActionGrids[dim], actionData, 0.0, 0.0);
    }
  }
}


double ParametrizedFreeNodalActionClass::GetAction(int slice, int sliceDiff, int refPtcl, int ptcl)
{
  dVec diff;
  double dist;
  Path.RefDistDisp (slice, refPtcl, ptcl, dist, diff);
  double action = 0.0;
  for (int dim=0; dim<NDIM; dim++)
    action += ActionSplines(model,sliceDiff)[dim](diff[dim]);
  return action;
}


double ParametrizedFreeNodalActionClass::GetAction(int slice, int sliceDiff, int refPtcl, int ptcl, Array<dVec,1> &tempPath)
{
  dVec diff;
  double dist;
  Path.RefDistDisp (slice, refPtcl, ptcl, dist, diff, tempPath);
  double action = 0.0;
  for (int dim=0; dim<NDIM; dim++)
    action += ActionSplines(model,sliceDiff)[dim](diff[dim]);
  return action;
}


void ParametrizedFreeNodalActionClass::GetActionDeriv(int slice, int sliceDiff, int refPtcl, int ptcl, dVec &gradPhi, Array<double,2> &detMatrix)
{
  dVec diff;
  double dist;
  Path.RefDistDisp (slice, refPtcl, ptcl, dist, diff);
  for (int dim=0; dim<NDIM; dim++)
    gradPhi[dim] = -ActionSplines(model,sliceDiff)[dim].Deriv(diff[dim]) * detMatrix(refPtcl, ptcl);
}


void ParametrizedFreeNodalActionClass::GetActionDeriv(int slice, int sliceDiff, int refPtcl, int ptcl, dVec &gradPhi, Array<double,2> &detMatrix, Array<dVec,1> &tempPath)
{
  dVec diff;
  double dist;
  Path.RefDistDisp (slice, refPtcl, ptcl, dist, diff, tempPath);
  for (int dim=0; dim<NDIM; dim++)
    gradPhi[dim] = -ActionSplines(model,sliceDiff)[dim].Deriv(diff[dim]) * detMatrix(refPtcl, ptcl);
}


// Change the model
void ParametrizedFreeNodalActionClass::ChangeModel(int tmpModel)
{
  model = tmpModel;
}


// Return current model
int ParametrizedFreeNodalActionClass::GetModel()
{
  return model;
}


// Return number of models
int ParametrizedFreeNodalActionClass::GetNumModels()
{
  return NumModels;
}


NodeType ParametrizedFreeNodalActionClass::Type()
{
  return FREE_PARTICLE;
}


void ParametrizedFreeNodalActionClass::WriteInfo (IOSectionClass &out)
{
  //out.WriteVar ("Type", "PARAMETRIZED_FREE_PARTICLE");
}


string ParametrizedFreeNodalActionClass::GetName()
{
  return "ParametrizedFreeNodal";
}
