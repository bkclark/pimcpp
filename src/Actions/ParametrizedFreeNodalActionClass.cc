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
    lambdanTau /= ParamList(model,0); // !!!!!!!!!!!!!!!!!!!!!!HACHACKACHAKFHLAKJFLKASJ
  } else if (ModelType == 1) {
    double lambda = Path.Species(SpeciesNum).lambda;
    double nTau = lambdanTau/lambda;
    double lambdaStar = lambda * ParamList(model,0) * (1. + pow(nTau,ParamList(model,1)))
                        /(1. + pow(nTau,ParamList(model,2)));
    lambdanTau = lambdaStar*nTau;
  }

  return (lambdanTau!=0.0) ? (1.0/(4.0*lambdanTau)) : 0.0;
}

double ParametrizedFreeNodalActionClass::ActionImageSum (double L, double fourLambdaBetaInv, double disp, int periodic)
{
  // If the images won't contributed anything, let's not worry about image sums.
  if ((disp*disp*fourLambdaBetaInv) > 50.0 || !periodic) {
    return (disp*disp*fourLambdaBetaInv);
  }

  // Sum over images
  double sum = 0.0;
  int numImages = 10;
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
  if (ModelType == 0 || ModelType == 1 || ModelType == 2) {
    const int nPoints = 10000;
    // Setup grids
    for (int i=0; i<NDIM; i++)
      ActionGrids[i].Init (-0.5*Path.GetBox()[i], 0.5*Path.GetBox()[i], nPoints);

    Array<double,1> actionData(nPoints);
    int nSplines = Path.TotalNumSlices/2 + (Path.TotalNumSlices%2)+1;
    double lambdaTau = Path.tau * Path.Species(SpeciesNum).lambda;

    // Now, setup up actions
    ActionSplines.resize(NumModels,nSplines);
    dVec periodic = Path.GetPeriodic();
    for (int spline=0; spline<nSplines; spline++) {
      double lambdaBeta = lambdaTau * (double)spline;
      double fourLambdaBetaInv = GetFourLambdanTauInv(lambdaBeta);
      for (int dim=0; dim<NDIM; dim++) {
        double L = Path.GetBox()[dim];
        double ActionImageSum0 = ActionImageSum(L, fourLambdaBetaInv, 0.0, periodic[dim]); // Normalizing by madelung constant
        ActionImageSum0 = 0.; // No real need to normalize
        for (int i=0; i<nPoints; i++) {
          double disp = ActionGrids[dim](i);
          actionData(i) = ActionImageSum (L, fourLambdaBetaInv, disp, periodic[dim]) - ActionImageSum0;
          //if (spline == 120 && dim == 0) {
            //cout << actionData(i) << " ";
            //double sum = 0.0;
            //int numImages = 10;
            //for (int image=-numImages; image<=numImages; image++) {
            //  double x = disp + (double)image*L;
            //  sum += exp (-(x*x)*fourLambdaBetaInv);
            //  cout << image << " " << x << " " << -log(sum) << endl;
            //}
          //}
        }
        //if (spline == 120 && dim == 0)
        //  cout << endl;
        // Since the action is periodic, the slope should be zero
        // at the boundaries
        ActionSplines(model,spline)[dim].Init (&ActionGrids[dim], actionData, 0.0, 0.0);
      }
      //if (spline == 120) {
      //  for (int i=0; i<nPoints; i++)
      //    cout << ActionSplines(model,spline)[0](ActionGrids[0](i)) << " ";
      //  cout << endl;
      //}
    }
  }

}

double ParametrizedFreeNodalActionClass::GetRhoij(int slice, int sliceDiff, int refPtcl, int ptcl)
{
  dVec diff;
  double dist;
  Path.RefDistDisp (slice, refPtcl, ptcl, dist, diff);
  double nTau = sliceDiff*Path.tau;
  double lambdanTau = nTau*Path.Species(SpeciesNum).lambda;

  // Free nodes
  double freeAction = 0.;
  if (ModelType == 0 || ModelType == 1 || ModelType == 2) {
    for (int dim=0; dim<NDIM; dim++)
      freeAction += ActionSplines(model,sliceDiff)[dim](diff[dim]);
    //cout << slice << " " << sliceDiff << " " << diff << " " << freeAction << endl;
    if (ModelType == 0 || ModelType == 1)
      return exp(-freeAction);
  }

  // Hydrogen nodes
  if (ModelType == 2) {
    double ionAction = 0.;
    for (int iSpecies=0; iSpecies<Path.IonSpeciesNums.size(); ++iSpecies) {
      double weight = ParamList(model,0);
      double Z = ParamList(model,1);
      double Z6 = pow(Z,6);
      for (int iPtcl=Path.Species(Path.IonSpeciesNums(iSpecies)).FirstPtcl; iPtcl<=Path.Species(Path.IonSpeciesNums(iSpecies)).LastPtcl; iPtcl++) {
        double idist0, idist1;
        Path.RefDistDisp (slice, refPtcl, iPtcl, idist0, diff);
        Path.DistDisp (slice, ptcl, iPtcl, idist1, diff);
        ionAction += weight*exp(-Z*idist0)*exp(-Z*idist1)/Z6; // 1s orbital
      }
    }
    return pow(4.*lambdanTau,-NDIM/2.)*exp(-freeAction) + ionAction;
  }

  // Harmonic ions
  if (ModelType == 3) {
    double w = ParamList(model,0);
    double harmonicAction = 0.;
    for (int iSpecies=0; iSpecies<Path.IonSpeciesNums.size(); ++iSpecies) {
      for (int iPtcl=Path.Species(Path.IonSpeciesNums(iSpecies)).FirstPtcl; iPtcl<=Path.Species(Path.IonSpeciesNums(iSpecies)).LastPtcl; iPtcl++) {
        double idist0, idist1;
        Path.RefDistDisp (slice, refPtcl, iPtcl, idist0, diff);
        Path.DistDisp (slice, ptcl, iPtcl, idist1, diff);
        harmonicAction += ((idist0*idist0 + idist1*idist1)/tanh(nTau*w)) - (2.*idist0*idist1/sinh(nTau*w));
      }
    }
    harmonicAction *= (w/(4.*lambdanTau));
    return exp(-harmonicAction);
  }

  // Harmonic trap
  if (ModelType == 4) {
    double w = ParamList(model,0);
    double dist1 = sqrt(dot(Path(slice,ptcl),Path(slice,ptcl)));
    double dist0 = sqrt(dot(Path.RefPath(refPtcl),Path.RefPath(refPtcl)));
    double harmonicAction = (w/(4.*lambdanTau))*(((dist0*dist0 + dist1*dist1)/tanh(nTau*w)) - (2.*dist0*dist1/sinh(nTau*w)));
    return exp(-harmonicAction);
  }
}

double ParametrizedFreeNodalActionClass::GetRhoij(int slice, int sliceDiff, int refPtcl, int ptcl, Array<dVec,1> &tempPath)
{
  dVec diff;
  double dist;
  Path.RefDistDisp (slice, refPtcl, ptcl, dist, diff, tempPath);
  double nTau = sliceDiff*Path.tau;
  double lambdanTau = nTau*Path.Species(SpeciesNum).lambda;

  // Free nodes
  double freeAction = 0.;
  if (ModelType == 0 || ModelType == 1 || ModelType == 2) {
    for (int dim=0; dim<NDIM; dim++)
      freeAction += ActionSplines(model,sliceDiff)[dim](diff[dim]);
    if (ModelType == 0 || ModelType == 1)
      return exp(-freeAction);
  }

  // Hydrogen nodes
  if (ModelType == 2) {
    double ionAction = 0.;
    for (int iSpecies=0; iSpecies<Path.IonSpeciesNums.size(); ++iSpecies) {
      double weight = ParamList(model,0);
      double Z = ParamList(model,1);
      double Z6 = pow(Z,6);
      for (int iPtcl=Path.Species(Path.IonSpeciesNums(iSpecies)).FirstPtcl; iPtcl<=Path.Species(Path.IonSpeciesNums(iSpecies)).LastPtcl; iPtcl++) {
        double idist0, idist1;
        Path.RefDistDisp (slice, refPtcl, iPtcl, idist0, diff, tempPath);
        Path.DistDisp (slice, ptcl, iPtcl, idist1, diff);
        ionAction += weight*exp(-Z*idist0)*exp(-Z*idist1)/Z6; // 1s orbital
      }
    }
    return pow(4.*lambdanTau,-NDIM/2.)*exp(-freeAction) + ionAction;
  }

  // Harmonic nodes
  if (ModelType == 3) {
    double w = ParamList(model,0);
    double harmonicAction = 0.;
    for (int iSpecies=0; iSpecies<Path.IonSpeciesNums.size(); ++iSpecies) {
      for (int iPtcl=Path.Species(Path.IonSpeciesNums(iSpecies)).FirstPtcl; iPtcl<=Path.Species(Path.IonSpeciesNums(iSpecies)).LastPtcl; iPtcl++) {
        double idist0, idist1;
        Path.RefDistDisp (slice, refPtcl, iPtcl, idist0, diff);
        Path.DistDisp (slice, ptcl, iPtcl, idist1, diff);
        harmonicAction += ((idist0*idist0 + idist1*idist1)/tanh(nTau*w)) - (2.*idist0*idist1/sinh(nTau*w));
      }
    }
    harmonicAction *= (w/(4.*lambdanTau));
    return sqrt(w/(4.*M_PI*Path.Species(SpeciesNum).lambda*sinh(nTau*w)))*exp(-harmonicAction);
  }
}


void ParametrizedFreeNodalActionClass::GetActionDeriv(int slice, int sliceDiff, int refPtcl, int ptcl, dVec &gradPhi, Array<double,2> &detMatrix)
{
  double nTau = sliceDiff*Path.tau;
  double lambdanTau = nTau*Path.Species(SpeciesNum).lambda;

  dVec diff;
  double dist;
  Path.RefDistDisp (slice, refPtcl, ptcl, dist, diff);
  if (ModelType == 0 || ModelType == 1 || ModelType == 2) {
    for (int dim=0; dim<NDIM; dim++)
      gradPhi[dim] = -ActionSplines(model,sliceDiff)[dim].Deriv(diff[dim]) * detMatrix(refPtcl, ptcl);
  }

  if (ModelType == 2) {
    double dIonAction = 0.;
    for (int iSpecies=0; iSpecies<Path.IonSpeciesNums.size(); ++iSpecies) {
      double weight = ParamList(model,0);
      double Z = ParamList(model,1);
      double Z6 = pow(Z,6);
      for (int iPtcl=Path.Species(Path.IonSpeciesNums(iSpecies)).FirstPtcl; iPtcl<=Path.Species(Path.IonSpeciesNums(iSpecies)).LastPtcl; iPtcl++) {
        double idist0, idist1;
        Path.RefDistDisp (slice, refPtcl, iPtcl, idist0, diff);
        Path.DistDisp (slice, ptcl, iPtcl, idist1, diff);
        dIonAction += -Z*weight*exp(-Z*idist0)*exp(-Z*idist1)/Z6; // 1s orbital
      }
    }
    for (int dim=0; dim<NDIM; dim++)
      gradPhi[dim] += dIonAction * detMatrix(refPtcl, ptcl);
  }

  if (ModelType == 3) {
    double w = ParamList(model,0);
    double harmonicAction = 0.;
    double dHarmonicAction = 0.;
    for (int iSpecies=0; iSpecies<Path.IonSpeciesNums.size(); ++iSpecies) {
      for (int iPtcl=Path.Species(Path.IonSpeciesNums(iSpecies)).FirstPtcl; iPtcl<=Path.Species(Path.IonSpeciesNums(iSpecies)).LastPtcl; iPtcl++) {
        double idist0, idist1;
        Path.RefDistDisp (slice, refPtcl, iPtcl, idist0, diff);
        Path.DistDisp (slice, ptcl, iPtcl, idist1, diff);
        harmonicAction += ((idist0*idist0 + idist1*idist1)/tanh(nTau*w)) - (2.*idist0*idist1/sinh(nTau*w));
        dHarmonicAction += ((-2.*idist0/tanh(nTau*w)) + (2*idist1/sinh(nTau*w)));
      }
    }
    harmonicAction *= (w/(4.*lambdanTau));
    dHarmonicAction *= (w/(4.*lambdanTau))*exp(-harmonicAction);
    for (int dim=0; dim<NDIM; dim++)
      gradPhi[dim] = dHarmonicAction * detMatrix(refPtcl, ptcl);
  }

}


void ParametrizedFreeNodalActionClass::GetActionDeriv(int slice, int sliceDiff, int refPtcl, int ptcl, dVec &gradPhi, Array<double,2> &detMatrix, Array<dVec,1> &tempPath)
{
  double nTau = sliceDiff*Path.tau;
  double lambdanTau = nTau*Path.Species(SpeciesNum).lambda;

  dVec diff;
  double dist;
  Path.RefDistDisp (slice, refPtcl, ptcl, dist, diff, tempPath);
  if (ModelType == 0 || ModelType == 1 || ModelType == 2) {
    for (int dim=0; dim<NDIM; dim++)
      gradPhi[dim] = -ActionSplines(model,sliceDiff)[dim].Deriv(diff[dim]) * detMatrix(refPtcl, ptcl);
  }

  if (ModelType == 2) {
    double dIonAction = 0.;
    for (int iSpecies=0; iSpecies<Path.IonSpeciesNums.size(); ++iSpecies) {
      double weight = ParamList(model,0);
      double Z = ParamList(model,1);
      double Z6 = pow(Z,6);
      for (int iPtcl=Path.Species(iSpecies).FirstPtcl; iPtcl<=Path.Species(iSpecies).LastPtcl; iPtcl++) {
        double idist0, idist1;
        Path.RefDistDisp (slice, refPtcl, iPtcl, idist0, diff, tempPath);
        Path.DistDisp (slice, ptcl, iPtcl, idist1, diff);
        dIonAction += -Z*weight*exp(-Z*idist0)*exp(-Z*idist1)/Z6; // 1s orbital
      }
    }
    for (int dim=0; dim<NDIM; dim++)
      gradPhi[dim] += dIonAction * detMatrix(refPtcl, ptcl);
  }

  if (ModelType == 3) {
    double w = ParamList(model,0);
    double harmonicAction = 0.;
    double dHarmonicAction = 0.;
    for (int iSpecies=0; iSpecies<Path.IonSpeciesNums.size(); ++iSpecies) {
      for (int iPtcl=Path.Species(iSpecies).FirstPtcl; iPtcl<=Path.Species(iSpecies).LastPtcl; iPtcl++) {
        double idist0, idist1;
        Path.RefDistDisp (slice, refPtcl, iPtcl, idist0, diff, tempPath);
        Path.DistDisp (slice, ptcl, iPtcl, idist1, diff);
        harmonicAction += ((idist0*idist0 + idist1*idist1)/tanh(nTau*w)) - (2.*idist0*idist1/sinh(nTau*w));
        dHarmonicAction += ((-2.*idist0/tanh(nTau*w)) + (2*idist1/sinh(nTau*w)));
      }
    }
    harmonicAction *= (w/(4.*lambdanTau));
    dHarmonicAction *= (w/(4.*lambdanTau))*exp(-harmonicAction);
    for (int dim=0; dim<NDIM; dim++)
      gradPhi[dim] = dHarmonicAction * detMatrix(refPtcl, ptcl);
  }

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
  out.WriteVar ("ModelType", ModelType);
  out.WriteVar ("NumModels", NumModels);
  out.WriteVar ("NumParams", NumParams);
  out.WriteVar ("ParamList", ParamList);
}


string ParametrizedFreeNodalActionClass::GetName()
{
  return "ParametrizedFreeNodal";
}
