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


#include "SHONodalActionClass.h"
#include "../PathDataClass.h"
#include "../MatrixOps/MatrixOps.h"


SHONodalActionClass::SHONodalActionClass (PathDataClass &pathData,int speciesNum) :
  NodalActionClass (pathData), Path (pathData.Path)
{
  SpeciesNum = speciesNum;
  int N = Path.Species(speciesNum).LastPtcl - Path.Species(speciesNum).FirstPtcl+1;
  DetMatrix.resize(N,N);
  Cofactors.resize(N,N);
  GradVec.resize(N);
  SavePath.resize(N);
}


void SHONodalActionClass::Read (IOSectionClass &in)
{
  if (!in.ReadVar ("omega",omega))
    omega = 1.0;

  NodalActionClass::Read(in);
}


void SHONodalActionClass::SetupActions()
{
  double tau = Path.tau;
  c1.resize(Path.TotalNumSlices);
  c2.resize(Path.TotalNumSlices);
  c3.resize(Path.TotalNumSlices);
  for (int i = 1; i < Path.TotalNumSlices; ++i) {
    //c1(i) = pow((omega/(2.0*pi*sinh(i*tau*omega))),NDIM/2.0);
    //c2(i) = omega/(2.0*sinh(i*tau*omega));
    //c3(i) = cosh(i*tau*omega);
    c1(i) = omega/(4.*Path.Species(SpeciesNum).lambda*i*tau);
    c2(i) = 1./tanh(i*tau*omega);
    c3(i) = 2./sinh(i*tau*omega);
  }
}


double SHONodalActionClass::GetRhoij(int slice, int sliceDiff, int refPtcl, int ptcl)
{
  double dist1 = sqrt(dot(Path(slice,ptcl),Path(slice,ptcl)));
  double dist0 = sqrt(dot(Path.RefPath(refPtcl),Path.RefPath(refPtcl)));
  double rhoij = exp(-c1(sliceDiff)*(c2(sliceDiff)*(dist0*dist0+dist1*dist1) - c3(sliceDiff)*dist0*dist1));

  return rhoij;
}


double SHONodalActionClass::GetRhoij(int slice, int sliceDiff, int refPtcl, int ptcl, Array<dVec,1> &tempPath)
{
  double dist1 = sqrt(dot(tempPath(slice,ptcl),tempPath(slice,ptcl)));
  double dist0 = sqrt(dot(Path.RefPath(refPtcl),Path.RefPath(refPtcl)));
  double rhoij = exp(-c1(sliceDiff)*(c2(sliceDiff)*(dist0*dist0+dist1*dist1) - c3(sliceDiff)*dist0*dist1));

  return rhoij;
}


void SHONodalActionClass::GetActionDeriv(int slice, int sliceDiff, int refPtcl, int ptcl, dVec &gradPhi, Array<double,2> &detMatrix)
{
  for (int dim=0; dim<NDIM; dim++)
    gradPhi[dim] = -GetRhoij(slice,sliceDiff,refPtcl,ptcl) * detMatrix(refPtcl, ptcl);
}


void SHONodalActionClass::GetActionDeriv(int slice, int sliceDiff, int refPtcl, int ptcl, dVec &gradPhi, Array<double,2> &detMatrix, Array<dVec,1> &tempPath)
{
  for (int dim=0; dim<NDIM; dim++)
    gradPhi[dim] = -GetRhoij(slice,sliceDiff,refPtcl,ptcl,tempPath) * detMatrix(refPtcl, ptcl);
}


NodeType SHONodalActionClass::Type()
{
  return FREE_PARTICLE;
}


void SHONodalActionClass::WriteInfo (IOSectionClass &out)
{
  out.WriteVar ("Type", "SHO");
}


string SHONodalActionClass::GetName()
{
  return "SHONodal";
}
