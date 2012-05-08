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

#include "WindingNumber.h"




///This currently tells you the winding number of all the species. It
///might make sense to fix this so it tells only the winding number of
///a specific species
void 
WindingNumberClass::Accumulate()
{
  SamplesInBlock++;
  // Move the join to the end so we don't have to worry about
  // permutation
  PathData.MoveJoin(PathData.NumTimeSlices()-1);

  dVec totalDisp;

  totalDisp = 0.0;

  int numLinks=PathData.Path.NumTimeSlices()-1;
  //  for (int si=0; si<SpeciesList.size(); si++) {
  //    SpeciesClass &species = PathData.Path.Species(SpeciesList(si));
  //    int first = species.FirstPtcl;
  //    int last  = species.LastPtcl;
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
    for (int slice=0;slice<numLinks;slice++) {
       dVec disp;
       disp = PathData.Path.Velocity(slice,slice+1,ptcl);
       totalDisp += disp;
    }
  }

  for (int i=0; i<NDIM; i++) {
    totalDisp[i] *= PathData.Path.GetBoxInv()[i];
  }

  WNVec.push_back(totalDisp);
}

void 
WindingNumberClass::Read(IOSectionClass& in)
{
  ObservableClass::Read(in);
  Array<string,1> speciesStrings;
  assert(in.ReadVar("SpeciesList", speciesStrings));
  SpeciesList.resize(speciesStrings.size());
  for (int i=0; i<speciesStrings.size(); i++) {
    SpeciesList(i) = PathData.Path.SpeciesNum(speciesStrings(i));
    if (SpeciesList(i) == -1) {
      perr << "Unrecognized species name " << speciesStrings(i) << "\n";
      abort();
    }
  }
  WN2Array.resize(NDIM);
  WN2ArrayLowVariance.resize(NDIM);
}


void 
WindingNumberClass::CalcWN2()
{
  if (SamplesInBlock == 0) {
    cerr << "WindingNumberClass::CalcWN2 called before there is any data available.\n";
    abort();
  }

  // Better variance version
  WN2ArrayLowVariance = 0.0;
  for (int dim = 0; dim < NDIM; dim++){
    double W2avg = 0.0; 
    double Wavg = 0.0;
    for (int i = 0; i < WNVec.size(); i++){
      W2avg += WNVec[i][dim] * WNVec[i][dim];
      Wavg += WNVec[i][dim];
    }
    double Wavg2 = Wavg*Wavg;
    double W2avgsum;
    double Wavgsum;
    double Wavg2sum;
    W2avgsum = PathData.Path.Communicator.Sum(W2avg);
    Wavgsum = PathData.Path.Communicator.Sum(Wavg);
    Wavg2sum = PathData.Path.Communicator.Sum(Wavg2);
    if (PathData.Path.Communicator.MyProc() == 0){
      double N = (double)(PathData.Path.Communicator.NumProcs());
      double k = (double)(WNVec.size());
      WN2ArrayLowVariance(dim) = W2avgsum/k + ((Wavgsum*Wavgsum)-Wavg2sum)/(k*k);
    }
  }
  //end better variance version

  Array<dVec,1> sendVec(WNVec.size()), recvVec(WNVec.size());
  recvVec = 0.0;
  for (int i = 0; i < WNVec.size(); i++)
    sendVec(i) = WNVec[i];
  PathData.Path.Communicator.Sum(sendVec, recvVec);

  if (PathData.Path.Communicator.MyProc() == 0) {
    WN2Array = 0.0;
    for (int i = 0; i < recvVec.size(); i++) {
      for (int j = 0; j < NDIM; j++) {
        WN2Array(j) += recvVec(i)[j] * recvVec(i)[j];
      }
    }
    WN2Array /= (double)SamplesInBlock;
  }
  WNVec.clear();
  SamplesInBlock = 0;
}

void
WindingNumberClass::WriteBlock()
{
  CalcWN2();
  // Only processor 0 writes.
  if (PathData.Path.Communicator.MyProc()==0) {
    if (FirstTime) {
      FirstTime = false;
      WriteInfo();
      IOSection.WriteVar("Type",string("Vector"));
    }
    WNVar.Write(WN2Array);
    WNVarLowVariance.Write(WN2ArrayLowVariance);
  }
}


//   Array<double,1> dummy(3);
//   double beta=PathData.Path.tau*(PathData.Path.NumTimeSlices()-1);
//   double norm=(double)NumSamples*(2*PathData.Path.Species(0).lambda*
// 				  beta*PathData.Path.NumParticles());
//   for (int dim=0;dim<NDIM;dim++)
//     dummy(dim)=TotalW2[dim]/norm;
//   WNVar.Write(dummy);
//   NumSamples = 0;
//   TotalW2=0.0;

