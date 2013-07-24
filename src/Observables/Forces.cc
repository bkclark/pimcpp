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

#include "Forces.h"

void
ForcesClass::SetSpecies(int speciesNum)
{
  TimesCalled = 0;
  SpeciesNum = speciesNum;
  SpeciesClass &species = PathData.Path.Species(speciesNum);
  Forces.resize(species.NumParticles);
  SumTmp.resize(species.NumParticles);
  dVec zero;
  for (int i=0; i<NDIM; i++)
    zero[i] = 0.0;
  Forces = zero;
  Ptcls.resize(species.NumParticles);
  for (int i=0; i<Ptcls.size(); i++)
    Ptcls(i) = i + species.FirstPtcl;
  ForcesArray.resize(Ptcls.size(),NDIM);
}


void
ForcesClass::Accumulate()
{
#if NDIM==3
//     Array<dVec,1> Fanalytic(Ptcls.size()), FFD(Ptcls.size());
//     dVec zero;
//     for (int i=0; i<NDIM; i++) zero[i] = 0.0;
//     Fanalytic = zero;
//     FFD = zero;
//     PathData.Actions.GetForces(Ptcls, Fanalytic);
//     PathData.Actions.GetForcesFD(Ptcls, FFD);
//     cerr << "Forces:\n";
//     for (int i=0; i<Ptcls.size(); i++) {
//       for (int dim=0; dim<NDIM; dim++)
// 	fprintf (stderr, "%12.6f ", Fanalytic(i)[dim]);
//       for (int dim=0; dim<NDIM; dim++)
// 	fprintf (stderr, "%12.6f ", FFD(i)[dim]);
//       fprintf (stderr, "\n");
//     }

  PathData.Actions.GetForces(Ptcls, Forces, Forces);
  Counts++;
#else
  cerr<<"2D not implemented in forces"<<endl;
  assert(NDIM==3);
#endif
}


void
ForcesClass::WriteBlock()
{
  double norm = 1.0/(double)Counts;
  /// Sum up over all processors to get the total for the path.
  PathData.Path.Communicator.AllSum(Forces,SumTmp);
  Forces = SumTmp;
  for (int i=0; i<Forces.size(); i++)
    for (int j=0; j<NDIM; j++) 
      ForcesArray(i,j) = norm * Forces(i)[j];
  ForcesVar.Write(ForcesArray);
  dVec zero;
  for (int i=0; i<NDIM; i++)
    zero[i] = 0.0;
  Forces = zero;
  Counts = 0;
}


void
ForcesClass::WriteInfo()
{
  IOSection.WriteVar("SpeciesNum", SpeciesNum);
  Array<string,1> actionTypes(2);
  actionTypes(0) = "ShortRange";
  actionTypes(1) = "LongRange";
  IOSection.WriteVar("ActionTypes", actionTypes);
}


void
ForcesClass::Read(IOSectionClass &in)
{
  string speciesString;
  ObservableClass::Read(in);
  assert (in.ReadVar("Species", speciesString));
  SpeciesNum = PathData.Path.SpeciesNum(speciesString);
  SetSpecies(SpeciesNum);
  WriteInfo();
}
