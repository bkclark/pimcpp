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

#include "Coupling.h"

////////////////////////////////////////
///Pair Correlation Class           ///
///////////////////////////////////////

void CouplingClass::Read(IOSectionClass& in)
{
  
  ObservableClass::Read(in);
  TotalCounts=0;
  Coupling.resize(110);
  Coupling=0;
  /// Now write the one-time output variables
  if (PathData.Path.Communicator.MyProc()==0)
    WriteInfo();
}



void CouplingClass::WriteInfo()
{
  ObservableClass::WriteInfo();
  IOSection.WriteVar("ylabel", "g(r)");
  IOSection.WriteVar("Cumulative", false);
}



void CouplingClass::WriteBlock()
{
  PathClass &Path= PathData.Path;
  Array<int,1> CouplingSum(Coupling.size());
  double norm=0.0;
  norm=TotalCounts;
  Path.Communicator.Sum(Coupling, CouplingSum);
  Array<double,1> couplingArray(CouplingSum.size());
  for (int i=0; i<couplingArray.size(); i++){
    couplingArray(i) = (double) CouplingSum(i) / (norm);
  }
  couplingVar.Write(couplingArray);
  couplingVar.Flush();
  Coupling = 0;
  TotalCounts = 0;
}

/// Fix me to accumulate data only between the two species I'm
/// interested in.
void CouplingClass::Accumulate()
{

  TotalCounts++;
  int couplingBin=(int)((floor)((PathData.Path.ExistsCoupling*100)+0.2));
  assert(couplingBin>=0);
  assert(couplingBin<Coupling.size());
  Coupling(couplingBin)++;
}


void CouplingClass::Initialize()
{
  TotalCounts = 0;
  TimesCalled=0;
}


