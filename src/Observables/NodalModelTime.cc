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

#include "NodalModelTime.h"


/// Writes the data relevant for this classes output including
/// its name, axis to plot, etc.
void NodalModelTimeClass::WriteInfo()
{
  ObservableClass::WriteInfo();
}


/// Writes a block of histogram data.
void NodalModelTimeClass::WriteBlock()
{
  PathClass &Path = PathData.Path;
  int myProc = Path.Communicator.MyProc();

  if (myProc == 0) {
    counts = counts/NumSamples;
    NodalModelTimeVar.Write(counts);
    NodalModelTimeVar.Flush();
  }
  counts = 0;

  NumSamples = 0;
}


void NodalModelTimeClass::Read(IOSectionClass &in)
{
  ObservableClass::Read(in);

  assert(in.ReadVar("Species",Species));
  int SpeciesNum = PathData.Path.SpeciesNum(Species);
  int NumModels = PathData.Actions.NodalActions(SpeciesNum) -> GetNumModels();

  counts.resize(NumModels);
  counts = 0;

  /// Now write the one-time output variables
  if (PathData.Path.Communicator.MyProc() == 0) {
    WriteInfo();
  }

}

// Check out MetaMoves.cc! Shift was changed to 0 b/c for some reason it's messing this procedure up!
void NodalModelTimeClass::Accumulate()
{
  NumSamples++;
  int myProc = Path.Communicator.MyProc();
  int SpeciesNum = PathData.Path.SpeciesNum(Species);
  int model = PathData.Actions.NodalActions(SpeciesNum) -> GetModel();
  counts(model) += 1;
}
