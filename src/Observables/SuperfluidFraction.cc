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

#include "SuperfluidFraction.h"
#include <Common/MPI/Communication.h>

void 
SuperfluidFractionClass::Read(IOSectionClass& IO)
{
  WindingNumberClass::Read(IO);
  assert(SpeciesList.size() == 1);
}

void
SuperfluidFractionClass::WriteBlock()
{
  int species = SpeciesList(0);
  double beta = PathData.Path.tau * PathData.Path.TotalNumSlices;
  int numParticles = PathData.Path.Species(species).LastPtcl - PathData.Path.Species(species).FirstPtcl + 1;
  double factor = 1.0/(2.0 * PathData.Path.Species(species).lambda * beta * numParticles);

  CalcWN2();
  // Only processor 0 writes.
  if (PathData.Path.Communicator.MyProc() == 0) {
    if (FirstTime) {
      FirstTime = false;
      WriteInfo();
      IOSection.WriteVar("Type",string("Vector"));
    }
    for (int dim = 0; dim < NDIM; dim++)
      WN2Array(dim) = WN2Array(dim) * factor * PathData.Path.GetBox()[dim] * PathData.Path.GetBox()[dim];
    SFVar.Write(WN2Array);
  }
}

