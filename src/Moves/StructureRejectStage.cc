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

#include <Common/MPI/Communication.h>
#include "StructureRejectStage.h"


double StructureRejectStageClass::Sample(int &slice1,int &slice2,
				      Array<int,1> &changedParticles)
{
  //do nothing for now
}




void StructureRejectStageClass::Read (IOSectionClass &in)
{
  assert(in.ReadVar("MaxValue",MaxValue));
  StructureFactor.Read(in);
}

bool StructureRejectStageClass::Attempt (int &slice1, int &slice2, 
			       Array<int,1> &activeParticles, 
			       double &prevActionChange)
{
  bool toReturn=false;
  //  double tcounts=(double)StructureFactor.TotalCounts;
  //  double nslices=(double)PathData.Path.NumTimeSlices();
  double norm=PathData.Path.NumParticles();
  StructureFactor.Calculate();
  for (int counter=0;counter<StructureFactor.Sk.size();counter++){
    cerr<<StructureFactor.Sk(counter)/norm<<endl;
    if (StructureFactor.Sk(counter)/norm>4){
      toReturn=true;
    }
  }
  StructureFactor.Clear();
  return toReturn;
}
