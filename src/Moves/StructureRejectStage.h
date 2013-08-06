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

#ifndef STRUCTURE_REJECT_STAGE_CLASS_H
#define STRUCTURE_REJECT_STAGE_CLASS_H

#include "MultiStage.h"
#include "../Observables/StructureFactor.h"

class StructureRejectStageClass : public LocalStageClass
{
private:
  int Species1;
  int Species2;
  double MaxValue;
  StructureFactorClass StructureFactor;
public:
  void Read(IOSectionClass &in);
  double Sample (int &slice1,int &slice2,
		 Array<int,1> &changedParticles); 
  bool Attempt(int &slice1, int &slice2,
	       Array<int,1> &activeParticles,
	       double &prevActionChange);

  StructureRejectStageClass (PathDataClass &pathData, IOSectionClass &in,
			     IOSectionClass &out) 
    : LocalStageClass(pathData,out),
    StructureFactor(pathData,in)
  {
    // do nothing for now
  }

};

#endif
