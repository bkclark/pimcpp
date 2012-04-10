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

#ifndef SHORT_RANGE_PRIMITIVE_CLASS_H
#define SHORT_RANGE_PRIMITIVE_CLASS_H

#include "ActionBase.h"
#include "../PairAction/PAFit.h"

/// The ShortRangePrimitiveClass is an action class which stores the part of
/// the potential pair action that is summed in real space.  If the
/// potential is short range, it contains the whole potential action.
/// This action, in general, contains off diagaonal contributions.
class ShortRangePrimitiveClass : public ActionBaseClass
{
protected:
  Array<PairActionFitClass*,2> &PairMatrix;
  Array<bool,1> DoPtcl;
public:
  void Read (IOSectionClass &in);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  string GetName();
  ShortRangePrimitiveClass (PathDataClass &pathData,
		   Array<PairActionFitClass*, 2> &pairMatrix);
};

#endif
