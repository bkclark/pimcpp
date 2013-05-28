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

#ifndef BLEND_ACTIONS_H
#define BLEND_ACTIONS_H

#include "ActionBase.h"

class BlendActionsClass : public ActionBaseClass
{
  // store two actions to blend in the form:
  // S = (1-lambda_blend) * S0 + lambda_blend * S1
  Array<ActionBaseClass*,1> S0;
  Array<ActionBaseClass*,1> S1;
  // blending parameter (must be betwen 0 and 1)
  double lambda_blend;

  public:
  double SingleAction(int slice1, int slice2,
			      const Array<int,1> &activeParticles,
			      int level);

  double d_dBeta (int slice1, int slice2,
			  int level);

  void Read (IOSectionClass &in);
  string GetName();
  BlendActionsClass(PathDataClass &pathData);				   
};

#endif
