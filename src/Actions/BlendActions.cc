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

#include "BlendActions.h"
#include "../PathDataClass.h"

BlendActionsClass::BlendActionsClass(PathDataClass &pathData) : 
  ActionBaseClass(pathData)
{
}

void BlendActionsClass::Read (IOSectionClass &in)
{
  Array<string,1> sname;
  assert(in.ReadVar("FirstAction",sname));
  S0.resize(sname.size());
  for(int s=0; s<sname.size(); s++)
    S0(s) = PathData.Actions.GetAction(sname(s));
  assert(in.ReadVar("SecondAction",sname));
  S1.resize(sname.size());
  for(int s=0; s<sname.size(); s++)
    S1(s) = PathData.Actions.GetAction(sname(s));
  assert(in.ReadVar("Lambda",lambda_blend));
}


double BlendActionsClass::SingleAction(int slice1, int slice2, const Array<int,1> &activeParticles, int level)
{
  double Ublend = 0.0;
  for(int s=0; s<S0.size(); s++)
    Ublend += (1 - lambda_blend) * S0(s)->SingleAction(slice1, slice2, activeParticles, level);
  for(int s=0; s<S1.size(); s++)
    Ublend += lambda_blend * S1(s)->SingleAction(slice1, slice2, activeParticles, level);
  return Ublend;
}

double BlendActionsClass::d_dBeta (int slice1, int slice2, int level)
{
  double Ublend = 0.0;
  for(int s=0; s<S0.size(); s++)
    Ublend += (1 - lambda_blend) * S0(s)->d_dBeta(slice1, slice2, level);
  for(int s=0; s<S1.size(); s++)
    Ublend += lambda_blend * S1(s)->d_dBeta(slice1, slice2, level);
  return Ublend;
}

string BlendActionsClass::GetName()
{
  return("BlendActionsClass");
}
