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

#ifndef HARMONIC_POTENTIAL_CLASS_H
#define HARMONIC_POTENTIAL_CLASS_H

#include "ActionBase.h"

/// The HarmonicPotentialClass is an action class which stores the part of
/// the potential pair action that is summed in real space.  If the
/// potential is short range, it contains the whole potential action.
/// This action, in general, contains off diagaonal contributions.
class HarmonicPotentialClass : public ActionBaseClass
{
protected:
  int TotalTime;
  double omega;
public:
  void Read (IOSectionClass &in);
  double dUdR(int slice, int ptcl1, int level);
  double d2UdR2(int slice, int ptcl1, int ptcl2, int level);
  double SingleAction (int slice1, int slice2, const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  string GetName();
  HarmonicPotentialClass (PathDataClass &pathData);
};

#endif
