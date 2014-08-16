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

#ifndef Al_EMBEDDED_ATOM_POTENTIAL_CLASS_H
#define Al_EMBEDDED_ATOM_POTENTIAL_CLASS_H

#include "ActionBase.h"

/// Implementation of an Embeded Atom potential
/// (EAM) for Al.
/// See Mei and Davenport, PRB 46, 21 (1992)

class EAMPotentialClass : public ActionBaseClass
{
  protected:
  double Ec, phi0, r0, alpha, beta, gamma, delta;
  double rn, rc;
  double conversion;
  Array<double, 1> c;
  Array<double, 1> s;
  int mMax, lMax;
  int Level, ki;
  double aob;

  Array<double, 2> rho;

  public:
  void Read (IOSectionClass &in);
  double ComputeEnergy (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);

  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);

  double d_dBeta (int slice1, int slice2,  int level);

  string GetName();
  double F(double r);
  double f(double r);
  double q(double r);
  double phi(double r);
  void UpdateRho(int slice1, int slice2, const Array<int,1>& activeP);

  EAMPotentialClass (PathDataClass &pathData);

};

#endif
