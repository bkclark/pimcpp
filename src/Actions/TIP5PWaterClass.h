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

#ifndef TIP5PWATER_CLASS_H
#define TIP5PWATER_CLASS_H

#include "ActionBase.h"

/// The KineticClass calculates the kinetic part of the action.  This
/// is the "spring term", of the form
/// \f$ K \equiv \left(4\pi\lambda\tau\right)^{-\frac{ND}{2}} 
/// \exp\left[-\frac{(R-R')^2}{4\lambda\tau}\right] \f$
class TIP5PWaterClass : public ActionBaseClass
{
public:
  void Read (IOSectionClass &in);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  string GetName();
  double OOSeparation (int slice,int ptcl1,int ptcl2);
  double RotationalKinetic(int startSlice, int endSlice, const Array<int,1> &activeParticles,int level);
  double RotationalEnergy(int startSlice, int endSlice, int level);
  void GetAngles(dVec disp, double &theta, double &phi);
  double CalcEnergy(double reftheta,double dtheta, double dphi);
  double SineOfPolar(dVec coords);
  dVec COMVelocity (int slice1,int slice2,int ptcl);
  dVec COMCoords (int slice, int ptcl);
  dVec Displacement(int slice1, int slice2, int ptcl1, int ptcl2);
  int FindCOM(int ptcl);
  double ProtonKineticAction (int slice1, int slice2, const Array<int,1> &changedParticles, int level);
  double ProtonKineticEnergy (int slice1, int slice2, int level);
  double SecondProtonKineticAction(int startSlice, int endSlice, const Array<int,1> &activeParticles,int level);
  double SecondProtonKineticEnergy(int startSlice, int endSlice, int level);
  double dotprod(dVec vec1, dVec vec2, double mag);
  int FindOtherProton(int ptcl);
  double NewRotKinAction(int startSlice, int endSlice, const Array<int,1> &activeParticles, int level);
  double NewRotKinEnergy(int startSlice, int endSlice, int level);
  double FixedAxisAction(int startSlice, int endSlice, const Array<int,1> &activeParticles, int level);
  double FixedAxisEnergy(int startSlice, int endSlice, int level);
  dVec CrossProd(dVec v1, dVec v2);
  double Mag(dVec v);
  double CalcCutoff(int ptcl1, int ptcl2, int slice,double Rcmag);
  double GetAngle(dVec v1, dVec v2);
  dVec Normalize(dVec v);
  dVec Scale(dVec v, double scale);
  dVec Strip(dVec R, dVec u);
  dVec GetBisector(dVec v1, dVec v2);
  double CalcPsi(double theta);
  bool dVecsEqual(dVec u,dVec v);
  TIP5PWaterClass (PathDataClass &pathData);
};

const double O_H_moment_arm = 1.0;//st2
//const double O_H_moment_arm = 0.9572;//tip5p
const double lambda_p = 0.047848;
//const double HOH_angle = 104.52;
//const double HOH_angle = 1.8242; //tip5p
const double HOH_angle = 1.9106; //st2
const double HOH_half_angle = HOH_angle/2;
const double hbar = 9.5425095602294465*pow(10.0,-14);


#endif
