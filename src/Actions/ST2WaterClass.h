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

#ifndef ST2WATER_CLASS_H
#define ST2WATER_CLASS_H

#include "ActionBase.h"

/// The KineticClass calculates the kinetic part of the action.  This
/// is the "spring term", of the form
/// \f$ K \equiv \left(4\pi\lambda\tau\right)^{-\frac{ND}{2}} 
/// \exp\left[-\frac{(R-R')^2}{4\lambda\tau}\right] \f$
class ST2WaterClass : public ActionBaseClass
{
public:
  void Read (IOSectionClass &in);
	double Action (int startSlice, int endSlice, const Array<int,1> &activeParticles, int level);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  string GetName();
  double EField (Array<int,1> &activeMol, int startSlice, int endSlice,  int level);
  void EFieldVec (int molIndex, dVec &Efield, double &Emag, int slice, int level);
  void BackSub(Array<double,2> &A,Array<double,1> &X,int N, int C);
  void RowScale(Array<double,2> &A,int n,double scale,int C);
  void Diff(Array<double,2> &A,int n,int m, double s, int C);
  double OOSeparation (int slice,int ptcl1,int ptcl2);
  double S(double r);


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
  double GetAngle(dVec v1, dVec v2);
  dVec Normalize(dVec v);
  dVec Scale(dVec v, double scale);
  dVec GetBisector(dVec v1, dVec v2);
  double CalcPsi(double theta);
  ST2WaterClass (PathDataClass &pathData);
};

const double RL = 2.0160;
const double RU = 3.1287;

#endif
