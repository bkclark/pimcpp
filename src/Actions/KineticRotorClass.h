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

#ifndef KINETIC_ROTOR_CLASS_H
#define KINETIC_ROTOR_CLASS_H

#include "ActionBase.h"
#include <Common/Splines/Grid.h>

#include <Common/Splines/BicubicSpline.h>

// temporary to calculate explicitly...
// but very slow... only for testing!
//#include "../Rotor/RotorContainer.h"
//#include "../Rotor/WignerFunction.h"
//#include "../Rotor/RotorDensityMatrix.h"
////#include "../Rotor/RhoFixedAxis.h"
//#include "../Rotor/legendre/legendre.h"

class RotorActionBaseClass : public ActionBaseClass
{
  public:
  void getAngles(int slice1, int mol1, int slice2, int mol2, double& theta, double& phi, double& chi);

  RotorActionBaseClass(PathDataClass &pathData) :
    ActionBaseClass(pathData)
    {
    }
};

bool isEqualOpp(dVec u, dVec v);
bool isEqual(dVec u, dVec v);

/// The KineticRotorClass calculates the 
/// rotational kinetic part of the action 
/// for a rigid asymmetric rotor.
/// It was designed for water simulations.
/// This is the "spring term", of the form
/// \f$ K \equiv \left(4\pi\lambda\tau\right)^{-\frac{ND}{2}} 
/// \exp\left[-\frac{(R-R')^2}{4\lambda\tau}\right] \f$

class KineticRotorClass : public RotorActionBaseClass
{
  //moments of inertia
  double Ia, Ib, Ic;
  double A, B, C;

  Grid* ThetaGrid;
  Grid* PhiGrid;
  //Grid* ChiGrid;
  //TricubicSpline spline;
  BicubicSpline spline;

  Grid* ThetaEnergyGrid;
  Grid* PhiEnergyGrid;
  //Grid* ChiEnergyGrid;
  //TricubicSpline EnergySpline;
  BicubicSpline EnergySpline;

  // explicit calculatoin for testing
  //RotorRhoClass* rho; 

  // some things for debugging
  //ofstream err, tolrho, overX, overS;
  ofstream tout, pout, chiout;
  Array<double,1> thist, phist, chist;
  double DT, DP;
  int times;
  double maxRhoX, maxRhoS;
  double TotalZ, TotalE, TotalESq;
  double THETA, PHI, CHI;
  bool doESq;
  ofstream esqOut;
public:
  int NumImages;

  void Read (IOSectionClass &in);
  void ReadGridPoints(IOSectionClass& in);
  void ReadEnergyGridPoints(IOSectionClass& in);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &changedParticles, int level);
//  double SingleActionForcedTau (int slice1, int slice2,
//				const Array<int,1> &changedParticles, 
//				int level,
//				double forcedTau);
  double d_dBeta (int slice1, int slice2, int level);
//  double d_dBetaForcedTau (int slice1, int slice2,
//			   int level,
//			   double forcedTau);
  
  string GetName();
  inline void SetNumImages (int num) { NumImages = num; }
  KineticRotorClass (PathDataClass &pathData);
};


class FixedAxisRotorClass : public RotorActionBaseClass
{
  double Ia, Ib, Ic;
  double C;

  public:

  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &changedParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  double CalcRho(double phi, double theta, double psi);
  double CalcFAEnergy(double phi, double theta, double psi);
  double AltFixedAxis (int slice, int mol, double& phi, double& theta, double& psi);
  void Read (IOSectionClass &in);
  string GetName();
  
  FixedAxisRotorClass(PathDataClass &pathData);
};

#endif
