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

#ifndef KINETIC_INTRA_VIBRATION_CLASS_H
#define KINETIC_INTRA_VIBRATION_CLASS_H

#include "ActionBase.h"

/// The KineticVibrationClass calculates the kinetic part of the action.  This
/// is the "spring term", of the form
/// \f$ K \equiv \left(4\pi\lambda\tau\right)^{-\frac{ND}{2}} 
/// \exp\left[-\frac{(R-R')^2}{4\lambda\tau}\right] \f$
class KineticVibrationClass : public ActionBaseClass
{
  double lambda;
public:
  int NumImages;

  void Read (IOSectionClass &in);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double SingleActionForcedTau (int slice1, int slice2,
				const Array<int,1> &changedParticles, 
				int level,
				double forcedTau);
  double d_dBeta (int slice1, int slice2, int level);
  double d_dBetaForcedTau (int slice1, int slice2,
			   int level,
			   double forcedTau);
  

  string GetName();
  inline void SetNumImages (int num) { NumImages = num; }
  KineticVibrationClass (PathDataClass &pathData);
};

class KineticVibrationEigenFunctionClass : public ActionBaseClass
{
  double lambda, R_OH_0, R_HH_0, omega_OH, omega_HH;
  int maxH;
public:
  void Read (IOSectionClass &in);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double SingleActionForcedTau (int slice1, int slice2,
				const Array<int,1> &changedParticles, 
				int level,
				double forcedTau);
  double d_dBeta (int slice1, int slice2, int level);
  double d_dBetaForcedTau (int slice1, int slice2,
			   int level,
			   double forcedTau);
  
  string GetName();
  KineticVibrationEigenFunctionClass (PathDataClass &pathData);
};

#endif
