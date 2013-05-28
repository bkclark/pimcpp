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

#ifndef PA_ZERO_FIT_H
#define PA_ZERO_FIT_H
#include "PAFitBase.h"


class PAzeroFitClass : public PairActionFitClass
{
public:

  void ReadParams  (IOSectionClass &inSection) {};
  void WriteBetaIndependentInfo (IOSectionClass &outSection) {};
  /// Returns weighter RMS error
  //  void Error (Rho &rho, double &Uerror, double &dUerror) {};
  //  void DoFit (Rho &rho) {};
  void WriteFit(IOSectionClass &outSection) {};
  void ReadInput (IOSectionClass &inSection) {};
  
  bool Read  (IOSectionClass &inSection, double lowestBeta,
	      int NumBetas);
  double U(double q, double z, double s2, int level);
  double dU(double q, double z, double s2, int level);
  void Derivs (double q, double z, double s2, int level,
	       double &d_dq, double &d_dz, double &d_ds);
  double V(double r);
  bool IsLongRange();
  PAzeroFitClass()
  { 

  }
};

#endif
