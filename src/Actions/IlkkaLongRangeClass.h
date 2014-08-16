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

#ifndef ILKKA_LONG_RANGE_CLASS_H
#define ILKKA_LONG_RANGE_CLASS_H

#include "ActionBase.h"

/// The LongRangeClass is an action class responsible for the long
/// wavelength components of the action that are summed in k-space.
/// This class reads in and uses an optimized breakup from a file that
/// Ilkka supplies. 
class IlkkaLongRangeClass : public ActionBaseClass
{
protected:
  Array<PairActionFitClass*,2> &PairMatrix;
  Array<PairActionFitClass*,1> &PairArray;
  Array<int,2> &PairIndex;

  inline double mag2 (const complex<double> &z)
  {
    return (z.real()*z.real() + z.imag()*z.imag());
  }

  inline double mag2 (const complex<double> &z1, const complex<double> &z2)
  {
    return (z1.real()*z2.real() + z1.imag()*z2.imag());
  }

public:
  Array<double,1> specNum1, specNum2;
  Array<double,1> uk0, duk0, vk0, ur0, dur0, vr0;
  Array<double,2> uk, duk, Vlong_k;

  void Read (IOSectionClass &in);
  void Build_MultipleSpecies();
  double SingleAction (int slice1, int slice2, const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2,  int level);
  double V(int slice1,int slice2,int level);
  string GetName();
  bool fequals(double a,double b, double tol);
  bool vecEquals(dVec &a, dVec &b,double tol);
  void WriteInfo(IOSectionClass &out);
  IlkkaLongRangeClass::IlkkaLongRangeClass(PathDataClass &pathData, Array<PairActionFitClass*,2> &pairMatrix, Array<PairActionFitClass*,1> &pairArray, Array<int,2> &pairIndex)
    : ActionBaseClass (pathData), PairMatrix(pairMatrix), PairArray(pairArray), PairIndex(pairIndex)
  {}
};

#endif
