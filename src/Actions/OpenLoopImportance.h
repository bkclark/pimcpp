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

#ifndef OPENLOOP_IMPORTANCE_H
#define OPENLOOP_IMPORTANCE_H

#include "ActionBase.h"

typedef enum {NOIMP,DISTIMP,DISPXIMP,CONSTSHIFT,TUNEDFUNCTION,RETUNEDFUNCTION, POLYNOMIAL,EXPONENTIAL} SampleChoice;

typedef enum {MIN_IMAGE_DISP,MIN_IMAGE_DIST,DISP,DIST} Xvalue;

class OpenLoopImportanceClass : public ActionBaseClass
{
public:
  Array<double,1> Polynom;
  double a;
  double alpha;
  double s;
  void Read (IOSectionClass &in);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  string GetName();
  OpenLoopImportanceClass (PathDataClass &pathData);
  SampleChoice ImpChoice;
  Xvalue Xis;
  double Shift;
};


#endif
